import pandas as pd
import glob
import pysam
import os
import numpy as np
from collections import defaultdict

# def pairs_generator(file, refname=None):
    
#     read_dict = defaultdict(lambda: [None, None])
#     file.reset()
#     for read in file.fetch(refname):
#         qname = read.query_name
#         if qname in read_dict:
#             if read.is_read1:
#                 yield qname, read, read_dict[qname][1]
#             else:
#                 yield qname, read_dict[qname][0], read
#             del read_dict[qname]
#         else:
#             if read.is_read1:
#                 read_dict[qname][0] = read
#             else:
#                 read_dict[qname][1] = read


def get_pairs_from_pairend_alignment(directory, refname=None):
    """
    Retrieve pairs from paired-end alignment file and output pairs as a dictionary.
    """
    
    # check for alignment file and index file
    try:
        filename = glob.glob(os.path.join(directory, '*final.paired.aln_all.bam'))[0]
    except:
        print('Paired end alignment file not found!\nExiting from get_pairs_from_pairend_alignment')
        exit()
    if not os.path.isfile(filename+'.bai'):
        pysam.samtools.index(filename)
    else:
        if os.path.getctime(filename+'.bai') < os.path.getctime(filename):
            pysam.samtools.index(filename)
    
    output_dict = dict()
    read_dict = defaultdict(lambda: [None, None])
    
    alignmentfile = pysam.AlignmentFile(filename, 'rb')
    alignmentfile.reset()
    for read in alignmentfile.fetch(refname, until_eof=True):
        qname = read.query_name
        if qname in read_dict:
            # sanity check: no qname should have been in output_dict already
            if qname in output_dict:
                print('Duplicate read warning: {}, skip and continue'.format(qname))
                continue
            if read.is_unmapped or read.get_tag('AS') < 100:
                if read.is_read1:
                    r2 = read_dict[qname][1]
                    if r2 is 'unmapped':
                        output_dict[qname] = ['unmapped', 'unmapped', np.nan, np.nan,
                                              r2, r2, np.nan, np.nan,
                                              False]
                    else:
                        output_dict[qname] = ['unmapped', 'unmapped', np.nan, np.nan,
                                              alignmentfile.get_reference_name(r2.reference_id),
                                              r2.reference_start,
                                              r2.mapping_quality,
                                              r2.is_reverse,
                                              False]
                    
                else:
                    r1 = read_dict[qname][0]
                    if r1 is 'unmapped':
                        output_dict[qname] = [r1, r1, np.nan, np.nan,
                                              'unmapped', 'unmapped', np.nan, np.nan,
                                              False]
                    else:
                        output_dict[qname] = [alignmentfile.get_reference_name(r1.reference_id),
                                              r1.reference_start,
                                              r1.mapping_quality,
                                              r1.is_reverse,
                                              'unmapped', 'unmapped', np.nan, np.nan,
                                              False]
            else:
                if read.is_read1:
                    r2 = read_dict[qname][1]
                    if r2 is 'unmapped':
                        output_dict[qname] = [alignmentfile.get_reference_name(read.reference_id),
                                              read.reference_start,
                                              read.mapping_quality,
                                              read.is_reverse,
                                              r2, r2, np.nan, np.nan,
                                              False]
                    else:
                        output_dict[qname] = [alignmentfile.get_reference_name(read.reference_id),
                                              read.reference_start,
                                              read.mapping_quality,
                                              read.is_reverse,
                                              alignmentfile.get_reference_name(r2.reference_id),
                                              r2.reference_start,
                                              r2.mapping_quality,
                                              r2.is_reverse,
                                              read.is_proper_pair]
                    
                else:
                    r1 = read_dict[qname][0]
                    if r1 is 'unmapped':
                        output_dict[qname] = [r1, r1, np.nan, np.nan,
                                              alignmentfile.get_reference_name(read.reference_id),
                                              read.reference_start,
                                              read.mapping_quality,
                                              read.is_reverse,
                                              False]
                    else:
                        output_dict[qname] = [alignmentfile.get_reference_name(r1.reference_id),
                                              r1.reference_start,
                                              r1.mapping_quality,
                                              r1.is_reverse,
                                              alignmentfile.get_reference_name(read.reference_id),
                                              read.reference_start,
                                              read.mapping_quality,
                                              read.is_reverse,
                                              r1.is_proper_pair]
            # delete entry from dictionary as soon as its mate is found
            del read_dict[qname]
        else:
            if read.is_unmapped or read.get_tag('AS') < 100:
                if read.is_read1:
                    read_dict[qname][0] = 'unmapped'
                else:
                    read_dict[qname][1] = 'unmapped'
            else:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read

    # sanity check: read_dict should be empty if all single-end reads found their mates
    if bool(read_dict):
        print('Warning: pairing error!')
        for qname in read_dict:
            print('Mate is not found for '+qname)
        print('Skip and continue')
        
    alignmentfile.close()
    print('Done processing paired-end alignment file!')
    
    return output_dict


def get_pairs_from_singleend_alignment(directory, refname=None):
    """
    Retrieve pairs from single-end alignment file and improperly paired reads with high score.
    Output pairs as a dictionary.
    """
    
   # check for alignment + original files and index files
    try:
        files = [glob.glob(os.path.join(directory, '*final.single.aln_all.bam'))[0],
                 glob.glob(os.path.join(directory, '*final.improper_highq.bam'))[0],
                 glob.glob(os.path.join(directory, '*final.map_unmap.bam'))[0]]

    except:
        print('One or more missing files!\nExiting from get_pairs_from_singleend_alignment')
        exit()
    for file in files:
        if not os.path.isfile(file+'.bai'):
            pysam.samtools.index(file)
        else:
            if os.path.getctime(file+'.bai') < os.path.getctime(file):
                pysam.samtools.index(file)
    
    output_dict = dict()
    read_dict = defaultdict(lambda: [None, None])
    
    alignmentfile = pysam.AlignmentFile(files[0], 'rb')
    duplicates = 0
    for read in alignmentfile.fetch(refname, until_eof=True):
        qname = read.query_name
        # sanity check: no qname should have been in output_dict already
        if qname in read_dict:
            print('Duplicate read warning! {} is a duplicate'.format(qname))
            duplicates += 1
            continue
        # single-end reads don't have read1/read2 tags
        if read.is_unmapped or read.get_tag('AS') < 100:
            read_dict[qname][0] = 'unmapped'
        else:
            read_dict[qname][0] = read
    if duplicates > 0:
        print('{} duplicate reads found in single-end alignment file. Skip and continue'.format(duplicates))

    print('{} reads retrieved from single-end realignment'.format(len(read_dict)))
    originalfiles = []
    for file in files[1:]:
        originalfile = pysam.AlignmentFile(file, 'rb')
        originalfiles.append(originalfile)
        for read in originalfile.fetch(refname, until_eof=True):
            qname = read.query_name
            if qname in read_dict:
                # qname is first added by traversing alignmentfile
                # retrieve reference for first read from alignmentfile
                # retrieve reference for second read from originalfile
                if read.is_read1 and read_dict[qname][0] is not None:
                    r2 = read_dict[qname][0]
                    if read.is_unmapped:
                        if r2 is 'unmapped':
                            output_dict[qname] = ['unmapped', 'unmapped', np.nan, np.nan,
                                                  r2, r2, np.nan, np.nan,
                                                  False]
                        else:
                            output_dict[qname] = ['unmapped', 'unmapped', np.nan, np.nan,
                                                  alignmentfile.get_reference_name(r2.reference_id),
                                                  r2.reference_start,
                                                  r2.mapping_quality,
                                                  r2.is_reverse,
                                                  False]
                    else:
                        if r2 is 'unmapped':
                            output_dict[qname] = [originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  r2, r2, np.nan, np.nan,
                                                  False]
                        else:
                            output_dict[qname] = [originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  alignmentfile.get_reference_name(r2.reference_id),
                                                  r2.reference_start,
                                                  r2.mapping_quality,
                                                  r2.is_reverse,
                                                  read.is_proper_pair]
                    # delete entry from dictionary as soon as its mate is found
                    del read_dict[qname]
                elif read.is_read2 and read_dict[qname][0] is not None:
                    r1 = read_dict[qname][0]
                    if read.is_unmapped:
                        if r1 is 'unmapped':
                            output_dict[qname] = [r1, r1, np.nan, np.nan,
                                                  'unmapped', 'unmapped', np.nan, np.nan,
                                                  False]
                        else:
                            output_dict[qname] = [alignmentfile.get_reference_name(r1.reference_id), 
                                                  r1.reference_start,
                                                  r1.mapping_quality,
                                                  r1.is_reverse,
                                                  'unmapped', 'unmapped', np.nan, np.nan,
                                                  False]
                    else:
                        if r1 is 'unmapped':
                            output_dict[qname] = [r1, r1, np.nan, np.nan,
                                                  originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  False]
                        else:
                            output_dict[qname] = [alignmentfile.get_reference_name(r1.reference_id), 
                                                  r1.reference_start,
                                                  r1.mapping_quality,
                                                  r1.is_reverse,
                                                  originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  r1.is_proper_pair]
                    # delete entry from dictionary as soon as its mate is found
                    del read_dict[qname]
                    
                # qname is first added by traversing originalfile
                # retrieve all references from originalfile
                elif read.is_read1 and read_dict[qname][1] is not None:
                    r2 = read_dict[qname][1]
                    if read.is_unmapped:
                        if r2 is 'unmapped':
                            output_dict[qname] = ['unmapped', 'unmapped', np.nan, np.nan,
                                                  r2, r2, np.nan, np.nan,
                                                  False]
                        else:
                            output_dict[qname] = ['unmapped', 'unmapped', np.nan, np.nan,
                                                  originalfile.get_reference_name(r2.reference_id),
                                                  r2.reference_start,
                                                  r2.mapping_quality,
                                                  r2.is_reverse,
                                                  False]
                    else:
                        if r2 is 'unmapped':
                            output_dict[qname] = [originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  r2, r2, np.nan, np.nan,
                                                  False]
                        else:
                            output_dict[qname] = [originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  originalfile.get_reference_name(r2.reference_id),
                                                  r2.reference_start,
                                                  r2.mapping_quality,
                                                  r2.is_reverse,
                                                  read.is_proper_pair]
                    # delete entry from dictionary as soon as its mate is found
                    del read_dict[qname]
                elif read.is_read2 and read_dict[qname][1] is not None:
                    r1 = read_dict[qname][1]
                    if read.is_unmapped:
                        if r1 is 'unmapped':
                            output_dict[qname] = [r1, r1, np.nan, np.nan,
                                                  'unmapped', 'unmapped', np.nan, np.nan,
                                                  False]
                        else:
                            output_dict[qname] = [originalfile.get_reference_name(r1.reference_id), 
                                                  r1.reference_start,
                                                  r1.mapping_quality,
                                                  r1.is_reverse,
                                                  'unmapped', 'unmapped', np.nan, np.nan,
                                                  False]
                    else:
                        if r1 is 'unmapped':
                            output_dict[qname] = [r1, r1, np.nan, np.nan,
                                                  originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  False]
                        else:
                            output_dict[qname] = [originalfile.get_reference_name(r1.reference_id), 
                                                  r1.reference_start,
                                                  r1.mapping_quality,
                                                  r1.is_reverse,
                                                  originalfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse,
                                                  r1.is_proper_pair]
                    # delete entry from dictionary as soon as its mate is found
                    del read_dict[qname]
                    
            # should only execute else when file == improper_highq.bam
            else:
                if 'map_unmap' in file:
                    print('Warning: {} not found in realignment'.format(qname))
                    
                # if file == improper_highq.bam, also retrieve reads not found in read_dict
                # as improperly paired high score reads
                else:
                    read_dict[qname][1] = read
        originalfile.close()        

    # sanity check: read_dict should be empty if all single-end reads found their mates
    if bool(read_dict):
        print('Pairing error! {} single-end reads do not have mates'.format(len(read_dict)))
        for qname in read_dict:
            print('Mate is not found for {}. Skip and continue'.format(qname))
    
    alignmentfile.close()
    print('Done processing single-end alignment file!')
    
    return output_dict


def generate_alignment_dict(directory, mapq_threshold=10, refname=None):
    """
    Retrieve pairs from paired-end and single-end realignment and improperly paired with high score reads.
    Check all pairs against original alignment; if read is unmapped or mapped with a mapping quality below
    a certain set threshold (10 by default) in realignment, replace realignment with original alignment
    (improperly paired with low alignment score).
    Output a dictionary of read pairs.
    """
    
    # create dictionary for pairs from realignment
    alignment_dict = {**get_pairs_from_pairend_alignment(directory), **get_pairs_from_singleend_alignment(directory)}
    print('Total: {} pairs of reads found from realignment'.format(len(alignment_dict)))
    
    # check new alignment dict against original alignments and 
    # replace unmapped/low score alignments with original alignment
    # check for alignment file and index file
    try:
        filename = glob.glob(os.path.join(directory, '*final.improper_lowq.bam'))[0]
    except:
        print('Alignment file not found!\nExiting from generate_alignment_table')
        exit()
    if not os.path.isfile(filename+'.bai'):
        pysam.samtools.index(filename)
    else:
        if os.path.getctime(filename+'.bai') < os.path.getctime(filename):
            pysam.samtools.index(filename)
    
    alignmentfile = pysam.AlignmentFile(filename, 'rb')
    count = 0
    for read in alignmentfile.fetch(refname):
        qname = read.query_name
        if qname in alignment_dict:
            if read.is_read1:
                if alignment_dict[qname][0] == 'unmapped' or alignment_dict[qname][2] < mapq_threshold:
                    alignment_dict[qname][:4] = [alignmentfile.get_reference_name(read.reference_id),
                                                 read.reference_start,
                                                 read.mapping_quality,
                                                 read.is_reverse]
                    alignment_dict[qname][-1] = False
            else:
                if alignment_dict[qname][4] == 'unmapped' or alignment_dict[qname][6] < mapq_threshold:
                    alignment_dict[qname][4:8] = [alignmentfile.get_reference_name(read.reference_id),
                                                  read.reference_start,
                                                  read.mapping_quality,
                                                  read.is_reverse]
                    alignment_dict[qname][-1] = False
        # sanity check: all reads should have been realigned and be found in alignment_dict
        else:
            count += 1
            print('Warning: read {} not found in realignment. Skip and continue'.format(qname))
        
    if count > 0:
        print('Warning: {} reads not found in realignment!')
            
    return alignment_dict


