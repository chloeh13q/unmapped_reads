import pysam
import os
import re
from collections import Counter

def count_unmapped_reads(filepath, window=10000000, chromosomal=True):
    
    directory, filename = os.path.split(filepath)
    
    if 'bam' in filename:
        suffix = '.bai'
        readmode = 'rb'
    elif 'cram' in filename:
        suffix = '.crai'
        readmode = 'rc'
    else:
        print('File type not supported')
        exit()
        
    # samtools index if index file not found
    if not os.path.isfile(os.path.join(directory, filename+suffix)):
        pysam.samtools.index(filepath)
        
    samfile = pysam.AlignmentFile(filepath, readmode)
    
    ref = []
    all_unmapped_mate_counts, all_unmapped_self_counts = [], []
    
    # default: retrieve only chromosomal reads
    if chromosomal:
        num_ref = 24
    else:
        num_ref = len(samfile.references)
        
    # for each window size of base pairs in samfile references (@SQ)
    for i in range(num_ref):
        refname = samfile.references[i]
        seqlen = samfile.lengths[i]
        for j in range(1, seqlen, window):
            stop = j+window-1 if j+window-1 < samfile.lengths[i] else samfile.lengths[i]
            unmapped_mate = 0
            unmapped_self = 0
            
            if not re.match('chr', refname):
                ref.append('chr' + refname + ':' + str(j) + '-' + str(stop))
            else:
                ref.append(refname + ':' + str(j) + '-' + str(stop))
                
            # count unmapped reads
            for read in samfile.fetch(refname, j, stop):
                if not read.is_duplicate and not read.is_secondary and not read.is_supplementary:
                    if read.mate_is_unmapped and not read.is_unmapped:
                        unmapped_mate += 1
                    if read.is_unmapped and not read.mate_is_unmapped:
                        unmapped_self += 1
            all_unmapped_mate_counts.append(unmapped_mate)
            all_unmapped_self_counts.append(unmapped_self)

    samfile.close()
    print('Finished counting unmapped reads for ' + filepath)
    return ref, all_unmapped_mate_counts, all_unmapped_self_counts

def count_improper_pairs(filepath, window=100000000):
    
    directory, filename = os.path.split(filepath)
    
    if 'bam' in filename:
        suffix = '.bai'
        readmode = 'rb'
    elif 'cram' in filename:
        suffix = '.crai'
        readmode = 'rc'
    else:
        print('File type not supported')
        exit()
        
    # samtools index if index file not found
    if not os.path.isfile(os.path.join(directory, filename+suffix)):
        pysam.samtools.index(filepath)
    samfile = pysam.AlignmentFile(filepath, readmode)
    
    # build a list of chromosome and loci references by window size
    # must retrieve all references here to avoid index out of bounds when indexing into the reference dictionary
    ref = []
    for i in range(len(samfile.references)): 
        refname = samfile.references[i]
        seqlen = samfile.lengths[i]
        for j in range(1, seqlen, window):
            stop = j+window-1 if j+window-1 < samfile.lengths[i] else samfile.lengths[i]
            ref.append(refname + ',' + str(j) + ',' + str(stop))
    
    output = dict()
    
    for reference in ref:
        refname, start, stop = reference.split(',')
        improper_pairs_counter = Counter()
        counter_list = []
        for r in ref:
            improper_pairs_counter[r] = 0
        for read in samfile.fetch(refname, int(start), int(stop)):
            if not read.is_duplicate and not read.is_secondary and not read.is_supplementary:
                if not read.is_proper_pair:
                    index_start = int(read.next_reference_start) // window * window + 1
                    search_key = read.next_reference_name + ',' + str(index_start)
                    index = list(filter(lambda keys: search_key in keys, improper_pairs_counter.keys()))[0]
                    improper_pairs_counter[index] += 1
        for r in ref:
            counter_list.append(improper_pairs_counter[r])
        output[reference] = counter_list
    
    samfile.close()
    print('Finished counting improper pairs for ' + filepath)
    return output

def count_virus_alignments(filepath):
    
    directory, filename = os.path.split(filepath)
    
    if 'bam' in filename:
        suffix = '.bai'
        readmode = 'rb'
    elif 'cram' in filename:
        suffix = '.crai'
        readmode = 'rc'
    else:
        print('File type not supported')
        exit()

    # samtools index if index file not found
    if not os.path.isfile(os.path.join(directory, filename+suffix)):
        pysam.samtools.index(filepath)
    samfile = pysam.AlignmentFile(filepath, readmode)
    
    # build a list of chromosome and loci references by window size
    counter = Counter()
    for i in range(len(samfile.references)):
        refname = samfile.references[i]
        count = 0
        for read in samfile.fetch(refname):
            if not read.is_supplementary:
                count += 1
        if not count == 0:
            counter[refname] = count
        
    samfile.close()
            
    return counter
