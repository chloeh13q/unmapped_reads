#!/oak/stanford/groups/dpwall/computeEnvironments/miniconda/bin/python3

import sys
import urllib.request as urllib2
import re
import time

# input file should contain one column with a list of NCBI Accession #s

# open input file or stdin
if len(sys.argv) == 3:
    inpath = sys.argv[1]
    outpath = sys.argv[2]
    infile = open(inpath, "rb")
else:
    print('Insufficient number of arguments supplied to get_fasta_by_accn.py [inpath] [outfoler]')
    sys.exit()

# base url
base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

accns = []
# corresponding nuccore IDs, which are needed for query
nuccoreIDs = []

for line in infile:
    accession_number = line.strip().decode('ascii')
    accession_number += '[accn]'
    try:
        esearch_url = base + 'esearch.fcgi?db=nuccore&term=' + accession_number
        with urllib2.urlopen(esearch_url) as url:
            response = url.read()
            time.sleep(0.2)
        accns.append(re.sub('\\[accn\\]', '', accession_number))
        for field in response.split():
            if b'<Id>' in field:
                nuccoreIDs.append(re.sub('</*Id>', '', field.decode('ascii')))
    except urllib2.HTTPError:
        print('esearch error: Unable to retrieve genome for accession #' + accession_number)

infile.close()        
time.sleep(2)

for ii, ID in enumerate(nuccoreIDs):
    try:
        efetch_url = base + 'efetch.fcgi?db=nuccore&id=' + ID + '&rettype=fasta'
        with urllib2.urlopen(efetch_url) as url:
            fasta = url.read()
            time.sleep(0.2)
        with open(outpath + accns[ii] + ".fa", 'wb') as out:
            out.write(fasta)
    except urllib2.HTTPError:
        print('efetch error: Unable to retrieve genome for accession #' + accession_number)
        
