# alternative haplotype: NCBI GenBank nucleotide database, accession #s MH533022-MH534863

import sys
import re

if len(sys.argv) == 3:
    begin_accn = sys.argv[1]
    end_accn = sys.argv[2]
else:
    print('Insufficient number of arguments supplied to list_accn.py [begin_accn] [end_accn]')
    sys.exit()

accn = []
prefix = re.sub('[0-9]+', '', begin_accn)
begin_accn = int(re.sub('[a-zA-Z]+', '', begin_accn))
end_accn = int(re.sub('[a-zA-Z]+', '', end_accn))

with open('accn.txt', 'wb') as out:
    for i in range(begin_accn, end_accn+1):
        out.write(str.encode(prefix + str(i) + '\n'))
        
