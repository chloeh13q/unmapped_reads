import sys
from process_alignment_table import *

# open input file or stdin
if len(sys.argv) == 2:
    directory = sys.argv[1]
else:
    print('Insufficient number of arguments supplied to gen_batch_table [directory]')
    sys.exit()

batchid = os.path.basename(directory)
print(f'Generating tables for batch ID {batchid}')

if not os.path.isdir(f'/scratch/groups/dpwall/personal/chloehe/unmapped_reads/batch_tables/{batchid}'):
    os.mkdir(f'/scratch/groups/dpwall/personal/chloehe/unmapped_reads/batch_tables/{batchid}')
virus_df, bacteria_df = make_combined_table(directory)
virus_df.to_csv(f'/scratch/groups/dpwall/personal/chloehe/unmapped_reads/batch_tables/{batchid}/{batchid}.viral.csv')
bacteria_df.to_csv(f'/scratch/groups/dpwall/personal/chloehe/unmapped_reads/batch_tables/{batchid}/{batchid}.bacterial.csv')
# viral_hg38_df.to_csv(f'/scratch/groups/dpwall/personal/chloehe/unmapped_reads/batch_tables/{batchid}/{batchid}.viral_hg38.csv')
