import sys
from get_pairs import *

# open input file or stdin
if len(sys.argv) == 2:
    directory = sys.argv[1]
else:
    print('Insufficient number of arguments supplied to gen_alignment_table [directory]')
    sys.exit()

sampleid = os.path.basename(directory)
print('Generating alignment table for sample ID {}'.format(sampleid))
alignment_dict = generate_alignment_dict(directory)
alignment_df = pd.DataFrame.from_dict(alignment_dict, orient='index')
alignment_df.columns = ['R1_ref', 'R1_start', 'R1_MAPQ', 'R1_is_reverse', 
                        'R2_ref', 'R2_start', 'R2_MAPQ', 'R2_is_reverse',
                        'is_proper_pair']

alignment_df.to_csv(os.path.join(directory, sampleid+'.final_alignment_table.csv'))
