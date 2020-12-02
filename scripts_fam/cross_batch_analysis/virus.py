import pandas as pd
import os


def find_viral_integration_events(file, count_threshold=20, percent_threshold=0.2):
    
    df = pd.read_csv(file, index_col=0)
    df_copy = df.groupby('sampleID').sum()
    df_copy['sampleID'] = df_copy.index
    df_melt = df_copy.melt(id_vars='sampleID')
    df_stats = df_melt.groupby('variable').agg(['mean', 'max'])
    baseline_dict = df_stats[('value',  'mean')].to_dict()
    
    potential = []
    for virus in list(baseline_dict.keys()):
        potential.append(df_melt[(df_melt['variable']==virus) & (df_melt['value'] > max(count_threshold, baseline_dict[virus] * 1.5))])
    potential = pd.concat(potential)
    
    aligned_pairs = []
    # iterate through every sample + virus pair
    for row in potential.itertuples():
        # retreive alignments for mates
        counts = df.loc[(df['sampleID']==row.sampleID) & (df[row.variable] > 0), [row.variable, 'sampleID']]
        max_percent = counts[row.variable].max() / counts[row.variable].sum()
        max_id = counts[row.variable].idxmax()
        if max_percent > percent_threshold:
            new_row = counts.loc[max_id, :]
            aligned_pairs.append({'sampleID': row.sampleID, 
                                  'batch': os.path.basename(file).replace('.viral_hg38.csv', ''),
                                  'virus': row.variable, 
                                  'mate': new_row.name, 
                                  'count': new_row.loc[row.variable]})
    aligned_pairs = pd.DataFrame(aligned_pairs)
    
    return aligned_pairs