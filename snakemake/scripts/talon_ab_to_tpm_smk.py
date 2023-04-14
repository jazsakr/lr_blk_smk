import pandas as pd

df = pd.read_csv(snakemake.input['talon_af'], sep='\t', header=0)

def get_dataset_cols(df,
                    sample=None):
    """
    Get the names of the dataset columns from a TALON abundance file
    Parameters:
        df (pandas DataFrame): TALON ab dataframe
        sample (str): Tissue name of samples to pass
    Returns:
        dataset_cols (list of str): Names of dataset columns
    """

    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                           'annot_transcript_id', 'annot_gene_name',
                           'annot_transcript_name', 'n_exons', 'length',
                           'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_columns ]

    # subset by sample if necessary
    if sample != 'all':
        if sample == 'mouse_match':
            pass
        elif sample:
            dataset_cols = [x for x in dataset_cols if sample in x]
    return dataset_cols

datasets_cols = get_dataset_cols(df)

# calculate tpm
for d in dataset_cols:
    tpm_col = f'{d}_tpm'
    df['total'] = df[d].sum()
    df[tpm_col] =  (df[d]*(10**6))/df['total']
    df.drop('total', axis=1, inplace=True)

# replace your counts columns with the tpm columns
df.drop(dataset_cols, axis=1, inplace=True)

#print(df)

df.to_csv(snakemake.output['talon_tpm'], sep='\t')