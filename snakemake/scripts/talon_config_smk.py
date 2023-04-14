import pandas as pd
import glob
from snakemake.io import expand

file = glob.glob('*_samples.csv')
file = ' '.join(map(str,file))
metadata = pd.read_csv(file, sep=',')
metadata_df = pd.DataFrame(data = metadata)
sn_list = metadata_df['sample_file_name'].tolist()
exp_list = list()
for line in sn_list:
	exp_list.append(line.split('_')[0])

metadata_df1 = metadata_df[['sample_name', 'description', 'platform']]
metadata_df1 = metadata_df1.assign(path_file=expand(snakemake.config['sample']['labeled'],zip, exp=exp_list, sample=sn_list))

metadata_df1.to_csv(snakemake.output['talon_config'], sep=',', header=False, index=False)