source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

python -m snakemake -s ./snakemake/lr_blk.smk --configfile ./snakemake/lr_blk_config.yaml -j 9 \
--cluster 'sbatch -A ACCOUNT --mem {resources.mem_gb}G --cpus-per-task {threads} --output=$PWD/slurm_logs/slurm-%j.out' \
--shadow-prefix smk_tmp --use-conda  \
--latency-wait 60 \