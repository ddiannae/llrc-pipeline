from MultiAracne import Aracne
import sys
import pathlib
import shutil
import os

sys.stdout = open(snakemake.log[0], "w")

number_of_bootstraps = int(snakemake.params[0])
samples_dir = snakemake.params[1] + "/correlation/bootstrap_samples"
cond = snakemake.params[2]

for i in range(number_of_bootstraps):
    print(f"Processing bootstrap sample {i+1}")
    exp_matrix = f"{samples_dir}/matrix_{cond}_bootstrap_{i+1}.tsv"
    outmatrix = f"{samples_dir}/aracne_{cond}_bootstrap_{i+1}.adj"
    outdir = f"{samples_dir}/output_bootstrap"
    procs = int(snakemake.threads)
    
    if os.path.exists(outmatrix):
        print(f"{outmatrix} already exists, skipping.")
        continue

    print(f"Saving files in {outdir}")
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
        
    ma = Aracne(exp_matrix)
    ma.run(processes=procs, outdir=outdir, pval=1)
    ma.build_triu_missing_genes(outdir, outmatrix)

    shutil.rmtree(outdir)
    print(f"Completed processing bootstrap sample {i+1}")
