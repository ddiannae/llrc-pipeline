from MultiAracne import Aracne
import sys
import pathlib

sys.stdout = open(snakemake.log[0], "w")

exp_matrix = snakemake.input[0]
outmatrix = snakemake.output[0]
outdir = snakemake.params[0]
procs = int(snakemake.threads)

print(f"Saving files in {outdir}")
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    
ma = Aracne(exp_matrix)
ma.run(processes=procs, outdir=outdir, pval=1)
ma.build_triu_missing_genes(outdir, outmatrix)
