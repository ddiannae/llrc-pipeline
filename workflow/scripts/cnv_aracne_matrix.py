from multiprocessing import Pool
from functools import partial
from MultiAracne import Aracne
import sys
import pathlib
import pandas as pd
import numpy as np
import os
import sys

sys.stdout = open(snakemake.log[0], "w")

exp_matrix = snakemake.input[0]
cnv_matrix = snakemake.input[1]
# exp_genes = snakemake.input[2]
outmatrix = snakemake.output[0]
outdir = snakemake.params[0] + "/correlation/output_" + snakemake.params[1] + "_" + snakemake.params[2]+"_ascat"
procs = int(snakemake.threads)

print(f"Saving files in {outdir}")
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

try:
    df_exp = pd.read_csv(exp_matrix, sep="\t")
    df_cnv = pd.read_csv(cnv_matrix, sep="\t")
    df_cnv["gene_id"] = df_cnv["gene_id"]. apply(lambda x: x + "_cnv")
except:
    with open(outmatrix, 'a'):
        os.utime(outmatrix, None)
    sys.exit(f"Not enough columns")

if(df_cnv.shape[1] < 15 or df_cnv.shape[1] < 15):
    with open(outmatrix, 'a'):
        os.utime(outmatrix, None)
    sys.exit(f"Not enough columns")

n_cols_exp = df_exp.shape[1]
n_cols_cnv = df_cnv.shape[1]
min_cols = min(n_cols_exp, n_cols_cnv)
print(f"Selecting {min_cols} columns")
random_cols_exp = np.random.choice(range(1, n_cols_exp), min_cols-1, replace=False)
random_cols_cnv = np.random.choice(range(1, n_cols_cnv), min_cols-1, replace=False)
colnames = ["gene_id"] + ["sample_" + str(x) for x in range(1, min_cols)]

df_cnv = df_cnv.iloc[:, np.r_[0, random_cols_cnv]]
df_exp = df_exp.iloc[:, np.r_[0, random_cols_exp]]

df_cnv.columns = colnames
df_exp.columns = colnames

def run_aracne_for_gene(gene, expr_df, cnv_df, outdir):
    print(f"\nRunning for gene: {gene}")
    file_name = outdir + "/" + gene + ".tsv"
    print(f"\nSaving {file_name} matrix") 
    df_merged = pd.concat([expr_df[expr_df["gene_id"] == gene], cnv_df])
    df_merged.to_csv(file_name, sep="\t", index=False)
    print(f"\nRunning aracne")
    ma = Aracne(file_name)
    ma.run_gene(gene = gene, outdir=outdir, pval=1)
    print(f"\nDone for gene")

file_names = [outdir+"/"+fn for fn in os.listdir(outdir) if fn.endswith(".adj")]
existing_genes = [x.split("/")[-1].replace(".adj", "") for x in file_names]
remaining_genes = list(set(df_exp["gene_id"]) -set(existing_genes)) 

print(f"\nRunning exp+ascat for {len(remaining_genes)} genes")
p = Pool(procs)
fparam = partial(run_aracne_for_gene, expr_df=df_exp, cnv_df=df_cnv, outdir=outdir)
result = p.map_async(fparam, remaining_genes)
p.close()
p.join()
print(f"\nAracne computations done")

ma = Aracne(snakemake.input[0])
ma.build_nm_matrix(outdir, outmatrix, df_cnv["gene_id"])
