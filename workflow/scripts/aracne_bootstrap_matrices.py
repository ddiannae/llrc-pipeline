
import sys
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")

matrix_cancer = snakemake.input[0]
matrix_normal = snakemake.input[1]
number_of_bootstraps = int(snakemake.params[0])
outdir = snakemake.params[1] + "/correlation/bootstrap_samples"
cond = snakemake.params[2]


print(f"Reading matrices")

if cond == "cancer":
    matrix_normal = pd.read_csv(matrix_normal, sep="\t", index_col=0, nrows=1)
    num_columns = matrix_normal.shape[1]
    for_bootstrap = pd.read_csv(matrix_cancer, sep="\t", index_col=0)
else:
    matrix_cancer = pd.read_csv(matrix_cancer, sep="\t", index_col=0, nrows=1)
    num_columns = matrix_cancer.shape[1]
    for_bootstrap = pd.read_csv(matrix_normal, sep="\t", index_col=0)
print(f"Number of columns: {num_columns}")

print("Generating bootstrap samples")

for i in range(number_of_bootstraps):
    sample_indices = pd.Series(range(num_columns)).sample(n=num_columns, replace=True).values
    bootstrap_matrix = for_bootstrap.iloc[:, sample_indices]
    out_matrix = f"{outdir}/matrix_{cond}_bootstrap_{i+1}.tsv"
    bootstrap_matrix.to_csv(out_matrix, sep="\t")
    print(f"Saved bootstrap sample {i+1} to {out_matrix}")

