
import sys
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")

matrix_cancer = snakemake.input[0]
matrix_normal = snakemake.input[1]
number_of_bootstraps = int(snakemake.params[0])
outdir = snakemake.params[1] + "/correlation/bootstrap_samples"


print(f"Reading matrices")
matrix_cancer = pd.read_csv(matrix_cancer, sep="\t", index_col=0)
matrix_normal = pd.read_csv(matrix_normal, sep="\t", index_col=0, nrows=1)

num_columns = matrix_normal.shape[1]
print(f"Number of columns: {num_columns}")

print("Generating bootstrap samples")

for i in range(number_of_bootstraps):
    sample_indices = pd.Series(range(num_columns)).sample(n=num_columns, replace=True).values
    bootstrap_cancer = matrix_cancer.iloc[:, sample_indices]
    out_cancer = f"{outdir}/matrix_cancer_bootstrap_{i+1}.tsv"
    bootstrap_cancer.to_csv(out_cancer, sep="\t")
    print(f"Saved bootstrap sample {i+1} to {out_cancer}")

