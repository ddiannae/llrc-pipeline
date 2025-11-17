## Snakemake rules for correlation analysis using ARACNe

## This rule generates a mutual information matrix from aracne for gene expression - gene copy number data.
rule aracne_ascat:
    input:
        config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_{type}.tsv",
        config["datadir"]+"/{tissue}/results/{type}-ascat-matrix.tsv",
        config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_genes.tsv",
    output:
        config["datadir"]+"/{tissue}/correlation/{arsyn}{data_format}_{gene_id}_{type}_ascat.adj",
    singularity:
        config["aracne_singularity"]
    params:
        get_tissue_dir,
        "{arsyn}{data_format}", 
        "{type}" 
    threads: 16
    log:
        config["datadir"]+"/{tissue}/log/{arsyn}{data_format}_{gene_id}_{type}_ascat_aracne.log"
    script:
        "../scripts/cnv_aracne_matrix.py"

## This rule generates a mutual information matrix from aracne for gene expression - gene expression data. 
rule aracne:
    input:
        config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_{type}.tsv",
    output:
        config["datadir"]+"/{tissue}/correlation/{arsyn}{data_format}_{gene_id}_{type}.adj"
    singularity:
        config["aracne_singularity"]
    params:
        get_tissue_dir,
        "{arsyn}{data_format}", 
        "{type}" 
    threads: 16
    log:
        config["datadir"]+"/{tissue}/log/{arsyn}{data_format}_{gene_id}_{type}_aracne.log"
    script:
        "../scripts/aracne_matrix.py"

rule aracne_bootsrap:
    input:
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples/matrix_cancer_bootstrap_"+str(config["bootstrap_samples"])+".tsv"
    output:
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples/aracne_cancer_bootstrap_"+str(config["bootstrap_samples"])+".adj"
    threads: 8
    singularity:
        config["aracne_singularity"]
    params:
        config["bootstrap_samples"],
        get_tissue_dir
    log:
        config["datadir"]+"/{tissue}/log/aracne_bootstrap.log"
    script:
        "../scripts/aracne_bootstrap.py"

rule generate_aracne_bootstrap:
    input:
        config["datadir"]+"/{tissue}/results/deseq2_ensembl_cancer.tsv",
        config["datadir"]+"/{tissue}/results/deseq2_ensembl_normal.tsv"
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples/done.txt"
    output:
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples/matrix_cancer_bootstrap_"+str(config["bootstrap_samples"])+".tsv"
    params:
        config["bootstrap_samples"],
        get_tissue_dir,
    threads: 4
    log:
        config["datadir"]+"/{tissue}/log/generate_aracne_bootstrap.log"
    script:
        "../scripts/aracne_bootstrap_matrices.py"

rule setup_bootstrap:
    output:
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples/done.txt"
    params:
        config["datadir"]+"/{tissue}/correlation/bootstrap_samples"
    shell:
        """
        mkdir -p {params}
        touch {output}
        """