## Snakemake rules for normalization tasks

## This rule runs arsyn for batch correction.
## Maria j. Nueda, Alberto Ferrer, Ana Conesa, ARSyN: a method for the identification and removal of
## systematic noise in multifactorial time course microarray experiments,
## Biostatistics, Volume 13, Issue 3, July 2012, Pages 553–566, https://doi.org/10.1093/biostatistics/kxr042
rule arsyn:
    input:
        config["datadir"]+"/{tissue}/rdata/{data_format}.RData",
    output:
        config["datadir"]+"/{tissue}/rdata/arsyn_{data_format}.RData",
    threads: 8
    log:
        config['datadir']+"/{tissue}/log/arsyn_{data_format}.log"
    script:
        "../scripts/runArsyn.R"

## This rule performs DESeq2 normalization.
rule deseq_normalization:
    input:
        config["datadir"]+"/{tissue}/rdata/raw.RData"
    output:
        config["datadir"]+"/{tissue}/rdata/deseq2_full.RData",
        dds=config["datadir"]+"/{tissue}/rdata/dds.RDS"
    log:
        config["datadir"]+"/{tissue}/log/deseq2_normalization.log"
    threads: 8
    script:
        "../scripts/deseq2Normalization.R"

## This rule obtains the TPM values from the raw counts.
rule tpm_normalization:
    input:
        config["datadir"]+"/{tissue}/rdata/raw.RData"
    output:
        config["datadir"]+"/{tissue}/rdata/tpm_full.RData",
    log:
        config['datadir']+"/{tissue}/log/get_tpms.log"
    threads: 8
    script:
        "../scripts/tpmNormalization.R"

## This rule saves tsv files with the normalized matrices for cancer and normal samples.
rule get_norm_matrix:
    input:
        config["datadir"]+"/{tissue}/rdata/{arsyn}{data_format}.RData",
        config["datadir"]+"/{tissue}/plots/{data_format}/{arsyn}pca_score.png"
    output:
        cancer=config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_cancer.tsv",
        normal=config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_normal.tsv",
        genes=config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_genes.tsv",
    params:
        gene_id="{gene_id}"
    threads: 4
    log:
        config['datadir']+"/{tissue}/log/{arsyn}{data_format}_{gene_id}_matrices.log"
    script:
        "../scripts/getNormMatrix.R"

