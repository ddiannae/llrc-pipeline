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

