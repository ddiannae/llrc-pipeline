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
