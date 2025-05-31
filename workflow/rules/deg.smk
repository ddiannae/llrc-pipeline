### Snakemake rule for differential gene expression analysis
### with deseq2
rule get_deg_deseq:
    input:
        raw=config["datadir"]+"/{tissue}/rdata/raw.RData",
        deseq2=config["datadir"]+"/{tissue}/rdata/deseq2.RData"
    output:
        deg_results=config["datadir"]+"/{tissue}/deg/deg_results.tsv"
    threads: 4
    log:
        config['datadir']+"/{tissue}/log/deg.log"
    script:
        "../scripts/deseq2DEG.R"

