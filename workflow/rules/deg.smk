rule get_deg:
    input:
        full=config["datadir"]+"/{tissue}/rdata/{step1}_{step2}_{step3}_{arsyn}_full.RData",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    output:
        deg_results=config["datadir"]+"/{tissue}/deg/{step1}_{step2}_{step3}_{arsyn}_deg_results.tsv"
    params:
        tissue="{tissue}",
        deg_dir=config["datadir"]+"/{tissue}/deg"
    log:
        config['datadir']+"/{tissue}/log/{step1}_{step2}_{step3}_{arsyn}_deg.log"
    script:
        "../scripts/deg.R"
     
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

