## Snakemake rules for downloading GDC data
rule get_gdc_files:
    input:
        config["datadir"]+"/{tissue}/manifests/{type}-rna_counts.txt"
    output:
        config["datadir"]+"/{tissue}/raw/{type}/downloads_done.txt"
    params:
        rawdir=get_raw_dir
    shell:
        """
        mkdir -p {params.rawdir};
        ./bin/gdc-client download -d {params.rawdir} -m {input} --retry-amount 3;
        touch {output}
        """
    
rule get_manifest:
    input:
        config["datadir"]+"/{tissue}/log/done.txt"
    output:
        config["datadir"]+"/{tissue}/manifests/{type}-files.tsv",
        config["datadir"]+"/{tissue}/manifests/{type}-rna_counts.txt",
    params:
        tissue="{tissue}",
        type="{type}"
    log:
        config['datadir']+"/{tissue}/log/query_gdc_{type}.log"
    script:
        "../scripts/queryGDC.py"

rule setup_log:
    output:
        config["datadir"]+"/{tissue}/log/done.txt"
    params:
        config["datadir"]+"/{tissue}/log"
    shell:
        """
        mkdir -p {params}
        touch {output}
        """
