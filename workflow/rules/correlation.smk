rule aracne_ascat:
    input:
        config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_{type}.tsv",
        config["datadir"]+"/{tissue}/results/{type}-ascat-matrix.tsv",
        config["datadir"]+"/{tissue}/results/{arsyn}{data_format}_{gene_id}_genes.tsv",
    output:
        config["datadir"]+"/{tissue}/correlation/{arsyn}{data_format}_{gene_id}_{type}_ascat.adj",
        #config["datadir"]+"/{tissue}/correlation/{arsyn}{data_format}_{gene_id}_{type}_ascat_matrix.tsv"
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

