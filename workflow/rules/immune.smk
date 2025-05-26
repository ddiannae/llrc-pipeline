rule cibersort:
    input:
        config["datadir"]+"/{tissue}/results/{type}_tpm_genesymbols.tsv" 
    output:
        config["datadir"]+"/{tissue}/immune/{type}_{matrix}_cibersort.tsv",
    params:
        signature_matrix=get_signature_matrix,
        cibersort=config["cibersort"],
        outdir=config["datadir"]+"/{tissue}/immune"
    log:
        config['datadir']+"/{tissue}/log/{type}_xcell_{matrix}_cibersort.log"
    threads: 18
    script:
        "../scripts/cibersort.R"
     
rule cibersort_xcell_all:
    input:
        cancer=config["datadir"]+"/{tissue}/results/cancer_tpm_genesymbols.tsv" ,
        normal=config["datadir"]+"/{tissue}/results/normal_tpm_genesymbols.tsv" 
    output:
        xcell=config["datadir"]+"/{tissue}/immune/all_{matrix}_xcell.tsv",
        cibersort=config["datadir"]+"/{tissue}/immune/all_{matrix}_cibersort.tsv",
    params:
        signature_matrix=get_signature_matrix,
        cibersort=config["cibersort"],
        outdir=config["datadir"]+"/{tissue}/immune"
    log:
        config['datadir']+"/{tissue}/log/xcell_{matrix}_cibersort.log"
    threads: 18
    script:
        "../scripts/cibersort_xcell_all.R"

rule cibersort_xcell_heatmaps:
    input:
        xcell=config["datadir"]+"/{tissue}/immune/all_{matrix}_xcell.tsv",
        cibersort=config["datadir"]+"/{tissue}/immune/all_{matrix}_cibersort.tsv"
    output:
        xcell=config["datadir"]+"/{tissue}/immune/xcell_{matrix}_heatmap.png",
        cibersort=config["datadir"]+"/{tissue}/immune/cibersort_{matrix}_heatmap.png"
    params:
        tissue="{tissue}"
    log:
        config['datadir']+"/{tissue}/log/xcell_{matrix}_cibersort_heatmap.log"
    threads: 7
    script:
        "../scripts/cibersort_xcell_heatmap.R"
     
rule set_cibersort_files:
    input:
        config["datadir"]+"/{tissue}/results/{type}_tpm_genesymbols.tsv" 
    output:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/input/tpm_genesymbols.tsv",
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/input/signature_matrix.txt"
    params:
        input_dir=config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/input",
        output_dir=config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output",
        signature_matrix=get_signature_matrix
    log:
        config['datadir']+"/{tissue}/log/{type}_{matrix}_set_cibersort_files.log" 
    shell:
        """
        mkdir -p {params.input_dir} {params[1]}
        cp {input[0]} {output[0]}
        cp {params.signature_matrix} {output[1]}
        """

rule run_cibersort_hires:
    input:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/input/tpm_genesymbols.tsv",
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/input/signature_matrix.txt"
    output:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/done.txt"
    params:
        input_dir=config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/input",
        output_dir=config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output",
        username=config["cibersort_username"],
        token=config["cibersort_token"]
    log:
        config['datadir']+"/{tissue}/log/{type}_{matrix}_run_cibersort.log" 
    threads: 30
    shell:
        """
        docker run -v {params.input_dir}:/src/data -v {params.output_dir}:/src/outdir cibersortx/hires --username {params.username} --token {params.token} --mixture tpm_genesymbols.tsv --sigmatrix signature_matrix.txt --threads {threads}
        touch {params.output_dir}/done.txt
        """

rule filter_cibersort_matrix:
    input:
        get_cibersort_matrix,
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/done.txt"
    output:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/cibersort_{cell}.tsv"
    log:
        config['datadir']+"/{tissue}/log/{type}_{matrix}_{cell}_filter.log"
    script:
        "../scripts/filter_cibersort.R"
        

rule get_cibersort_aracne:
    input:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/cibersort_{cell}.tsv"
    output:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/cibersort_{cell}.adj"
    params:
        config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/aracne_{cell}_output"
    singularity:
        config["aracne_singularity"]
    log:
        config['datadir']+"/{tissue}/log/{type}_{matrix}_{cell}_aracne.log"
    threads: 39
    script:
        "../scripts/cibersort_aracne.py"


