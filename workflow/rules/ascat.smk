## Snakefile for ASCAT2 files from GDC
rule get_cnv_genes:
  input:
    ascat_matrix=config["datadir"]+"/{tissue}/results/{type}-ascat-matrix.tsv"
  output:
    cnv_genes=config["datadir"]+"/{tissue}/results/{type}-ascat-cnv-genes.tsv"
  log: 
    config["datadir"]+"/{tissue}/log/{tissue}-{type}-ascat-get-cnv-genes.log"
  script:
    "../scripts/getCNVGenes.R"

rule get_heatmap:
  input: 
    config["datadir"]+"/{tissue}/results/{type}-ascat-matrix_all.tsv"
  output:
    heatmap=config["datadir"]+"/{tissue}/plots/{type}-ascat-heatmap.png",
    pcmatrix=config["datadir"]+"/{tissue}/results/{type}-ascat-matrix.tsv"
  params:
    biomart=config["biomart"]
  threads:32
  log: 
    config["datadir"]+"/{tissue}/log/{tissue}-{type}-ascat-heatmap.log"
  script:
    "../scripts/getAscatHeatmap.R" 

## We need to run these two together because the output of the download_files
## tasks depends on the manifest and there is no easy way to specify this on 
## snakemake
rule get_ascat_matrix:
  input:
    ## Manifest file
    config["datadir"]+"/{tissue}/raw_ascat/{type}/downloads_done.txt"
  output: 
    config["datadir"]+"/{tissue}/results/{type}-ascat-matrix_all.tsv",
    config["datadir"]+"/{tissue}/results/{type}-ascat-files.tsv"
  params: 
    rawdir=get_raw_ascat_dir,
    type="{type}",
    tissue="{tissue}"
  threads: 32
  log: 
    config["datadir"]+"/{tissue}/log/{tissue}-{type}-ascat-matrix.log"
  script:
    "../scripts/getAscatMatrix.R"

rule get_ascat_files:
  input:
    config["datadir"]+"/{tissue}/manifests/{type}-ascat-manifest.txt"
  output: 
    config["datadir"]+"/{tissue}/raw_ascat/{type}/downloads_done.txt"
  params:
    rawdir=get_raw_ascat_dir
  shell:
    """
    mkdir -p {params.rawdir} 
    ./bin/gdc-client download -d {params.rawdir} -m {input} --retry-amount 3
    touch {output}
    """

rule get_ascat_manifest:
  input:
    config["datadir"]+"/{tissue}/cnvs/done.txt"
  output:
    config["datadir"]+"/{tissue}/manifests/{type}-ascat-files.tsv",
    config["datadir"]+"/{tissue}/manifests/{type}-ascat-manifest.txt",
    ## Example: data/breast/manifests/breast-tumor-ascat.txt"
  params:
    tissue="{tissue}",
    type="{type}",
  log:
    config["datadir"]+"/{tissue}/log/query_gdc_ascat_{type}.log"
  script:
    "../scripts/queryGDCascat.py" 

rule setup_dirs:
  output:
    config["datadir"]+"/{tissue}/cnvs/done.txt"
  params:
    config["datadir"]+"/{tissue}/cnvs"
  shell:
    """
    mkdir -p {params} 
    touch {output}
    """ 
