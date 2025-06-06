## Snakemake rules to get a raw expression matrix

## This rule gets gene annotations in the expression matrix for a given tissue.
rule join_and_annotate:
    input: unpack(get_annot_input)
    output:
        annot_tsv = config["datadir"]+"/{tissue}/results/{tissue}-annotation.tsv",
        annot_rdata = config["datadir"]+"/{tissue}/rdata/annot.RData",
        raw_rdata=config["datadir"]+"/{tissue}/rdata/raw_full_unfiltered.RData"
    log:
        config['datadir']+"/{tissue}/log/join_annotate.log"
    params:
        tissue_dir=get_tissue_dir,
        is_xena=is_xena_tissue,
        new_annot=config["new_annot"]
    threads: 3
    script:
        "../scripts/addAnnotations.R"

## This rule gets the raw expression matrix and a samples dataframe for a given tissue.
## It concatenates all the downloaded files from GDC or Xena
rule get_raw_matrix:
    input: get_raw_matrix_input
    output: 
        matrix = config["datadir"]+"/{tissue}/results/{type}-matrix.tsv",
        samples = config["datadir"]+"/{tissue}/results/{type}-samples.tsv",
    params:
        raw_dir=get_raw_dir,
        tissue=get_tissue_name,
        type="{type}",
        is_xena=is_xena_tissue, 
        primary=get_xena_primary,
    threads: 3
    log:
        config['datadir']+"/{tissue}/log/raw_matrix_{type}.log"
    script:
            "../scripts/getRawMatrix.R"
   
