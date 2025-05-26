import glob

def get_output_files(wildcards):
    files = []
    if config["end"] == "raw":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/rdata/raw_full_unfiltered.RData")
    if config["end"] == "qc":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/plots/raw/pca_score.png")
    if config["end"] == "norm":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/results/tpm_name_cancer.tsv")
            files.append(config["datadir"]+"/"+t["name"]+"/results/deseq2_name_cancer.tsv")
            files.append(config["datadir"]+"/"+t["name"]+"/results/tpm_ensembl_cancer.tsv")
            files.append(config["datadir"]+"/"+t["name"]+"/results/deseq2_ensembl_cancer.tsv")
    elif config["end"] == "correlation":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/tpm_ensembl_cancer.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/tpm_ensembl_normal.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/deseq2_ensembl_cancer.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/deseq2_ensembl_normal.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/arsyn_tpm_ensembl_normal.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/arsyn_tpm_ensembl_cancer.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/arsyn_deseq2_ensembl_cancer.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/arsyn_deseq2_ensembl_normal.adj")
    elif config["end"] == "deg":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/deg/deg_results.tsv")
    elif config["end"] == "ascat":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/results/cancer-ascat-matrix.tsv")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/arsyn_tpm_ensembl_normal_ascat.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/correlation/arsyn_tpm_ensembl_cancer_ascat.adj")
            files.append(config["datadir"]+"/"+t["name"]+"/results/cancer-ascat-cnv-genes.tsv")
            files.append(config["datadir"]+"/"+t["name"]+"/results/normal-ascat-cnv-genes.tsv")
    elif config["end"] == "immune":
        for t in config["tissues"]:
            files.append(config["datadir"]+"/"+t["name"]+"/immune/cibersort_tr4_heatmap.png")
            for j in ["CD31", "CD10", "CD45", "EPCAM"]:
                files.append(config["datadir"]+"/"+t["name"]+"/cancer_cibersort_tr4/output/cibersort_"+j+".adj")
                files.append(config["datadir"]+"/"+t["name"]+"/normal_cibersort_tr4/output/cibersort_"+j+".adj")
    return files

def get_cibersort_matrix(wildcards):
    if wildcards.tissue == "bladder" and wildcards.type == "normal":
        return config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/CIBERSORTxHiRes_NA_{cell}_Window10.txt"
    else:
        return config["datadir"]+"/{tissue}/{type}_cibersort_{matrix}/output/CIBERSORTxHiRes_NA_{cell}_Window16.txt"

def get_signature_matrix(wildcards):
    if wildcards.matrix == "lm22":
        return config["lm22"]
    else:
        return config["tr4"]

## raw.smk
def get_annot_input(wildcards):
    preffix = f'{config["datadir"]}/{wildcards.tissue}/results/'
    input = {
    "normal_targets": f'{preffix}normal-samples.tsv',
    "cancer_targets": f'{preffix}cancer-samples.tsv',
    "normal_matrix":  f'{preffix}normal-matrix.tsv',
    "cancer_matrix":  f'{preffix}cancer-matrix.tsv'
    }
    
    if is_xena_tissue(wildcards):
        xena_dir = get_xena_dir(wildcards)
        input["xena_annot"] =  f'{xena_dir}/annot.tsv'
        return input
    else:
        input["gdc_annot"]=config["gdc_annot"]
        return input

def get_raw_matrix_input(wildcards):
    if is_xena_tissue(wildcards):
        xena_dir = get_xena_dir(wildcards)
        tissue_dir = get_tissue_dir(wildcards)
        return [xena_dir +"/counts.gz", xena_dir +"/samples.txt.gz", tissue_dir+"/log/done.txt"]
    else:
        return [f'{config["datadir"]}/{wildcards.tissue}/raw/{wildcards.type}/downloads_done.txt']

def get_raw_ascat_dir(wildcards):
    return f'{get_tissue_dir(wildcards)}/raw-ascat/{wildcards.type}'

def get_raw_dir(wildcards):
    return f'{get_tissue_dir(wildcards)}/raw/{wildcards.type}'

def get_plots_dir(wildcards):
    return f'{get_tissue_dir(wildcards)}/plots/{wildcards.data_format}'

def get_xena_dir(wildcards):
    return f'{config["datadir"]}/{config["xenadir"]}'

def is_xena_tissue(wildcards):
    return wildcards.tissue in config["xena"]
     
def get_xena_primary(wildcards):
    if is_xena_tissue(wildcards):
        if wildcards.type == "normal":
            return get_normal_tissue(wildcards)
        else:
            return get_cancer_tissue(wildcards)
    else:
        return ""

def get_normal_tissue(wildcards):
    return [x["normal"] for x in config["tissues"] if x["name"] == wildcards.tissue][0]

def get_cancer_tissue(wildcards):
    return [x["cancer"] for x in config["tissues"] if x["name"] == wildcards.tissue][0]

def get_xena_extended_type(wildcards):
    return [x[f'sample_type_{wildcards.type}'] if f'sample_type_{wildcards.type}' in x else None for x in config["tissues"] if x["name"] == wildcards.tissue][0]

def get_tissue_name(wildcards):
    return [x[f'tissue_name_{wildcards.type}'] if f'tissue_name_{wildcards.type}' in x else x["name"] for x in config["tissues"] if x["name"] == wildcards.tissue][0]

## gdc.smk
def get_tissue_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}'

