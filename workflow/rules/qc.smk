## Snakemake rules for qc steps

## This rule removes genes with mean expression below 10
rule filter_low_expression:
    input:
        config["datadir"]+"/{tissue}/rdata/raw_full_unfiltered.RData",
    output:
        config["datadir"]+"/{tissue}/rdata/raw_full.RData"
    log:
        config["datadir"]+"/{tissue}/log/filter_low.log"
    script:
        "../scripts/filterLowExpression.R"

## This rule creates QC plots with NOISeq
rule qc:
    input:
        config["datadir"]+"/{tissue}/rdata/{data_format}_full.RData"
    output:
        config["datadir"]+"/{tissue}/plots/{data_format}/rna_composition.png"
    params:
        plots_dir=get_plots_dir
    log:
        config["datadir"]+"/{tissue}/log/{data_format}_qc.log"
    script:
        "../scripts/NOISeqPlots.R"

## This rule creates a density plot of the data and identifies outliers
rule density:
    input:
        config["datadir"]+"/{tissue}/rdata/{data_format}_full.RData",
        config["datadir"]+"/{tissue}/plots/{data_format}/rna_composition.png"
    output:
        outliers=config["datadir"]+"/{tissue}/results/{data_format}_outliers.tsv",
        density=config["datadir"]+"/{tissue}/plots/{data_format}/density.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/log/{data_format}_density.log"
    script:
        "../scripts/densityPlot.R"

## This rule filters outliers based on the density plot
rule filter_outliers:
    input:
        config["datadir"]+"/{tissue}/rdata/{data_format}_full.RData",
        outliers=config["datadir"]+"/{tissue}/results/{data_format}_outliers.tsv"
    output:
        config["datadir"]+"/{tissue}/rdata/{data_format}.RData"
    log:
        config["datadir"]+"/{tissue}/log/{data_format}_remove_outliers.log"
    script:
        "../scripts/filterOutliers.R"

## This rule performs PCA on the filtered data
rule pca:
    input:
        config["datadir"]+"/{tissue}/rdata/{arsyn}{data_format}.RData"
    output:
        score=config["datadir"]+"/{tissue}/plots/{data_format}/{arsyn}pca_score.png",
        loading=config["datadir"]+"/{tissue}/plots/{data_format}/{arsyn}pca_loading.png",
        variance=config["datadir"]+"/{tissue}/plots/{data_format}/{arsyn}pca_variance.png"
    params:
        plots_dir=get_plots_dir
    log:
        config["datadir"]+"/{tissue}/log/{arsyn}{data_format}_pca.log"
    script:
        "../scripts/PCA.R"
