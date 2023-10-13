rule filter_low_expression:
    input:
        config["datadir"]+"/{tissue}/rdata/raw_full_unfiltered.RData",
    output:
        config["datadir"]+"/{tissue}/rdata/raw_full.RData"
    log:
        config["datadir"]+"/{tissue}/log/filter_low.log"
    script:
        "../scripts/filterLowExpression.R"

rule qc:
    input:
        config["datadir"]+"/{tissue}/rdata/{plots_type}_full.RData"
    output:
        config["datadir"]+"/{tissue}/plots/{plots_type}/rna_composition.png"
    params:
        plots_dir=get_plots_dir
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_qc.log"
    script:
        "../scripts/NOISeqPlots.R"


rule density:
    input:
        config["datadir"]+"/{tissue}/rdata/{plots_type}_full.RData",
        config["datadir"]+"/{tissue}/plots/{plots_type}/rna_composition.png"
    output:
        outliers=config["datadir"]+"/{tissue}/results/{plots_type}_outliers.tsv",
        density=config["datadir"]+"/{tissue}/plots/{plots_type}/density.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_density.log"
    script:
        "../scripts/densityPlot.R"

rule filter_outliers:
    input:
        config["datadir"]+"/{tissue}/rdata/{plots_type}_full.RData",
        outliers=config["datadir"]+"/{tissue}/results/{plots_type}_outliers.tsv"
    output:
        config["datadir"]+"/{tissue}/rdata/{plots_type}.RData"
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_remove_outliers.log"
    script:
        "../scripts/filterOutliers.R"

rule pca:
    input:
        config["datadir"]+"/{tissue}/rdata/{plots_type}.RData"
    output:
        score=config["datadir"]+"/{tissue}/plots/{plots_type}/pca_score.png",
        loading=config["datadir"]+"/{tissue}/plots/{plots_type}/pca_loading.png",
        variance=config["datadir"]+"/{tissue}/plots/{plots_type}/pca_variance.png"
    params:
        plots_dir=get_plots_dir
    log:
        config["datadir"]+"/{tissue}/log/{plots_type}_pca.log"
    script:
        "../scripts/PCA.R"
