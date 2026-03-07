#' Simple analysis of single-cell data
#'
#' Simple analysis of scRNA-seq or CITE-seq data, from quality control to clustering and marker gene detection. 
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object where each column represents a single cell.
#' Rows are usually expected to contain genes or antibody-derived tags, see \code{rna.experiment=} and \code{adt.experiment=} for details.
#' @param rna.experiment Identity of the experiment containing the RNA data, when \code{x} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' If \code{TRUE}, the main experiment is assumed to contain the RNA data (see \code{\link[SingleCellExperiment]{mainExpName}}).
#' If a string is supplied, it is treated as the name of the alteranative experiment (see \code{\link[SingleCellExperiment]{altExps}}) containing the RNA data.
#' If \code{FALSE} or \code{NULL}, it is assumed that no RNA data is available.
#' @param adt.experiment Identity of the experiment containing the ADT data, when \code{x} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' If \code{TRUE}, the main experiment is assumed to contain the ADT data (see \code{\link[SingleCellExperiment]{mainExpName}}).
#' If a string is supplied, it is treated as the name of the alteranative experiment (see \code{\link[SingleCellExperiment]{altExps}}) containing the ADT data.
#' If \code{FALSE} or \code{NULL}, it is assumed that no ADT data is available.
#' @param subset.factor String specifying the name of the \code{\link[SummarizedExperiment]{colData}(se)} column containing a factor with which to subset the cells.
#' @param subset.levels Vector containing the subset of levels to retain in the factor specified by \code{subset.factor}.
#' @param block.field String specifying the name of the \code{\link[SummarizedExperiment]{colData}(se)} column containing the block assignment for each cell.
#' @param qc.mito.regex String containing a regular expression to identify mitochondrial genes from the row names.
#' If \code{symbol.field=} is specified, this is applied to the gene symbols instead.
#' If \code{NULL}, the mitochondrial genes are identified from the sequence names of \code{\link[SummarizedExperiment]{rowRanges}} matching \code{qc.mito.seqnames}. 
#' Only used if \code{rna.experiment=} indicates that RNA data is available.
#' @param qc.mito.seqnames Character vector containing the sequence names of the mitochondrial chromosome.
#' Only used if \code{rna.experiment=} indicates that RNA data is available and \code{qc.mito.regex=} is not specified. 
#' @param qc.igg.regex String containing a regular expression to identify IgG controls from the row names.
#' If \code{symbol.field=} is specified, this is applied to the gene symbols instead.
#' Only used if \code{adt.experiment=} indicates that ADT data is available.
#' @param qc.num.mads Integer specifying the number of median absolute deviations (MADs) with which to define a quality control (QC) filtering threshold.
#' Smaller values increase the stringency of the filter.
#' @param qc.filter Boolean indicating whether putative low-quality cells should be removed.
#' If \code{FALSE}, low-quality cells are identified but not removed prior to further cells.
#' @param num.hvgs Integer specifying the number of highly variable genes (HVGs) to retain for downstream analyses.
#' More HVGs capture more biological signal at the cost of capturing more technical noise and increasing computational work.
#' Only used if \code{rna.experiment=} indicates that RNA data is available.
#' @param num.pcs Integer specifying the number of top principal components (PCs) to retain for downstream analyses.
#' More PCs capture more biological signal at the cost of capturing more technical noise and increasing computational work.
#' @param mnn.num.neighbors Integer specifying the number of neighbors for batch correction with mutual nearest neighbors.
#' Larger values improve stability but reduce resolution for rare subpopulations.
#' @param mnn.num.steps Integer specifying the number of steps for the center of mass calculation during batch correction with mutual nearest neighbors.
#' Larger values improve intermingling of cells from different batches but increase the risk of merging the wrong subpopulations.
#' @param cluster.method Character vector specifying the clustering methods to run on the top PCs,
#' namely graph-based clustering (\code{"graph"}) or k-means clustering (\code{"kmeans"}).
#' Multiple methods may be specified here to generate multiple clusterings, but only the clusters from the first method are used for marker detection.
#' @param cluster.kmeans.k Integer specifying the number of clusters to generate from k-means clustering.
#' Only used if \code{"kmeans"} is in \code{cluster.method}.
#' @param cluster.graph.method String naming the community detection method to use in graph-based clustering.
#' These are roughly equivalent to the functions of the same name from the \pkg{igraph} R package.
#' (Note that \code{"multilevel"} is a synonym for the Louvain method.)
#' Only used if \code{"graph"} is in \code{cluster.method}.
#' @param cluster.graph.num.neighbors Integer specifying the number of nearest neighbors to use during construction of the shared-nearest neighbor graph. 
#' Larger values increase graph connectivity and decrease cluster resolution. 
#' Only used if \code{"graph"} is in \code{cluster.method}.
#' @param cluster.graph.resolution Number specifying the resolution to use in the multi-level or Leiden community detection algorithms.
#' Larger values usually result in a greater number of smaller clusters.
#' Only used if \code{"graph"} is in \code{cluster.method}.
#' @param reduced.dimensions Character vector specifying the dimensionality reduction algorithms to use for visualization.
#' This can be zero, one or both of \code{"tsne"} (for t-SNEs) and \code{"umap"} (for UMAPs).
#' @param tsne.perplexity Number specifying the perplexity to use in t-SNE.
#' Larger values increase the size of each cell's neighborhood and focus more on global structure. 
#' Only used if \code{"tsne"} is in \code{reduced.dimensions}.
#' @param umap.num.neighbors Integer specifying the number of neighbors to use in the UMAP.
#' Larger values increase connectivity and focus more on global structure. 
#' Only used if \code{"umap"} is in \code{reduced.dimensions}.
#' @param umap.min.dist Number specifying the minimum distance between points in the UMAP.
#' Larger values favor a more even distribution of cells throughout the low-dimensional space.
#' Only used if \code{"umap"} is in \code{reduced.dimensions}.
#' @param marker.effect.size String naming the effect size to use when ranking marker genes.
#' This should be one of:
#' \itemize{
#' \item \code{"cohens.d"}: Cohen's d, i.e., the standardized difference in log-expression.
#' This is analogous to the t-statistic in a two-group t-test and accounts for the variance within each group.
#' \item \code{"auc"}: AUC, i.e., the area under the curve.
#' This is closely related to the U statistic in a Wilcoxon rank sum test, and is a more robust/less sensitive alternative Cohen's d.
#' \item \code{"delta.mean"}: difference in the mean log-expression.
#' This is an estimate of the log-fold change.
#' \item \code{"delta.detected"}: difference in the proportion of cells with detected expression.
#' Large differences correspond to genes that are silent in one group and activated in another.
#' }
#' This is combined with \code{marker.summary} to determine the actual statistic for ranking, see \code{?\link[scrapper]{scoreMarkers}}.
#' @param marker.summary String specifying the summary statistic to use when ranking marker genes.
#' This should be one of:
#' \itemize{
#' \item \code{"min.rank"}: the minimum rank of each gene across all pairwise comparisons for a particular cluster.
#' This creates a ranking where the top genes are guaranteed to separate the cluster of interest from every other cluster.
#' \item \code{"mean"}: the mean effect size of each gene across all pairwise comparisons for a particular cluster.
#' This creates a ranking where the top genes are upregulated in the cluster of interest against the average of all other clusters.
#' \item \code{"median"}: the median effect size of each gene across all pairwise comparisons for a particular cluster.
#' This creates a ranking where the top genes are upregulated in the cluster of interest against most other clusters.
#' \item \code{"min"}: the minimum effect size of each gene across all pairwise comparisons for a particular cluster.
#' This creates a ranking where the top genes must be upregulated in the cluster of interest against all other clusters.
#' }
#' This is combined with \code{marker.effect} to determine the actual statistic for ranking, see \code{?\link[scrapper]{scoreMarkers}}.
#' @param marker.lfc.threshold Non-negative number specifying the log-fold change threshold to test against.
#' Larger values focus on marker genes with larger log-fold changes at the expense of those with smaller variances.
#' @param assay String or integer specifying the assay in \code{x} containing the count matrix.
#' For \link[SingleCellExperiment]{SingleCellExperiment} objects, the same assay is used for the main and alternative experiments.
#' @param num.threads Integer specifying the number of threads to use in the various computations.
#' @param output.dir String containing the path to an output directory in which to write the Rmarkdown report and save results.
#' @param metadata Named list of additional fields to add to each result's metadata.
#' @param author Character vector of authors.
#' @param dry.run Boolean indicating whether to perform a dry run.
#' This will write the Rmarkdown report without evaluating it.
#' @param symbol.field String specifying the name of the \code{\link[SummarizedExperiment]{rowData}} column that contains gene symbols.
#' This is added to the various gene-based results, e.g., marker detection tables.
#' If \code{NULL}, no symbols are added.
#' @param save.results Boolean indicating whether the results should also be saved to file.
#'
#' @return 
#' A Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir} that contains the analysis commands.
#'
#' If \code{dry.run=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{sce}, a \link[SingleCellExperiment]{SingleCellExperiment} object with dimensionality reduction and clustering results.
#' This may have fewer cells than the input \code{x} if \code{qc.filter=TRUE}.
#' \item \code{qc.rna}, a \link[S4Vectors]{DataFrame} object with RNA-based QC metrics for each cell in the input \code{x}.
#' Only returned if \code{rna.experiment=} indicates that RNA data is available.
#' \item \code{qc.adt}, a \link[S4Vectors]{DataFrame} object with ADT-based QC metrics for each cell in the input \code{x}.
#' Only returned if \code{adt.experiment=} indicates that ADT data is available.
#' \item \code{markers.rna}, a named list of \link[S4Vectors]{DataFrame} results containing marker gene statistics for each cluster.
#' Only returned if \code{rna.experiment=} indicates that RNA data is available.
#' \item \code{markers.adt}, a named list of \link[S4Vectors]{DataFrame} results containing marker tag statistics for each cluster.
#' Only returned if \code{adt.experiment=} indicates that ADT data is available.
#' }
#'
#' If \code{save.results=TRUE}, the results are saved in a \code{results} directory inside \code{output}.
#' 
#' If \code{dry.run=TRUE}, \code{NULL} is returned.
#' Only the Rmarkdown report is saved to file.
#'
#' @author Aaron Lun
#' 
#' @examples
#' # Getting an example dataset from the scRNAseq package.
#' library(scRNAseq)
#' se <- ZeiselBrainData()
#'
#' # Running the pipeline; omitting the save for speed here.
#' tmp <- tempfile()
#' output <- runSolo(se, qc.mito.regex="^mt-", output.dir=tmp, save.results=FALSE)
#'
#' # Checking the outputs.
#' list.files(tmp, recursive=TRUE)
#' output$sce
#' output$markers
#' 
#' @export
#' @import augere.core
#' @importFrom scrapper analyze.se
#' @importFrom scater plotXY
runSolo <- function(
    x, 
    rna.experiment = TRUE,
    adt.experiment = NULL,
    subset.factor = NULL, 
    subset.levels = NULL,
    block.field = NULL,
    qc.mito.seqnames = c("MT", "M", "chrM", "chrMT"),
    qc.mito.regex = NULL,
    qc.igg.regex = "IgG|igg|IGG",
    qc.num.mads = 3, 
    qc.filter = TRUE, 
    num.hvgs = 2000, 
    num.pcs = 25,
    mnn.num.neighbors = 15,
    mnn.num.steps = 1,
    cluster.method = 'graph',
    cluster.kmeans.k = 10,
    cluster.graph.method = c("multilevel", "leiden", "walktrap"),
    cluster.graph.num.neighbors = 10, 
    cluster.graph.resolution = NULL, 
    reduced.dimensions = c("tsne", "umap"),
    tsne.perplexity = 30,
    umap.num.neighbors = 15,
    umap.min.dist = 0.1,
    marker.effect.size = c("cohens.d", "auc", "delta.mean", "delta.detected"),
    marker.summary = c("min.rank", "mean", "median", "min"),
    marker.lfc.threshold = 0,
    assay = 1,
    symbol.field = NULL,
    metadata = NULL,
    output.dir = "solo", 
    author = NULL,
    dry.run = FALSE, 
    save.results = TRUE, 
    num.threads =1
) {
    restore.fun <- resetInputCache()
    on.exit(restore.fun(), after=FALSE, add=TRUE)

    if (is.null(author)) {
        author <- Sys.info()[["user"]]
    }
    author.txt <- deparseToString(as.list(author))

    dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
    fname <- file.path(output.dir, "report.Rmd")

    template <- system.file("templates", "solo.Rmd", package="augere.solo", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    ##################
    ### Setting up ###
    ##################

    .map_feature_type <- function(modality) {
        if (modality == "RNA" || modality == "rna") {
            "gene"
        } else {
            "tag"
        }
    }

    rna.obj <- adt.obj <- NULL
    main.modality <- NULL
    if (isTRUE(rna.experiment)) {
        main.modality <- "RNA"
        rna.obj <- "sce"
    } else if (isTRUE(adt.experiment)) {
        main.modality <- "ADT"
        adt.obj <- "sce"
    }
    if (is.null(main.modality)) {
        parsed[["mainexp-text"]] <- NULL
    } else {
        parsed[["mainexp-text"]] <- replacePlaceholders(
            parsed[["mainexp-text"]],
            list(
                MAIN_FEATURE_TYPE = .map_feature_type(main.modality[1]),
                MAIN_MODALITY = main.modality 

            )
        )
    }

    altexp.text <- list()
    altexp.names <- character()
    altexp.template <- "SingleCellExperiment::altExp(sce, %s)"

    if (is.character(rna.experiment)) {
        altexp.text[["RNA"]]  <- replacePlaceholders(
            parsed[["altexp-text"]],
            list(
                ALTEXP_SE = rna.experiment,
                ALTEXP_FEATURE_TYPE = .map_feature_type("RNA"),
                ALTEXP_MODALITY = "RNA"
            )
        )
        rna.obj <- sprintf(altexp.template, deparseToString(rna.experiment))
        altexp.names <- c(altexp.names, rna.experiment)
    }

    if (is.character(adt.experiment)) {
        altexp.text[["ADT"]]  <- replacePlaceholders(
            parsed[["altexp-text"]],
            list(
                ALTEXP_SE = adt.experiment,
                ALTEXP_FEATURE_TYPE = .map_feature_type("ADT"),
                ALTEXP_MODALITY = "ADT"
            )
        )
        adt.obj <- sprintf(altexp.template, deparseToString(adt.experiment))
        altexp.names <- c(altexp.names, adt.experiment)
    }

    use.rna <- !is.null(rna.obj)
    use.adt <- !is.null(adt.obj)

    parsed[["altexp-text"]] <- altexp.text
    parsed[["initialize-se"]] <- processInputCommands(x, name="se")

    if (!is.null(subset.factor)) {
        parsed[["subset-data"]] <- replacePlaceholders(
            parsed[["subset-data"]],
            list(
                SUBSET_FACTOR = deparseToString(subset.factor),
                SUBSET_LEVELS = deparseToString(subset.levels)
            )
        )
    } else {
        parsed[["subset-data"]] <- NULL
    }

    if (!is.null(block.field)) {
        parsed[["block-setup"]] <- replacePlaceholders(
            parsed[["block-setup"]],
            list(BLOCK = deparseToString(block.field))
        )
    } else {
        parsed[["block-setup"]] <- NULL
    }

    if (use.rna) {
        parsed[["rna-force-load"]] <- replacePlaceholders(
            parsed[["rna-force-load"]],
            list(
                RNA_SE = rna.obj,
                ASSAY = deparseToString(assay)
            )
        )
    } else {
        parsed[["rna-force-load"]] <- NULL
    }

    if (use.adt) {
        parsed[["adt-force-load"]] <- replacePlaceholders(
            parsed[["adt-force-load"]],
            list(
                ADT_SE = adt.obj,
                ASSAY = deparseToString(assay)
            )
        )
    } else {
        parsed[["adt-force-load"]] <- NULL
    }

    force.altexp <- list()
    for (ae in altexp.names) {
        force.altexp[[ae]] <- replacePlaceholders(
            parsed[["altexp-force-sce"]],
            list(ALTEXP = deparseToString(ae))
        )
    }
    parsed[["altexp-force-sce"]] <- force.altexp

    if (use.adt) {
        analysis.type <- "CITE-seq"
    } else {
        analysis.type <- "scRNA-seq"
    }
    parsed[[1]] <- replacePlaceholders(
        parsed[[1]],
        list(
            AUTHOR = author.txt,
            ANALYSIS_TYPE = analysis.type
        )
    )

    #######################
    ### Quality control ###
    #######################

    if (use.rna) {
        rna.qc.replacements <- list(
            RNA_SE = rna.obj,
            NUM_MADS = deparseToString(qc.num.mads),
            ASSAY = deparseToString(assay),
            NUM_THREADS = num.threads
        )

        rna.metrics.parsed <- parsed[["rna-qc-metrics"]]
        if (is.null(qc.mito.regex)) {
            rna.qc.replacements$MITO_SEQNAMES <- deparseToString(qc.mito.seqnames)
            rna.metrics.parsed[["define-mito-from-rownames"]] <- NULL
            rna.metrics.parsed[["define-mito-from-symbols"]] <- NULL
        } else {
            rna.metrics.parsed[["define-mito-from-ranges"]] <- NULL
            rna.qc.replacements$MITO_REGEX <- deparseToString(qc.mito.regex)
            if (is.null(symbol.field)) {
                rna.metrics.parsed[["define-mito-from-symbols"]] <- NULL
            } else {
                rna.metrics.parsed[["define-mito-from-rownames"]] <- NULL
                rna.qc.replacements$SYMBOL_FIELD <- deparseToString(symbol.field)
            }
        }

        if (!is.null(block.field)) {
            rna.qc.replacements$BLOCK <- deparseToString(block.field)
            rna.metrics.parsed[["plot-simple"]] <- NULL
        } else {
            rna.metrics.parsed[["plot-blocked"]] <- NULL
            rna.metrics.parsed[["block-args"]] <- NULL
            rna.metrics.parsed[["block-text"]] <- NULL
        }

        parsed[["rna-qc-metrics"]] <- replacePlaceholders(rna.metrics.parsed, rna.qc.replacements)
    } else {
        parsed[["rna-qc-metrics"]] <- NULL
    }

    if (use.adt) {
        adt.qc.replacements <- list(
            ADT_SE = adt.obj,
            IGG_REGEX = deparseToString(qc.igg.regex),
            NUM_MADS = deparseToString(qc.num.mads),
            ASSAY = deparseToString(assay),
            NUM_THREADS = num.threads
        )

        adt.metrics.parsed <- parsed[["adt-qc-metrics"]]
        if (is.null(symbol.field)) {
            adt.metrics.parsed[["define-igg-from-symbols"]] <- NULL
        } else {
            adt.metrics.parsed[["define-igg-from-rownames"]] <- NULL
            adt.qc.replacements$SYMBOL_FIELD <- deparseToString(symbol.field)
        }

        if (!is.null(block.field)) {
            adt.qc.replacements$BLOCK <- deparseToString(block.field)
            adt.metrics.parsed[["plot-simple"]] <- NULL
        } else {
            adt.metrics.parsed[["plot-blocked"]] <- NULL
            adt.metrics.parsed[["block-args"]] <- NULL
            adt.metrics.parsed[["block-text"]] <- NULL
        }

        parsed[["adt-qc-metrics"]] <- replacePlaceholders(adt.metrics.parsed, adt.qc.replacements)
    } else {
        parsed[["adt-qc-metrics"]] <- NULL
    }

    if (qc.filter) {
        filter.parsed <- parsed[["filter-qc"]]
        filter.replacements <- list()

        subset.by <- character()
        if (use.rna) {
            subset.by <- c(subset.by, "rna.qc.df$keep")
        }
        if (use.adt) {
            subset.by <- c(subset.by, "adt.qc.df$keep")
        }
        filter.replacements$QC_SUBSET <- paste(subset.by, collapse=" & ")

        if (!is.null(block.field)) {
            filter.replacements$BLOCK <- deparseToString(block.field)
        } else {
            filter.parsed[["with-block"]] <- NULL
        }

        parsed[["filter-qc"]] <- replacePlaceholders(filter.parsed, filter.replacements)
        parsed[["no-filter-qc"]] <- NULL
    } else {
        parsed[["filter-qc"]] <- NULL
    }

    #####################
    ### Normalization ###
    #####################

    if (is.null(block.field)) {
        parsed[["norm-block-text"]] <- NULL
    }

    if (use.rna) {
        norm.parsed <- parsed[["rna-norm"]]
        norm.replacements <- list(
            RNA_SE = rna.obj,
            ASSAY = deparseToString(assay)
        )

        if (!is.null(block.field)) {
            norm.replacements$BLOCK <- deparseToString(block.field)
            norm.replacements$PLOT_X <- deparseToString(block.field) 
        } else {
            norm.replacements$PLOT_X <- "NULL"
            norm.parsed[["block-args"]] <- NULL
        }

        parsed[["rna-norm"]] <- replacePlaceholders(norm.parsed, norm.replacements)
    } else {
        parsed[["rna-norm"]] <- NULL
    }

    if (use.adt) {
        norm.parsed <- parsed[["adt-norm"]]
        norm.replacements <- list(
            ADT_SE = adt.obj,
            ASSAY = deparseToString(assay),
            NUM_THREADS = num.threads
        )

        if (!is.null(block.field)) {
            norm.replacements$BLOCK <- deparseToString(block.field)
            norm.replacements$PLOT_X <- deparseToString(block.field) 
        } else {
            norm.replacements$PLOT_X <- "NULL"
            norm.parsed[["block-args"]] <- NULL
        }

        parsed[["adt-norm"]] <- replacePlaceholders(norm.parsed, norm.replacements)
    } else {
        parsed[["adt-norm"]] <- NULL
    }

    ##########################
    ### Variance modelling ###
    ##########################

    if (use.rna) {
        hvg.parsed <- parsed[["rna-hvg"]]
        hvg.replacements <- list(
            RNA_SE = rna.obj,
            TOP_HVGS = num.hvgs,
            NUM_THREADS = num.threads
        )

        if (is.null(symbol.field)) {
            hvg.replacements$HVG_SYMBOL_FIELD <- ""
        } else {
            hvg.replacements$HVG_SYMBOL_FIELD <- paste0(deparseToString(symbol.field), ", ")
        }

        if (is.null(block.field)) {
            hvg.parsed$blocked <- NULL
            hvg.parsed[["block-args"]] <- NULL
            hvg.parsed[["block-text"]] <- NULL
        } else {
            hvg.parsed$unblocked <- NULL
            hvg.replacements$BLOCK <- deparseToString(block.field)
        }

        parsed[["rna-hvg"]] <- replacePlaceholders(hvg.parsed, hvg.replacements)
    } else {
        parsed[["rna-hvg"]] <- NULL
    }

    ###########
    ### PCA ###
    ###########

    if (is.null(block.field)) {
        parsed[["pca-block-text"]] <- NULL
    }

    if (use.rna) {
        pca.parsed <- parsed[["rna-pca"]]
        pca.replacements <- list(
            RNA_SE = rna.obj,
            TOP_PCS = deparseToString(num.pcs),
            NUM_THREADS = num.threads
        )

        if (is.null(block.field)) {
            pca.parsed[["block-args"]] <- NULL
        } else {
            pca.replacements$BLOCK <- deparseToString(block.field)
        }

        parsed[["rna-pca"]] <- replacePlaceholders(pca.parsed, pca.replacements)
    } else {
        parsed[["rna-pca"]] <- NULL
    }

    if (use.adt) {
        pca.parsed <- parsed[["adt-pca"]]
        pca.replacements <- list(
            ADT_SE = adt.obj,
            TOP_PCS = deparseToString(num.pcs),
            NUM_THREADS = num.threads
        )

        if (is.null(block.field)) {
            pca.parsed[["block-args"]] <- NULL
            pca.parsed[["block-text"]] <- NULL
        } else {
            pca.replacements$BLOCK <- deparseToString(block.field)
        }

        parsed[["adt-pca"]] <- replacePlaceholders(pca.parsed, pca.replacements)
    } else {
        parsed[["adt-pca"]] <- NULL
    }

    if (use.adt && use.rna) {
        chosen.pcs <- "combined"
        if (is.null(main.modality)) {
            main.name <- NULL
        } else {
            main.name <- "PCA" 
        }

        combine.parsed <- parsed[["combined-pca"]]
        combine.replacements <- list(
            MAIN_NAME = deparseToString(main.name),
            ALTEXP_NAMES = deparseToString(altexp.names),
            NUM_THREADS = num.threads
        )

        if (is.null(block.field)) {
            combine.parsed[["block-args"]] <- NULL
            combine.parsed[["block-text"]] <- NULL
        } else {
            combine.replacements$BLOCK <- deparseToString(block.field)
        }

        parsed[["combined-pca"]] <- replacePlaceholders(combine.parsed, combine.replacements)
    } else {
        parsed[["combined-pca"]] <- NULL
        chosen.pcs <- "PCA"
        if (length(altexp.names) == 1L) {
            names(chosen.pcs) <- altexp.names
        }
    }

    if (!is.null(block.field)) {
        parsed[["mnn"]] <- replacePlaceholders(
            parsed[["mnn"]],
            list(
                BLOCK = deparseToString(block.field),
                MNN_NUM_NEIGHBORS = mnn.num.neighbors,
                MNN_NUM_STEPS = mnn.num.steps,
                NUM_THREADS = num.threads,
                CHOSEN_PCS = deparseToString(chosen.pcs)
            )
        )
        chosen.pcs <- "MNN"
    } else {
        parsed[["mnn"]] <- NULL
    }

    ##########################
    ### k-means clustering ###
    ##########################

    chosen.cluster <- NULL
    if (length(cluster.method)) {
        cluster.method <- match.arg(cluster.method, c("kmeans", "graph"), several.ok=TRUE)
        if (length(cluster.method) > 0) {
            if (cluster.method[1] == "kmeans") {
                chosen.cluster <- "kmeans.cluster"
            } else {
                chosen.cluster <- "graph.cluster"
            }
        }
    }

    if ("kmeans" %in% cluster.method) {
        kmeans.parsed <- parsed[["cluster-kmeans"]]
        kmeans.replacements <- list(
            CLUSTER_KMEANS_K = cluster.kmeans.k,
            NUM_THREADS = num.threads,
            CHOSEN_PCS = deparseToString(chosen.pcs)
        )

        if (is.null(block.field)) {
            kmeans.parsed[["block"]] <- NULL
        } else {
            kmeans.replacements$BLOCK <- deparseToString(block.field)
        }

        parsed[["cluster-kmeans"]] <- replacePlaceholders(kmeans.parsed, kmeans.replacements)
    } else {
        parsed[["cluster-kmeans"]] <- NULL 
    }

    ##########################
    ### Neighbor detection ###
    ##########################

    parsed.nn <- parsed[["neighbor-steps"]]
    nn.replacements <- list(
        CHOSEN_PCS = deparseToString(chosen.pcs),
        NUM_THREADS = num.threads
    )

    if ("graph" %in% cluster.method) {
        cluster.graph.method <- match.arg(cluster.graph.method)
        cluster.args <- list(method = cluster.graph.method)
        if (cluster.graph.method == "leiden") {
            cluster.args$leiden.resolution <- cluster.graph.resolution
        } else if (cluster.graph.method == "multilevel") {
            cluster.args$multilevel.resolution <- cluster.graph.resolution
        }

        nn.replacements$BUILD_SNN_GRAPH_ARGS <- deparseToString(list(num.neighbors = cluster.graph.num.neighbors))
        nn.replacements$CLUSTER_GRAPH_ARGS <- deparseToString(cluster.args)
        graph.in.use <- TRUE

        graph.parsed <- parsed.nn[["cluster-graph"]]
        if (is.null(block.field)) {
            graph.parsed[["block"]] <- NULL
        } else {
            nn.replacements$BLOCK <- deparseToString(block.field)
        }
        parsed.nn[["cluster-graph"]] <- graph.parsed
    } else {
        parsed.nn[["cluster-graph"]] <- NULL 
        nn.replacements$BUILD_SNN_GRAPH_ARGS <- "NULL"
        nn.replacements$CLUSTER_GRAPH_ARGS <- "NULL"
        graph.in.use <- FALSE 
    }

    parsed.plot <- parsed.nn[["plot-any"]]
    if (is.null(block.field)) {
        if (!is.null(chosen.cluster)) {
            colour.by.choice <- sprintf(", colour_by = %s", deparseToString(chosen.cluster))
        } else {
            colour.by.choice <- ""
        }
    } else {
        if (is.null(chosen.cluster)) {
            colour.by.choice <- sprintf(", colour_by = %s", deparseToString(block.field))
        } else {
            colour.by.choice <- NA
        }
    }

    if (length(reduced.dimensions)) {
        reduced.dimensions <- match.arg(reduced.dimensions, several.ok=TRUE)
        plot.in.use <- TRUE
    } else {
        plot.in.use <- FALSE
    }

    if ("tsne" %in% reduced.dimensions) {
        nn.replacements$TSNE_ARGS <- deparseToString(list(perplexity = tsne.perplexity))

        if (is.null(block.field) || is.null(chosen.cluster)) {
            parsed.plot[["tsne-simple"]] <- replacePlaceholders(
                parsed.plot[["tsne-simple"]],
                list(COLOUR_BY_CHOICE = colour.by.choice)
            )
            parsed.plot[["tsne-multi"]] <- NULL
        } else {
            parsed.plot[["tsne-multi"]] <- replacePlaceholders(
                parsed.plot[["tsne-multi"]],
                list(
                    CHOSEN_CLUSTER = deparseToString(chosen.cluster),
                    BLOCK = deparseToString(block.field)
                )
            )
            parsed.plot[["tsne-simple"]] <- NULL
        }

    } else {
        parsed.plot[["tsne-simple"]] <- NULL
        parsed.plot[["tsne-multi"]] <- NULL
        nn.replacements$TSNE_ARGS <- "NULL"
    }

    if ("umap" %in% reduced.dimensions) {
        nn.replacements$UMAP_ARGS <- deparseToString(list(num.neighbors = umap.num.neighbors, min.dist = umap.min.dist))

        if (is.null(block.field) || is.null(chosen.cluster)) {
            parsed.plot[["umap-simple"]] <- replacePlaceholders(
                parsed.plot[["umap-simple"]],
                list(COLOUR_BY_CHOICE = colour.by.choice)
            )
            parsed.plot[["umap-multi"]] <- NULL
        } else {
            parsed.plot[["umap-multi"]] <- replacePlaceholders(
                parsed.plot[["umap-multi"]],
                list(
                    CHOSEN_CLUSTER = deparseToString(chosen.cluster),
                    BLOCK = deparseToString(block.field)
                )
            )
            parsed.plot[["umap-simple"]] <- NULL
        }

    } else {
        parsed.plot[["umap-simple"]] <- NULL
        parsed.plot[["umap-multi"]] <- NULL
        nn.replacements$UMAP_ARGS <- "NULL"
    }

    if (plot.in.use) {
        parsed.nn[["plot-any"]] <- parsed.plot
    } else {
        parsed.nn[["plot-any"]] <- NULL
    }

    if (plot.in.use || graph.in.use) {
        parsed[["neighbor-steps"]] <- replacePlaceholders(parsed.nn, nn.replacements)
    } else {
        parsed[["neighbor-steps"]] <- NULL
    }

    ###############
    ### Markers ###
    ###############

    if (is.null(chosen.cluster)) {
        parsed[["markers"]] <- NULL

    } else {
        marker.summary <- match.arg(marker.summary)
        marker.effect.size <- match.arg(marker.effect.size)
        marker.field <- paste0(marker.effect.size, ".", marker.summary)

        marker.config <- list()
        if (use.rna) {
            marker.config$RNA <- list(lower="rna", feature=.map_feature_type("RNA"), se=rna.obj)
        }
        if (use.adt) {
            marker.config$ADT <- list(lower="adt", feature=.map_feature_type("ADT"), se=adt.obj)
        }

        if (marker.lfc.threshold) {
            more.args <- sprintf("\n    more.marker.args = list(threshold = %s),", marker.lfc.threshold)
        } else {
            more.args <- ""
        }

        parsed.markers <- parsed[["markers"]]
        parsed.markers[["overview"]] <- replacePlaceholders(
            parsed.markers[["overview"]],
            list(CHOSEN_CLUSTER = deparseToString(chosen.cluster))
        )
        if (is.null(block.field)) {
            parsed.markers[["overview"]][["block-text"]] <- NULL
        }

        all.markers <- list()
        for (mod in names(marker.config)) {
            payload <- marker.config[[mod]]
            current <- parsed.markers[["modality"]]

            if (is.null(block.field)) {
                current[["block-args"]] <- NULL
            } else {
                current[["block-args"]] <- replacePlaceholders(
                    current[["block-args"]],
                    list(BLOCK = deparseToString(block.field))
                )
            }

            if (is.null(symbol.field)) {
                current[["show-marker-symbols"]] <- NULL
                current[["extra-symbol-column"]] <- NULL
            } else {
                current[["extra-symbol-column"]] <- replacePlaceholders(
                    current[["extra-symbol-column"]],
                    list(SYMBOL_FIELD = deparseToString(symbol.field))
                )
                current[["show-marker-symbols"]] <- replacePlaceholders(
                    current[["show-marker-symbols"]],
                    list(
                        MARKER_PREFIX = payload$lower,
                        SYMBOL_FIELD = deparseToString(symbol.field)
                    )
                )
                current[["show-marker-ids"]] <- NULL
            }

            all.markers[[mod]] <- replacePlaceholders(
                current,
                list(
                    MARKER_MODALITY = mod,
                    CHOSEN_CLUSTER = deparseToString(chosen.cluster),
                    MARKER_FEATURE_TYPE = payload$feature,
                    MARKER_SE = payload$se,
                    MARKER_PREFIX = payload$lower,
                    MARKER_ORDER = deparseToString(marker.field),
                    MARKER_ARGS = more.args,
                    NUM_THREADS = num.threads
                )
            )
        }

        parsed.markers[["modality"]] <- all.markers
        parsed$markers <- parsed.markers
    }

    ##############
    ### Saving ###
    ##############

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        parsed[["create-common-metadata"]] <- replacePlaceholders(
            parsed[["create-common-metadata"]],
            list(COMMON_METADATA = deparseToString(metadata))
        )
    } else {
        parsed[["create-common-metadata"]] <- NULL
    }

    parsed[["save-sce"]] <- replacePlaceholders(
        parsed[["save-sce"]],
        list(
            AUTHOR = author.txt,
            ANALYSIS_TYPE = analysis.type
        )
    )
    if (merge.metadata) {
        parsed[["save-sce"]][["no-merge-metadata"]] <- NULL
    } else {
        parsed[["save-sce"]][["merge-metadata"]] <- NULL
    }

    modalities.to.save <- c("RNA", "ADT")[c(use.rna, use.adt)]

    qc.save.code <- list()
    for (nm in modalities.to.save) {
        current <- replacePlaceholders(
            parsed[["save-qc-body"]],
            list(
                AUTHOR = author.txt,
                MODALITY_NAME = nm,
                MODALITY_LOWER = tolower(nm)
            )
        )
        if (merge.metadata) {
            current[["no-merge-metadata"]] <- NULL
        } else {
            current[["merge-metadata"]] <- NULL
        }
        qc.save.code[[nm]] <- current
    }
    parsed[["save-qc-body"]] <- qc.save.code

    if (is.null(chosen.cluster)) {
        parsed[["save-markers"]] <- NULL
    } else {
        marker.save.code <- list()
        for (nm in modalities.to.save) {
            current <- replacePlaceholders(
                parsed[["save-markers"]][["body"]],
                list(
                    AUTHOR = author.txt,
                    MARKER_ORDER = deparseToString(marker.field),
                    CHOSEN_CLUSTER = deparseToString(chosen.cluster),
                    MODALITY_NAME = nm,
                    MODALITY_FEATURE_TYPE = .map_feature_type(nm),
                    MODALITY_LOWER = tolower(nm)
                )
            )
            if (merge.metadata) {
                current[["no-merge-metadata"]] <- NULL
            } else {
                current[["merge-metadata"]] <- NULL
            }
            marker.save.code[[nm]] <- current
        }
        parsed[["save-markers"]][["body"]] <- marker.save.code
    }

    writeRmd(parsed, file=fname)
    if (dry.run) {
        return(NULL)
    }

    if (save.results) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- c("save-directory", "save-sce", "save-qc", "save-markers")
    }

    # Avoid spewing out cat() statements in the marker detection section.
    parsed$markers <- sub("^ *cat\\(", "list\\(", unlist(parsed$markers))

    env <- new.env()
    compileReport(fname, env=env, skip.chunks=skip.chunks, contents=parsed)

    output <- list(sce=env$sce)
    if (use.rna) {
        output$qc.rna <- env$rna.qc.df
        output$markers.rna <- env$rna.markers
    }
    if (use.adt) {
        output$qc.adt <- env$adt.qc.df
        output$markers.adt <- env$adt.markers
    }

    output
}
