#' Cell type annotation from scRNA-seq data
#'
#' Annotate cells in a scRNA-seq dataset by computing correlations against reference data with known labels.
#'
#' @param test A \link[SummarizedExperiment]{SummarizedExperiment} object containing cells in the test dataset to be assigned labels.
#' @param references A list created by \code{configureReferenceAnnotation}, containing a configuration for a reference dataset.
#' 
#' Alternatively, a list of these lists may be supplied to specify multiple references.
#' This list may be named, in which case the names will be used to identify each reference in both the report and output of this function.
#' @param test.id.field String specifying the name of the \code{\link[SummarizedExperiment]{rowData}(test)} column containing the common gene identifiers.
#' These identifiers should be consistent with those that are used in each entry of \code{references}.
#' If \code{NULL}, the row names are assumed to contain the relevant identifiers.
#' @param test.symbol.field String specifying the name of the \code{\link[SummarizedExperiment]{rowData}(test)} column that contains gene symbols.
#' This used to improve the interpretability of some of the gene-based diagnostics in the report.
#' If \code{NULL}, no symbols are added.
#' @param test.is.ensembl Boolean indicating whether \code{rownames(test)} contains Ensembl IDs.
#' Otherwise, it is assumed to contain gene symbols.
#' Only relevant for entries of \code{references} where \code{ref} is a \pkg{celldex} reference,
#' where it is used to choose between Ensembl IDs or gene symbols for the row names. 
#' @param test.assay Integer or string specifying the assay of \code{test} containing the expression values to use for classification.
#' This is expected to be raw counts or log-transformed normalized expression values.
#' @param test.is.lognorm Boolean indicating whether the assay at \code{test.assay} contains log-normalized values.
#' If \code{FALSE}, it is assumed to contain counts. 
#' @param test.block.field String specifying the name of the \code{\link[SummarizedExperiment]{colData}(test)} column containing the block assignment for each cell.
#' This is used to ensure that some marker-related diagnostics are not driven by uninteresting differences between blocks.
#' If \code{NULL}, all cells in \code{test} are assumed to originate from the same block.
#' @param cluster.field String specifying the column of \code{\link[SummarizedExperiment]{colData}(test)} with cluster assignments for each cell.
#' If provided, this is used to generate some diagnostics to compare the empirical clusters to the predicted labels.
#' @param reduced.dimensions Character vector of reduced dimensions in \code{\link[SingleCellExperiment]{reducedDims}(test)} to plot with the predicted labels in the report.
#' @param output.dir String containing the path to an output directory in which to write the Rmarkdown report and save results.
#' @param metadata Named list of additional fields to add to each result's metadata.
#' @param author Character vector of authors.
#' @param dry.run Boolean indicating whether to perform a dry run.
#' This will write the Rmarkdown report without evaluating it.
#' @param save.results Boolean indicating whether the results should also be saved to file.
#' @param num.threads Integer specifying the number of threads to use in the various computations.
#' @param ref A \link[SummarizedExperiment]{SummarizedExperiment} object containing reference samples with known labels in \code{colData(ref)[[label.field]]}.
#' 
#' Alternatively, a string can be provided containing the name of \pkg{celldex} reference dataset (e.g., \code{"HumanPrimaryCellAtlasData"}).
#' In such cases, \code{label.field} is typically \code{"label.main"} or \code{"label.fine"}.
#' @param ref.label.field String specifying the name of the \code{\link[SummarizedExperiment]{colData}(ref)} column containing the label for each reference sample.
#' For \pkg{celldex} references, this is usually either \code{"label.main"} or \code{"label.fine"}.
#' @param ref.marker.method String specifying the method for choosing the top markers from each pairwise comparison between labels,
#' see the \code{de.method=} argument in \code{\link[SingleR]{trainSingleR}} for more details.
#' If \code{NULL}, this defaults to \code{"classic"} for the \pkg{celldex} references and \code{"t"} otherwise.
#' @param ref.num.markers Integer specifying the number of markers to use from each pairwise comparison between labels.
#' See the \code{de.n=} argument in \code{\link[SingleR]{trainSingleR}} for more details.
#' @param ref.assay Integer or string specifying the assay of \code{ref} containing the expression values to use for creating references.
#' This should be a matrix-like object that contains counts or their log-normalized values.
#' @param ref.is.lognorm Boolean indicating whether the assay at \code{assay} contains log-normalized values.
#' If \code{FALSE}, it is assumed to contain counts.
#' @param ref.id.field String specifying the name of the \code{\link[SummarizedExperiment]{rowData}(ref)} column containing the common gene identifiers.
#' These identifiers should be consistent with those that are used in \code{test}.
#' If \code{NULL}, the row names are assumed to contain the relevant identifiers.
#' @param ref.block.field String specifying the name of the \code{\link[SummarizedExperiment]{colData}(ref)} column containing the block assignment for each sample.
#' This is used to ensure that marker calculations are not driven by uninteresting differences between blocks.
#' If \code{NULL}, all cells in \code{test} are assumed to originate from the same block.
#' @param ref.aggregate Boolean indicating that references should be aggregated inside \code{\link[SingleR]{trainSingleR}}.
#' This can be set to \code{TRUE} for faster classification of large single-cell references.
#' @param suppress.plots Boolean indicating whether to suppress the plots.
#'
#' @return
#' For \code{runAnnotate}, a Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir} that contains the analysis commands.
#'
#' If \code{dry.run=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{predictions}, a list of length equal to \code{references}.
#' Each entry is a \link[S4Vectors]{DataFrame} containing the classification results of \code{test} against the corresponding reference.
#' See \code{\link[SingleR]{classifySingleR}} for more details on the expected columns.
#' \item \code{combined}, a \link[S4Vectors]{DataFrame} containing the combined classification results across multiple references.
#' See \code{\link[SingleR]{combineRecomputedResults}} for more details on the expected format.
#' Only present if \code{references} contains more than one reference.
#' }
#'
#' If \code{save.results=TRUE}, the results are saved in a \code{results} directory inside \code{output}.
#'
#' If \code{dry.run=TRUE}, \code{NULL} is returned.
#' Only the Rmarkdown report is saved to file.
#'
#' For \code{configureReferenceAnnotation}, a list of class \code{"reference"} is returned containing the configuration details for each reference.
#'
#' @examples
#' library(scRNAseq)
#' hESCs <- LaMannoBrainData('human-es')
#' hESCs <- hESCs[,1:100] # subsetting to speed it up.
#'
#' tmp <- tempfile()
#' results <- runAnnotate(
#'     hESCs,
#'     configureReferenceAnnotation(
#'         "HumanPrimaryCellAtlasData",
#'         "label.main"
#'     ),
#'     output.dir = tmp,
#'     num.threads = 2 # speed it up a little.
#' )
#'
#' list.files(tmp, recursive=TRUE)
#' results$predictions
#'
#' @export
#' @import augere.core
runAnnotate <- function(
    test,
    references,
    test.assay = 1,
    test.id.field = NULL,
    test.block.field = NULL,
    test.is.lognorm = (test.assay == "logcounts"),
    test.is.ensembl = FALSE, 
    test.symbol.field = NULL,
    cluster.field = NULL,
    reduced.dimensions = NULL,
    output.dir = "annotate", 
    metadata = NULL,
    author = NULL,
    dry.run = FALSE, 
    save.results = TRUE, 
    suppress.plots = interactive(),
    num.threads = 1
) {
    restore.fun <- resetInputCache()
    on.exit(restore.fun(), after=FALSE, add=TRUE)

    if (is.null(author)) {
        author <- Sys.info()[["user"]]
    }
    author.txt <- deparseToString(as.list(author))

    dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
    fname <- file.path(output.dir, "report.Rmd")

    template <- system.file("templates", "annotate.Rmd", package="augere.solo", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    ###################
    ### Set up test ###
    ###################

    parsed[[1]] <- replacePlaceholders(
        parsed[[1]],
        list(AUTHOR = author.txt)
    )

    parsed[["test-setup"]] <- processInputCommands(test, name="test")

    if (is.null(test.id.field)) {
        parsed[["replace-test-rownames"]] <- NULL
    } else {
        parsed[["replace-test-rownames"]] <- replacePlaceholders(
            parsed[["replace-test-rownames"]],
            list(ALTERNATIVE_NAMES = deparseToString(test.id.field))
        )
    }

    if (test.is.lognorm) {
        parsed[["normalize-test"]] <- NULL
        test.norm.assay <- test.assay
    } else {
        norm.parsed <- parsed[["normalize-test"]]
        if (is.null(test.block.field)) {
            norm.parsed[["block"]] <- NULL
        }
        parsed[["normalize-test"]] <- replacePlaceholders(
            norm.parsed,
            list(
                ASSAY = deparseToString(test.assay),
                BLOCK = deparseToString(test.block.field),
                NUM_THREADS = deparseToString(num.threads)
            )
        )
        test.norm.assay <- "logcounts"
    }

    ##########################
    ### Process references ###
    ##########################

    if (inherits(references, "reference")) {
        references <- list(references)
    }
    if (length(references) < 1) {
        stop("expected at least one reference entry in 'references'")
    }
    multiple.references <- length(references) > 1

    all.ref.chunks <- list()
    all.anno.chunks <- list()
    all.save.chunks <- list()
    all.fig.chunks <- character(0)

    for (r in seq_along(references)) {
        curref <- references[[r]]
        short.name <- names(references)[r]
        if (is.null(short.name)) {
            short.name <- as.character(r)
        }

        ###############
        ### Loading ###
        ###############

        ref.parsed <- parsed[["ref-body"]]
        if (is.character(curref$ref)) {
            celldex.args <- list(CELLDEX_CMD = curref$ref)
            if (test.is.ensembl) {
                celldex.args$CELLDEX_ARGS <- "ensembl=TRUE"
            } else {
                celldex.args$CELLDEX_ARGS <- ""
            }
            ref.parsed[["celldex"]] <- replacePlaceholders(ref.parsed[["celldex"]], celldex.args)
            ref.parsed[["custom"]] <- NULL
        } else {
            ref.parsed[["custom"]][["ref-setup"]] <- processInputCommands(curref$ref, name="ref")
            ref.parsed[["celldex"]] <- NULL
        }

        if (curref$is.lognorm) {
            ref.parsed[["normalize"]] <- NULL
        } else {
            norm.parsed <- ref.parsed[["normalize"]]
            if (is.null(curref$block.field)) {
                norm.parsed[["block"]] <- NULL
            }
            ref.parsed[["normalize"]] <- replacePlaceholders(
                norm.parsed,
                list(
                    ASSAY = deparseToString(curref$assay),
                    BLOCK = deparseToString(curref$block.field),
                    NUM_THREADS = deparseToString(num.threads)
                )
            )
            curref$assay <- "logcounts"
        }

        if (is.null(curref$id.field)) {
            ref.parsed[["replace-rownames"]] <- NULL
        } else {
            ref.parsed[["replace-rownames"]] <- replacePlaceholders(
                ref.parsed[["replace-rownames"]],
                list(ALTERNATIVE_NAMES = deparseToString(curref$id.field))
            )
        }

        all.ref.chunks[[r]] <- replacePlaceholders(
            ref.parsed,
            list(
                LABEL_FIELD = deparseToString(curref$label.field),
                SHORT_NAME = deparseToString(short.name),
                SHORT_NAME_RAW = short.name
           )
        )

        ##################
        ### Annotation ###
        ##################

        de.method <- curref$marker.method
        if (is.null(de.method)) {
            if (is.character(curref$ref)) {
                de.method <- "classic"
            } else {
                de.method <- "t"
            }
        }

        extra.train.args <- character(0)
        if (!is.null(curref$num.markers)) {
            extra.train.args <- c(extra.train.args, paste("de.n =", curref$num.markers)) 
        }
        if (multiple.references) {
            extra.train.args <- c(extra.train.args, "restrict = universe")
        }
        if (!is.null(curref$aggregate)) {
            extra.train.args <- c(extra.train.args, "aggr.ref = TRUE")
        }
        if (!is.null(curref$block.field)) {
            extra.train.args <- c(extra.train.args, sprintf("de.args = list(block = SummarizedExperiment::colData(ref)[,%s])", deparseToString(curref$block.field)))
        }

        if (length(extra.train.args)) {
            extra.train.args <- paste(sprintf("\n    %s,", extra.train.args), collapse="")
        } else {
            extra.train.args <- ""
        }

        anno.parsed <- parsed[["annotate-body"]]

        if (!is.null(test.symbol.field)) {
            anno.parsed[["alt-row-names"]] <- replacePlaceholders(
                anno.parsed[["alt-row-names"]],
                list(DISPLAY_NAME_FIELD = deparseToString(test.symbol.field))
            )
            anno.parsed[["old-row-names"]] <- NULL 
        } else {
            anno.parsed[["alt-row-names"]] <- NULL 
            if (is.null(test.id.field)) {
                # If the test ID field was specified, the row names have been changed, so we need to retrieve the old row names.
                anno.parsed[["old-row-names"]] <- NULL 
            }
        }

        if (!is.null(test.block.field)) {
            anno.parsed[["block"]] <- replacePlaceholders(
                anno.parsed[["block"]],
                list(BLOCK = deparseToString(test.block.field))
            )
        } else {
            anno.parsed[["block"]] <- NULL 
        }

        all.anno.chunks[[r]] <- replacePlaceholders(
            anno.parsed,
            list(
                NUM_THREADS = num.threads,
                SHORT_NAME = deparseToString(short.name),
                SHORT_NAME_RAW = short.name,
                LABEL_FIELD = deparseToString(curref$label.field),
                DE_METHOD = deparseToString(de.method),
                TRAIN_ARGS = extra.train.args,
                ASSAY_NORM_REF = deparseToString(curref$assay),
                ASSAY_TEST = deparseToString(test.assay),
                ASSAY_NORM_TEST = deparseToString(test.norm.assay)
            )
        )

        all.fig.chunks <- c(
            all.fig.chunks,
            paste0("plot-scores-", short.name),
            paste0("plot-markers-", short.name)
        )

        ##############
        ### Saving ###
        ##############

        save.parsed <- parsed[["save-body"]]
        if (is.null(metadata)) {
            save.parsed[["merge-metadata"]] <- NULL
        } else {
            save.parsed[["no-merge-metadata"]] <- NULL
        }

        all.save.chunks[[r]] <- replacePlaceholders(
            save.parsed,
            list(
                AUTHOR = author.txt,
                SHORT_NAME_RAW = short.name,
                SHORT_NAME = deparseToString(short.name),
                SAVE_PATH = deparseToString(paste0("pred-", short.name))
            )
        )
    }

    parsed[["ref-body"]] <- all.ref.chunks
    parsed[["annotate-body"]] <- all.anno.chunks
    parsed[["save-body"]] <- all.save.chunks

    ###########################
    ### Combined annotation ###
    ###########################

    if (multiple.references) {
        parsed[["combining"]] <- replacePlaceholders(
            parsed[["combining"]],
            list(
                ASSAY_TEST = deparseToString(test.assay),
                NUM_THREADS = num.threads 
            )
        )

        save.parsed <- parsed[["save-combined"]]
        if (is.null(metadata)) {
            save.parsed[["merge-metadata"]] <- NULL
        } else {
            save.parsed[["no-merge-metadata"]] <- NULL
        }
        parsed[["save-combined"]] <- replacePlaceholders(
            save.parsed,
            list(AUTHOR = author.txt)
        )

    } else {
        parsed[["intersect-ref"]] <- NULL
        parsed[["combining"]] <- NULL
        parsed[["save-combined"]] <- NULL
    }

    ###########################
    ### More visualizations ###
    ###########################

    if (is.null(reduced.dimensions) && is.null(cluster.field)) {
        parsed[["more-comparisons"]] <- NULL
    } else {
        if (multiple.references){ 
            target <- "combined"
        } else {
            target <- "predictions[[1]]"
        }

        more.parsed <- parsed[["more-comparisons"]]
        if (is.null(cluster.field)) {
            more.parsed[["clustering"]] <- NULL
        } else {
            more.parsed[["clustering"]] <- replacePlaceholders(
                more.parsed[["clustering"]],
                list(
                    TARGET = target,
                    CLUSTER_FIELD = deparseToString(cluster.field)
                )
            )
        }

        if (is.null(reduced.dimensions)) {
            more.parsed[["visualization"]] <- NULL
        } else {
            more.parsed[["visualization"]] <- replacePlaceholders(
                more.parsed[["visualization"]],
                list(
                    TARGET = target,
                    REDUCED_DIMENSIONS = deparseToString(reduced.dimensions)
                )
            )
        }

        parsed[["more-comparisons"]] <- more.parsed
    }

    ###################
    ### Wrapping up ###
    ###################

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        parsed[["create-common-metadata"]] <- replacePlaceholders(
            parsed[["create-common-metadata"]],
            list(COMMON_METADATA = deparseToString(metadata))
        )
    } else {
        parsed[["create-common-metadata"]] <- NULL
    }

    writeRmd(parsed, file=fname)
    if (dry.run) {
        return(NULL)
    }

    if (save.results) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- c("save-directory", "save-pred", "save-combined")
    }

    # Avoid spewing out cat() statements in the marker detection section.
    parsed <- sub("^ *cat\\(", "list\\(", unlist(parsed, use.names=FALSE))

    # No point even running the plots if we're going to skip them.
    if (suppress.plots) {
        skip.chunks <- c(
            skip.chunks,
            all.fig.chunks,
            "plot-reddim"
        )
    }

    env <- new.env()
    compileReport(fname, env=env, skip.chunks=skip.chunks, contents=parsed, suppress.plots=suppress.plots)

    output <- list(predictions = env$predictions)
    if (!is.null(env$combined)) {
        output$combined <- env$combined
    }

    output
}

#' @export
#' @rdname runAnnotate
configureReferenceAnnotation <- function(
    ref, 
    ref.label.field, 
    ref.marker.method = NULL, 
    ref.num.markers = NULL,
    ref.assay = "logcounts",
    ref.id.field = NULL,
    ref.block.field = NULL,
    ref.aggregate = NULL,
    ref.is.lognorm = (ref.assay == "logcounts")
) {
    output <- list(
        ref = ref,
        label.field = ref.label.field,
        marker.method = ref.marker.method,
        num.markers = ref.num.markers,
        assay = ref.assay,
        id.field = ref.id.field,
        aggregate = ref.aggregate,
        is.lognorm = ref.is.lognorm,
        block.field = ref.block.field
    )
    class(output) <- "reference"
    output
}
