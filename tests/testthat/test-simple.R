# library(testthat); library(augere.solo); source("test-simple.R")

set.seed(9999)
ngenes <- 10000L
ncells <- 200L
mu <- 2^rnorm(10000)
y <- matrix(rpois(ngenes * ncells, lambda=mu), ncol=ncells)
se <- SummarizedExperiment::SummarizedExperiment(list(counts=y))
is.mito <- rbinom(ngenes, 1, 0.01) == 1
SummarizedExperiment::rowRanges(se) <- GenomicRanges::GRanges(c("chrA", "chrM")[is.mito + 1], IRanges::IRanges(1, 1))
rownames(se) <- sprintf("GENE-%s", seq_len(ngenes))

tmp <- tempfile()
pdf(file=NULL)
ref <- runSolo(se, output.dir=tmp, suppress.plots=FALSE)
dev.off()

test_that("runSolo works correctly by default", {
    expect_s4_class(ref$sce, "SingleCellExperiment")
    expect_identical(nrow(ref$sce), ngenes)
    expect_lt(ncol(ref$sce), ncells)

    expect_identical(SummarizedExperiment::assayNames(ref$sce), c("counts", "logcounts"))
    expect_identical(SingleCellExperiment::reducedDimNames(ref$sce), c("PCA", "TSNE", "UMAP"))
    expect_true(is.factor(ref$sce$graph.cluster))

    is.hvg <- SummarizedExperiment::rowData(ref$sce)$hvg
    expect_true(any(is.hvg))
    expect_false(all(is.hvg))
    expect_identical(rownames(S4Vectors::metadata(ref$sce)$PCA$rotation), rownames(ref$sce)[is.hvg])

    expect_s4_class(ref$qc.rna, "DFrame")
    expect_identical(nrow(ref$qc.rna), ncells)
    expect_gt(mean(ref$qc.rna$subset.proportion.mito), 0)
    expect_true(any(!ref$qc.rna$keep))

    expect_s4_class(SingleCellExperiment::counts(ref$sce), "DelayedMatrix")
    expect_s4_class(SingleCellExperiment::logcounts(ref$sce), "LogNormalizedMatrix")

    expect_identical(names(ref$markers.rna), levels(ref$sce$graph.cluster))
    expect_s4_class(ref$markers.rna[[1]], "DFrame")
    expect_identical(sort(rownames(ref$markers.rna[[1]])), sort(rownames(se)))

    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "sce"))$x, "SingleCellExperiment")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "qc-rna"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "markers-rna", "1"))$x, "DFrame")

    lines <- readLines(file.path(tmp, "report.Rmd"))
    expect_false(any(grepl("block", lines))) # there shouldn't be any references to blocking here.
    expect_false(any(grepl("ADT", lines, ignore.case=TRUE))) # there shouldn't be any references to ADTs here.
})

test_that("runSolo works with a subset", {
    set.seed(999991)
    copy <- se
    copy$blah <- sample(LETTERS, ncol(se), replace=TRUE)
    chosen <- LETTERS[1:20]

    tmp <- tempfile()
    output <- runSolo(copy, subset.factor="blah", subset.levels=chosen, output.dir=tmp, reduced.dimensions=NULL, suppress.plots=TRUE, save.results=FALSE)
    ref <- runSolo(copy[,copy$blah %in% chosen], output.dir=tmp, reduced.dimensions=NULL, suppress.plots=TRUE, save.results=FALSE)
    expect_equal(output, ref)
})

test_that("runSolo works without QC filtering", {
    tmp <- tempfile()
    output <- runSolo(se, qc.filter=FALSE, reduced.dimensions=character(0), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_s4_class(output$sce, "SingleCellExperiment")
    expect_identical(ncol(output$sce), ncells)
    expect_identical(output$sce$sum, output$qc.rna$sum)
    expect_identical(output$qc.rna, ref$qc.rna)
})

test_that("runSolo works with a GRL for mitochondrial identification", {
    copy <- se
    SummarizedExperiment::rowRanges(copy) <- as(SummarizedExperiment::rowRanges(se), "GRangesList")

    tmp <- tempfile()
    output <- runSolo(copy, reduced.dimensions=character(0), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(ref$qc.rna$subset.proportion.mito, output$qc.rna$subset.proportion.mito)
})

test_that("runSolo works with regex for mitochondrial genes", {
    copy <- se
    new.ids <- sprintf("%s%s", ifelse(is.mito, "mt-", ""), rownames(copy)) 
    rownames(copy) <- new.ids

    tmp <- tempfile()
    output <- runSolo(copy, qc.mito.regex="^mt-", reduced.dimensions=character(0), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(ref$qc.rna$subset.proportion.mito, output$qc.rna$subset.proportion.mito)

    # Also works if you stick it in a symbol.
    copy <- se
    SummarizedExperiment::rowData(copy)$SYMBOL <- new.ids
    output <- runSolo(copy, qc.mito.regex="^mt-", symbol.field="SYMBOL", reduced.dimensions=character(0), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(ref$qc.rna$subset.proportion.mito, output$qc.rna$subset.proportion.mito)
})

test_that("runSolo works with k-means", {
    tmp <- tempfile()
    kout <- runSolo(se, cluster.method="kmeans", output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(names(kout$markers.rna), levels(kout$sce$kmeans.cluster))
    expect_false(identical(kout$sce$kmeans.cluster, ref$sce$graph.cluster))

    # Other results are unchanged.
    expect_identical(kout$qc.rna, ref$qc.rna)
    expect_identical(SummarizedExperiment::rowData(kout$sce), SummarizedExperiment::rowData(ref$sce))
    expect_identical(SingleCellExperiment::reducedDim(kout$sce, "TSNE"), SingleCellExperiment::reducedDim(ref$sce, "TSNE"))

    # We can actually use multiple methods.
    mult <- runSolo(se, cluster.method=c("kmeans", "graph"), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(mult$sce$kmeans.cluster, kout$sce$kmeans.cluster)
    expect_identical(mult$sce$graph.cluster, ref$sce$graph.cluster)
    expect_identical(mult$markers.rna, kout$markers.rna)
    expect_identical(SingleCellExperiment::reducedDim(mult$sce, "TSNE"), SingleCellExperiment::reducedDim(ref$sce, "TSNE"))

    mult2 <- runSolo(se, cluster.method=c("graph", "kmeans"), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(kout$sce$kmeans.cluster, mult2$sce$kmeans.cluster)
    expect_identical(ref$sce$graph.cluster, mult2$sce$graph.cluster)
    expect_identical(ref$markers.rna, mult2$markers.rna)
    expect_identical(SingleCellExperiment::reducedDim(ref$sce, "TSNE"), SingleCellExperiment::reducedDim(mult2$sce, "TSNE"))
})

test_that("runSolo works with various disabled NN options", {
    tmp <- tempfile()
    graph.only <- runSolo(se, cluster.method="graph", reduced.dimensions=NULL, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(SingleCellExperiment::reducedDimNames(graph.only$sce), "PCA")
    expect_true(is.factor(graph.only$sce$graph.cluster))

    tsne.only <- runSolo(se, cluster.method="kmeans", reduced.dimensions="tsne", output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(SingleCellExperiment::reducedDimNames(tsne.only$sce), c("PCA", "TSNE"))
    expect_null(tsne.only$sce$graph.cluster)

    umap.only <- runSolo(se, cluster.method="kmeans", reduced.dimensions="umap", output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(SingleCellExperiment::reducedDimNames(umap.only$sce), c("PCA", "UMAP"))
    expect_null(umap.only$sce$graph.cluster)

    none <- runSolo(se, cluster.method="kmeans", reduced.dimensions=NULL, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(SingleCellExperiment::reducedDimNames(none$sce), "PCA")
    expect_null(none$sce$graph.cluster)
})

test_that("runSolo works with no clustering at all", {
    tmp <- tempfile()
    noclust <- runSolo(se, cluster.method=NULL, reduced.dimensions=NULL, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_null(noclust$sce$graph.cluster)
    expect_null(noclust$sce$kmeans.cluster)
    expect_identical(ncol(noclust$sce), ncol(ref$sce))
    expect_null(noclust$markers.rna)
})

test_that("runSolo works with a log-fold change threshold", {
    tmp <- tempfile()
    lfc <- runSolo(se, marker.lfc.threshold=1, reduced.dimensions=NULL, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(names(lfc$markers.rna), names(ref$markers.rna))
    delta <- lfc$markers.rna[[1]]$cohens.d.mean - ref$markers.rna[[1]]$cohens.d.mean
    expect_true(all(delta <= 0))
    expect_true(any(delta < 0))
})

test_that("runSolo works with symbols", {
    copy <- se
    SummarizedExperiment::rowData(copy)$symbol <- sprintf("SYMBOL-%s", seq_len(nrow(copy)))

    # Don't suppress plots here, we want to make sure symbols are correctly passed to the heatmap.
    tmp <- tempfile()
    output <- runSolo(copy, symbol.field="symbol", output.dir=tmp, save.results=FALSE, suppress.plots=FALSE)

    expect_identical(SummarizedExperiment::rowData(output$sce)$symbol, SummarizedExperiment::rowData(copy)$symbol)
    expect_identical(output$markers[[1]]$symbol, sub("GENE-", "SYMBOL-", rownames(output$markers[[1]])))
})

test_that("runSolo respects custom metadata", {
    tmp <- tempfile()
    custom <- runSolo(se, output.dir=tmp, metadata=list(custom=list(foo="bar")), suppress.plots=TRUE)

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "sce", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$custom$foo, "bar")
    expect_match(meta$title, "SingleCellExperiment")

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "qc-rna", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$custom$foo, "bar")
    expect_match(meta$title, "Quality control")

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "markers-rna", "1", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$custom$foo, "bar")
    expect_match(meta$title, "Marker gene")
})

test_that("runSolo can do a dry run", {
    tmp <- tempfile()
    output <- runSolo(se, output=tmp, dry.run=TRUE)
    expect_null(output)
    fname <- file.path(tmp, "report.Rmd")
    expect_gt(length(readLines(fname)), 0L)
})
