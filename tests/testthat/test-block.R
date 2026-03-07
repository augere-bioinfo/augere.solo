# library(testthat); library(augere.solo); source("test-block.R")

set.seed(88)
ngenes <- 10000L
ncells <- 200L
mu <- 2^rnorm(10000)
y <- matrix(rpois(ngenes * ncells, lambda=mu), ncol=ncells)
se <- SummarizedExperiment::SummarizedExperiment(list(counts=y))

is.mito <- rbinom(ngenes, 1, 0.01) == 1
SummarizedExperiment::rowRanges(se) <- GenomicRanges::GRanges(c("chrA", "chrM")[is.mito + 1], IRanges::IRanges(1, 1))
rownames(se) <- sprintf("GENE-%s", seq_len(ngenes))

se$block <- sample(LETTERS[1:3], ncells, replace=TRUE)

tmp <- tempfile()
ref <- runSolo(se, block.field="block", output.dir=tmp, save.results=FALSE)

test_that("runSolo works with blocking", {
    expect_s4_class(ref$sce, "SingleCellExperiment")
    expect_identical(nrow(ref$sce), ngenes)
    expect_lt(ncol(ref$sce), ncells)

    qc.thresh <- S4Vectors::metadata(ref$sce)$qc$thresholds
    expect_gt(length(qc.thresh$sum), 0L)
    expect_identical(qc.thresh$block.ids, sort(unique(se$block)))

    expect_identical(SummarizedExperiment::assayNames(ref$sce), c("counts", "logcounts"))
    expect_identical(SingleCellExperiment::reducedDimNames(ref$sce), c("PCA", "MNN", "TSNE", "UMAP"))
    expect_true(is.factor(ref$sce$graph.cluster))

    is.hvg <- SummarizedExperiment::rowData(ref$sce)$hvg
    expect_true(any(is.hvg))
    expect_false(all(is.hvg))
    expect_identical(rownames(S4Vectors::metadata(ref$sce)$PCA$rotation), rownames(ref$sce)[is.hvg])

    expect_s4_class(SingleCellExperiment::counts(ref$sce), "DelayedMatrix")
    expect_s4_class(SingleCellExperiment::logcounts(ref$sce), "LogNormalizedMatrix")

    expect_identical(names(ref$markers.rna), levels(ref$sce$graph.cluster))
    expect_s4_class(ref$markers.rna[[1]], "DFrame")
    expect_identical(sort(rownames(ref$markers.rna[[1]])), sort(rownames(se)))

    # Check that we actually specify block= in the various calls.
    lines <- readLines(file.path(tmp, "report.Rmd"))
    has.block <- grep("^ *block = ", lines)
    expect_gte(length(has.block), 6) # QC, norm, HVG, PCA, MNN, markers.

    # Check that we actually use MNN in all subsequent reddim specifications. 
    reddim <- grep("^ *reddim.type = ", lines)
    expect_gte(length(reddim), 2L)
    expect_match(lines[reddim[1]], "PCA")
    expect_true(all(grep("MNN", lines[reddim[-1]])))

    # Check that results are actually different compared to the unblocked case.
    utmp <- tempfile()
    unblocked <- runSolo(se, output.dir=utmp, save.results=FALSE)
    expect_identical(ref$qc.rna[,setdiff(colnames(ref$qc.rna), "keep")], unblocked$qc.rna[,setdiff(colnames(unblocked$qc.rna), "keep")])
    expect_false(identical(ref$markers.rna, unblocked$markers.rna))

    ulines <- readLines(file.path(utmp, "report.Rmd"))
    uhas.block <- grep("block", ulines)
    expect_identical(length(uhas.block), 0L) # there shouldn't be any references to blocking here.
})

test_that("runSolo's kmeans behaves with blocking", {
    tmp <- tempfile()
    kout <- runSolo(se, block.field="block", cluster.method="kmeans", output.dir=tmp, save.results=FALSE)
    expect_identical(kout$qc.rna, ref$qc.rna)
    expect_identical(SingleCellExperiment::reducedDimNames(kout$sce), c("PCA", "MNN", "TSNE", "UMAP"))
    expect_null(kout$sce$graph.cluster)
    expect_true(is.factor(kout$sce$kmeans.cluster))

    # Check that we actually use MNN in all reddim specifications. 
    lines <- readLines(file.path(tmp, "report.Rmd"))
    reddim <- grep("^ *reddim.type = ", lines)
    expect_gte(length(reddim), 3L)
    expect_match(lines[reddim[1]], "PCA")
    expect_true(all(grep("MNN", lines[reddim[-1]])))
})

test_that("runSolo works if no clustering is performed", {
    # Specifically test whether the t-SNE and UMAP plots are configured correctly,
    # namely that they color by block instead of cluster.
    noclust <- runSolo(se, block.field="block", cluster.method=NULL, output.dir=tmp, save.results=FALSE)
    expect_null(noclust$sce$graph.cluster)
    expect_null(noclust$sce$kmeans.cluster)
    expect_null(noclust$sce$markers.rna)

    expect_identical(SingleCellExperiment::reducedDimNames(noclust$sce), c("PCA", "MNN", "TSNE", "UMAP"))
    expect_identical(SingleCellExperiment::reducedDim(noclust$sce, "MNN"), SingleCellExperiment::reducedDim(noclust$sce, "MNN"))
})

