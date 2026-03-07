# library(testthat); library(augere.solo); source("test-adt.R")

set.seed(123123)
ngenes <- 10000L
ncells <- 200L
mu <- 2^rnorm(10000)
y <- matrix(rpois(ngenes * ncells, lambda=mu), ncol=ncells)
se <- SummarizedExperiment::SummarizedExperiment(list(counts=y))
is.mito <- rbinom(ngenes, 1, 0.01) == 1
SummarizedExperiment::rowRanges(se) <- GenomicRanges::GRanges(c("chrA", "chrM")[is.mito + 1], IRanges::IRanges(1, 1))
rownames(se) <- sprintf("GENE-%s", seq_len(ngenes))

ntags <- 100L
ay <- matrix(rpois(ntags * ncells, lambda=100), ncol=ncells)
ase <- SummarizedExperiment::SummarizedExperiment(list(counts=ay))
rownames(ase) <- sprintf("TAG-%s", seq_len(ntags))
rownames(ase)[1:5] <- sprintf("IgG-%s", 1:5)

test_that("runSolo works with only ADTs", {
    tmp <- tempfile()
    adt <- runSolo(ase, adt.experiment=TRUE, rna.experiment=NULL, output.dir=tmp, save.results=FALSE)

    expect_s4_class(adt$sce, "SingleCellExperiment")
    expect_identical(nrow(adt$sce), ntags)
    expect_lte(ncol(adt$sce), ncells)

    expect_gt(mean(adt$qc.adt$subset.sum.igg), 0)
    expect_null(SummarizedExperiment::rowData(adt$sce)$hvg) # feature selection was skipped.

    expect_identical(SingleCellExperiment::reducedDimNames(adt$sce), c("PCA", "TSNE", "UMAP"))
    expect_true(is.factor(adt$sce$graph.cluster))
    expect_identical(sort(rownames(adt$markers.adt[[1]])), sort(rownames(ase)))

    # Still works if the ADTs are tucked into an altexp.
    sce <- as(se, "SingleCellExperiment")
    SingleCellExperiment::altExp(sce, "protein") <- ase
    adt2 <- runSolo(sce, adt.experiment="protein", rna.experiment=NULL, output.dir=tmp, save.results=FALSE)

    expect_identical(adt2$sce$graph.cluster, adt$sce$graph.cluster)
    expect_identical(SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(adt2$sce), "PCA"), SingleCellExperiment::reducedDim(adt$sce, "PCA"))
    expect_identical(SingleCellExperiment::reducedDim(adt2$sce, "TSNE"), SingleCellExperiment::reducedDim(adt$sce, "TSNE"))
    expect_identical(adt2$qc.adt, adt$qc.adt)
    expect_identical(adt2$markers, adt$markers)
})

test_that("runSolo works with regex for IgG symbols", {
    copy <- ase
    tmp <- tempfile()
    adt <- runSolo(ase, adt.experiment=TRUE, rna.experiment=NULL, reduced.dimensions=character(0), output.dir=tmp, save.results=FALSE)

    SummarizedExperiment::rowData(copy)$SYMBOL <- rownames(copy)
    rownames(copy) <- sprintf("whee-%s", seq_len(nrow(copy)))
    symb <- runSolo(copy, symbol.field="SYMBOL", adt.experiment=TRUE, rna.experiment=NULL, reduced.dimensions=character(0), output.dir=tmp, save.results=FALSE)

    expect_identical(adt$qc.adt$subset.sum.igg, symb$qc.adt$subset.sum.igg)
})

test_that("runSolo works with RNA plus ADTs", {
    sce <- as(se, "SingleCellExperiment")
    SingleCellExperiment::altExp(sce, "protein") <- ase

    tmp <- tempfile()
    combined <- runSolo(sce, adt.experiment="protein", cluster.method=c("graph", "kmeans"), output.dir=tmp, save.results=FALSE)

    expect_identical(
        ncol(SingleCellExperiment::reducedDim(combined$sce, "combined")),
        ncol(SingleCellExperiment::reducedDim(combined$sce, "PCA")) + ncol(SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(combined$sce), "PCA"))
    )
    expect_identical(names(combined$markers.rna), levels(combined$sce$graph.cluster))
    expect_identical(names(combined$markers.rna), names(combined$markers.adt))

    lines <- readLines(file.path(tmp, "report.Rmd"))
    rd.lines <- lines[grep("reddim.type", lines)]
    expect_identical(length(rd.lines), 2L) # one for k-means, one for the neighbor-related steps.
    expect_true(all(grepl("combined", rd.lines)))
})

test_that("runSolo with ADTs and blocking", {
    ase$batch <- rep(1:4, length.out=ncol(ase))

    tmp <- tempfile()
    adt <- runSolo(ase, block="batch", adt.experiment=TRUE, rna.experiment=NULL, output.dir=tmp, save.results=FALSE)

    qc.thresh <- S4Vectors::metadata(adt$sce)$qc$thresholds
    expect_identical(qc.thresh$block.ids, 1:4)
    expect_identical(SingleCellExperiment::reducedDimNames(adt$sce), c("PCA", "MNN", "TSNE", "UMAP"))

    lines <- readLines(file.path(tmp, "report.Rmd"))
    has.block <- grep("^ *block = ", lines)
    expect_gte(length(has.block), 5) # QC, norm, PCA, MNN, markers.

    # Make sure results are different to the unblocked case.
    utmp <- tempfile()
    unblocked <- runSolo(ase, adt.experiment=TRUE, rna.experiment=NULL, output.dir=utmp, save.results=FALSE)
    expect_identical(adt$qc.rna[,setdiff(colnames(adt$qc.rna), "keep")], unblocked$qc.rna[,setdiff(colnames(unblocked$qc.rna), "keep")])
    expect_false(identical(adt$markers.adt, unblocked$markers.adt))

    ulines <- readLines(file.path(utmp, "report.Rmd"))
    uhas.block <- grep("block", ulines)
    expect_identical(length(uhas.block), 0L) # shouldn't be any references to blocking.
})

test_that("runSolo saves ADT data correctly", {
    tmp <- tempfile()
    saved <- runSolo(ase, adt.experiment=TRUE, rna.experiment=NULL, output.dir=tmp)
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "sce"))$x, "SingleCellExperiment")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "qc-adt"))$x, "DFrame")
    expect_false(file.exists(file.path(tmp, "results", "qc-rna")))
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "markers-adt", "1"))$x, "DFrame")
    expect_false(file.exists(file.path(tmp, "results", "markers-rna")))

    # Now saving with both RNA and ADTs together.
    sce <- as(se, "SingleCellExperiment")
    SingleCellExperiment::altExp(sce, "protein") <- ase
    tmp <- tempfile()
    combined <- runSolo(sce, adt.experiment="protein", output.dir=tmp)
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "sce"))$x, "SingleCellExperiment")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "qc-rna"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "markers-rna", "1"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "qc-adt"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "markers-adt", "1"))$x, "DFrame")
})
