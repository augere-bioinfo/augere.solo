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

deftmp <- tempfile()
pdf(file=NULL)
default <- runSolo(ase, adt.experiment=TRUE, rna.experiment=NULL, output.dir=deftmp)
dev.off()

test_that("runSolo works with only ADTs", {
    expect_s4_class(default$sce, "SingleCellExperiment")
    expect_identical(nrow(default$sce), ntags)
    expect_lte(ncol(default$sce), ncells)

    expect_gt(mean(default$qc.adt$subset.sum.igg), 0)
    expect_null(SummarizedExperiment::rowData(default$sce)$hvg) # feature selection was skipped.

    expect_identical(SingleCellExperiment::reducedDimNames(default$sce), c("PCA", "TSNE", "UMAP"))
    expect_true(is.factor(default$sce$graph.cluster))
    expect_identical(sort(rownames(default$markers.adt[[1]])), sort(rownames(ase)))

    expect_s4_class(augere.core::readResult(file.path(deftmp, "results", "sce"))$x, "SingleCellExperiment")
    expect_s4_class(augere.core::readResult(file.path(deftmp, "results", "qc-adt"))$x, "DFrame")
    expect_false(file.exists(file.path(deftmp, "results", "qc-rna")))
    expect_s4_class(augere.core::readResult(file.path(deftmp, "results", "markers-adt", "1"))$x, "DFrame")
    expect_false(file.exists(file.path(deftmp, "results", "markers-rna")))

    lines <- readLines(file.path(deftmp, "report.Rmd"))
    expect_false(any(grepl("block", lines, fixed=TRUE))) # shouldn't be any references to blocking.
    expect_false(any(grepl("normalizeRnaCounts.se", lines, fixed=TRUE))) # shouldn't be any references to RNA normalization.
    expect_false(any(grepl("chooseRnaHvgs.se", lines, fixed=TRUE))) # shouldn't be any references to RNA normalization.
})

test_that("runSolo works with ADTs only in an altexp", {
    sce <- as(se, "SingleCellExperiment")
    SingleCellExperiment::altExp(sce, "protein") <- ase

    tmp <- tempfile()
    default2 <- runSolo(sce, adt.experiment="protein", rna.experiment=NULL, output.dir=deftmp, suppress.plots=TRUE, save.results=FALSE)

    expect_identical(default2$sce$graph.cluster, default$sce$graph.cluster)
    expect_identical(SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(default2$sce), "PCA"), SingleCellExperiment::reducedDim(default$sce, "PCA"))
    expect_identical(SingleCellExperiment::reducedDim(default2$sce, "TSNE"), SingleCellExperiment::reducedDim(default$sce, "TSNE"))
    expect_identical(default2$qc.adt, default$qc.adt)
    expect_identical(default2$markers, default$markers)
})

test_that("runSolo works with regex for IgG symbols", {
    copy <- ase
    SummarizedExperiment::rowData(copy)$SYMBOL <- rownames(copy)
    rownames(copy) <- sprintf("whee-%s", seq_len(nrow(copy)))

    tmp <- tempfile()
    symb <- runSolo(copy, symbol.field="SYMBOL", adt.experiment=TRUE, rna.experiment=NULL, reduced.dimensions=character(0), output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)
    expect_identical(default$qc.adt, symb$qc.adt)
})

test_that("runSolo works with RNA plus ADTs", {
    sce <- as(se, "SingleCellExperiment")
    SingleCellExperiment::altExp(sce, "protein") <- ase

    # Don't suppress plots, make sure that the RNA plots are also produced.
    tmp <- tempfile()
    pdf(file=NULL)
    combined <- runSolo(sce, adt.experiment="protein", cluster.method=c("graph", "kmeans"), output.dir=tmp)
    dev.off()

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

    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "sce"))$x, "SingleCellExperiment")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "qc-rna"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "markers-rna", "1"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "qc-adt"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "markers-adt", "1"))$x, "DFrame")
})

test_that("runSolo with ADTs and blocking", {
    ase$batch <- rep(1:4, length.out=ncol(ase))

    tmp <- tempfile()
    blocked <- runSolo(ase, block="batch", adt.experiment=TRUE, rna.experiment=NULL, output.dir=tmp, suppress.plots=TRUE, save.results=FALSE)

    qc.thresh <- S4Vectors::metadata(blocked$sce)$qc$thresholds
    expect_identical(qc.thresh$block.ids, 1:4)
    expect_identical(SingleCellExperiment::reducedDimNames(blocked$sce), c("PCA", "MNN", "TSNE", "UMAP"))

    lines <- readLines(file.path(tmp, "report.Rmd"))
    has.block <- grep("^ *block = ", lines)
    expect_gte(length(has.block), 5) # QC, norm, PCA, MNN, markers.

    # Make sure results are different to the unblocked case.
    expect_identical(blocked$qc.rna[,setdiff(colnames(blocked$qc.rna), "keep")], default$qc.rna[,setdiff(colnames(default$qc.rna), "keep")])
    expect_false(identical(blocked$markers.adt, default$markers.adt))
})
