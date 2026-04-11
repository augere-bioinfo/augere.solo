# library(testthat); library(augere.solo); source("test-annotate.R")

tmat <- matrix(rpois(100000, lambda=5), ncol=100)
test <- SummarizedExperiment::SummarizedExperiment(list(counts = tmat))
rownames(test) <- sprintf("GENE-%s", seq_len(nrow(test)))

rmat <- matrix(rpois(200000, lambda=10), ncol=200)
ref <- SummarizedExperiment::SummarizedExperiment(list(counts = rmat))
rownames(ref) <- sprintf("GENE-%s", seq_len(nrow(ref)))
ref$label <- sample(LETTERS, ncol(ref), replace=TRUE)

deftmp <- tempfile()
pdf(file=NULL)
default <- runAnnotate(
    test,
    configureReferenceAnnotation(ref, "label", ref.assay=1),
    output.dir=deftmp,
    suppress.plots=FALSE
)
dev.off()

test_that("runAnnotate works in the basic case", {
    expect_null(default$combined)

    expect_identical(names(default$predictions), "1")
    expect_identical(nrow(default$predictions[["1"]]), ncol(test)) 
    expect_s4_class(default$predictions[[1]], "DataFrame")
    expect_type(default$predictions[[1]]$labels, "character")
    expect_true(all(default$predictions[[1]]$labels %in% LETTERS))
    expect_gt(length(unlist(S4Vectors::metadata(default$predictions[[1]])$de.genes)), 0) # make sure we actually found some genes.

    expect_s4_class(augere.core::readResult(file.path(deftmp, "results", "pred-1"))$x, "DFrame")
    expect_false(file.exists(file.path(deftmp, "results", "combined")))

    lines <- readLines(file.path(deftmp, "report.Rmd"))
    m <- grep("normalizeRnaCounts.se", lines, fixed=TRUE)
    expect_identical(length(m), 2L) # one for the test, one for the reference.

    expect_false(any(grepl("block", lines)), 0L) # no mentions of blocking should be around.
    expect_false(any(grepl("restrict = ", lines, fixed=TRUE))) # no mentions of intersecting if there's only a single reference.
    expect_false(any(grepl("display.row.names =", lines, fixed=TRUE))) # no need to specify the old rownames for marker diagnostics.
})

test_that("runAnnotate handles pre-normalized test and reference data", {
    tcopy <- scrapper::normalizeRnaCounts.se(test)
    rcopy <- scrapper::normalizeRnaCounts.se(ref)

    tmp <- tempfile()
    out <- runAnnotate(
        tcopy,
        configureReferenceAnnotation(rcopy, "label", ref.assay="logcounts"),
        test.assay="logcounts",
        output.dir=tmp,
        suppress.plots=TRUE,
        save.results=FALSE
    )

    # Ultimately the same results after the same normalization is applied.
    expect_identical(default$predictions[[1]], out$predictions[[1]])

    lines <- readLines(file.path(tmp, "report.Rmd"))
    m <- grep("normalizeRnaCounts.se", lines, fixed=TRUE)
    expect_identical(length(m), 0L) # there should be no mention of normalization.
})

test_that("runAnnotate works with alternative IDs", {
    tcopy <- test
    SummarizedExperiment::rowData(tcopy)$FOO <- rownames(tcopy)
    rownames(tcopy) <- sprintf("ID-%s", seq_len(nrow(tcopy)))

    rcopy <- ref
    SummarizedExperiment::rowData(rcopy)$BAR <- rownames(rcopy)
    rownames(rcopy) <- sprintf("WHEE-%s", seq_len(nrow(rcopy)))

    tmp <- tempfile()
    out <- runAnnotate(
        tcopy,
        configureReferenceAnnotation(rcopy, "label", ref.assay=1, ref.id.field="BAR"),
        test.id.field="FOO",
        output.dir=tmp,
        suppress.plots=FALSE, # don't suppress to make sure the heatmaps have the correct IDs.
        save.results=FALSE
    )

    # Choice of IDs should still end up giving the same results.
    expect_identical(default$predictions[[1]], out$predictions[[1]])

    lines <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grepl("FOO", lines)))
    expect_true(any(grepl("BAR", lines)))
    expect_true(any(grepl("display.row.names = old.test.rownames", lines, fixed=TRUE)))
})

test_that("runAnnotate handles blocking in the test", {
    copy <- test
    copy$BOO <- rep(1:3, length.out=ncol(test))

    tmp <- tempfile()
    out <- runAnnotate(
        copy,
        configureReferenceAnnotation(ref, "label", ref.assay=1),
        test.block.field = "BOO",
        output.dir=tmp,
        suppress.plots=FALSE, # don't suppress to make sure block is correctly passed in plotMarkerHeatmap.
        save.results=FALSE
    )

    # Blocking in the test has no effect on the results.
    expect_identical(default$predictions[[1]], out$predictions[[1]])

    lines <- readLines(file.path(tmp, "report.Rmd"))
    m <- grep("block = ", lines)
    expect_identical(length(m), 2L) # once for normalization, another for the marker heatmap.
})

test_that("runAnnotate works for symbol diagnostics", {
    tcopy <- test
    SummarizedExperiment::rowData(tcopy)$FOO <- sprintf("ID-%s", seq_len(nrow(tcopy)))

    tmp <- tempfile()
    out <- runAnnotate(
        tcopy,
        configureReferenceAnnotation(ref, "label", ref.assay=1),
        test.symbol.field = "FOO",
        output.dir=tmp,
        suppress.plots=FALSE, # don't suppress to make sure symbols are correctly used in plotMarkerHeatmap.
        save.results=FALSE
    )

    # Symbol only affects the marker heatmap diagnostics, so the results should be unchanged.
    expect_identical(out, default)

    lines <- readLines(file.path(tmp, "report.Rmd"))
    expect_true(any(grep("display.row.names = SummarizedExperiment::rowData(test)[,\"FOO\"]", lines, fixed=TRUE)))
})

test_that("runAnnotate constructs extra visualizations", {
    copy <- as(test, "SingleCellExperiment")
    copy$cluster <- rep(letters[1:5], length.out=ncol(test))
    SingleCellExperiment::reducedDim(copy, "TSNE") <- matrix(runif(ncol(test) * 2), ncol=2)

    tmp <- tempfile()
    out <- runAnnotate(
        copy,
        configureReferenceAnnotation(ref, "label", ref.assay=1),
        output.dir=tmp,
        cluster.field="cluster",
        reduced.dimensions="TSNE",
        suppress.plots=FALSE, # don't suppress to make plots are correctly made.
        save.results=FALSE
    )

    # Adding these extra fields has no effect on the results.
    expect_identical(default$predictions[[1]], out$predictions[[1]])

    lines <- readLines(file.path(tmp, "report.Rmd"))
    m <- grep("plotReducedDim", lines)
    expect_identical(length(m), 1L)
})

test_that("runAnnotate handles blocking in the reference", {
    copy <- ref
    copy$BAR <- rep(1:3, length.out=ncol(ref))

    tmp <- tempfile()
    out <- runAnnotate(
        test,
        configureReferenceAnnotation(copy, "label", ref.assay=1, ref.block.field="BAR"), 
        output.dir=tmp,
        suppress.plots=TRUE,
        save.results=FALSE
    )

    # Marker detection might be different when blocking on the reference, so we just do more cursory checks here.
    expect_s4_class(out$predictions[[1]], "DataFrame")
    expect_identical(nrow(out$predictions[["1"]]), ncol(test)) 

    lines <- readLines(file.path(tmp, "report.Rmd"))
    m <- grep("block = ", lines)
    expect_identical(length(m), 2L) # once for normalization, another for trainSingleR itself.
})

test_that("runAnnotate works with multiple references", {
    # Normalize first otherwise the results will be slightly different if we subset and then normalize.
    # We also cut down the number of samples for speed as combined annotations are currently very slow. 
    tcopy <- scrapper::normalizeRnaCounts.se(test[,1:20])
    rcopy <- scrapper::normalizeRnaCounts.se(ref)

    rmat2 <- matrix(rpois(150 * 1000, lambda=10), ncol=150)
    ref2 <- SummarizedExperiment::SummarizedExperiment(list(counts = rmat2))
    rownames(ref2) <- sprintf("GENE-%s", seq_len(nrow(ref2)))
    ref2$assigned <- sample(letters, ncol(ref2), replace=TRUE)
    ref2 <- scrapper::normalizeRnaCounts.se(ref2)

    tmp <- tempfile()
    out <- runAnnotate(
        tcopy[1:800,],
        list(
            FOO=configureReferenceAnnotation(rcopy[101:900,], "label", ref.assay="logcounts"),
            BAR=configureReferenceAnnotation(ref2[201:1000,], "assigned", ref.assay="logcounts")
        ),
        test.assay="logcounts",
        output.dir=tmp,
        suppress.plots=TRUE
    )

    expect_true(all(out$predictions$FOO$labels %in% LETTERS))
    expect_true(all(out$predictions$BAR$labels %in% letters))
    expect_true(all(out$combined$labels %in% c(LETTERS, letters)))
    expect_true(all(out$combined$reference %in% 1:2))

    all.genes.1 <- unlist(S4Vectors::metadata(out$predictions[[1]])$de.genes)
    expect_gt(length(all.genes.1), 0)
    expect_true(all(all.genes.1 %in% rownames(tcopy)[201:800]))
    all.genes.2 <- unlist(S4Vectors::metadata(out$predictions[[2]])$de.genes)
    expect_gt(length(all.genes.2), 0)
    expect_true(all(all.genes.2 %in% rownames(tcopy)[201:800]))

    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "pred-FOO"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "pred-BAR"))$x, "DFrame")
    expect_s4_class(augere.core::readResult(file.path(tmp, "results", "combined"))$x, "DFrame")

    # Comparing it to manual slicing of the feature space.
    manual <- runAnnotate(
        tcopy[201:800,],
        list(
            FOO=configureReferenceAnnotation(rcopy[201:800,], "label", ref.assay="logcounts"),
            BAR=configureReferenceAnnotation(ref2[201:800,], "assigned", ref.assay="logcounts")
        ),
        test.assay="logcounts",
        output.dir=tmp,
        suppress.plots=TRUE,
        save.results=FALSE
    )

    expect_identical(out, manual)
})

test_that("runAnnotate works with celldex references", {
    ref.inbuilt <- celldex::BlueprintEncodeData()

    tmat <- matrix(rpois(100000, lambda=5), ncol=100)
    test <- SummarizedExperiment::SummarizedExperiment(list(counts = tmat))
    rownames(test) <- sample(rownames(ref.inbuilt), nrow(test))

    tmp <- tempfile()
    out <- runAnnotate(
        test,
        configureReferenceAnnotation("BlueprintEncodeData", "label.main"),
        output.dir=tmp,
        suppress.plots=TRUE,
        save.results=FALSE
    )
    expect_gt(length(unlist(S4Vectors::metadata(out$predictions[[1]])$de.genes)), 0) # make sure we actually found some genes.

    tmp2 <- tempfile()
    manual <- runAnnotate(
        test,
        configureReferenceAnnotation(ref.inbuilt, "label.main", ref.marker.method="classic"),
        output.dir=tmp2,
        suppress.plots=TRUE,
        save.results=FALSE
    )
    expect_identical(out, manual)

    # Also works with Ensembl.
    suppressWarnings(ref.inbuilt <- celldex::BlueprintEncodeData(ensembl=TRUE))
    rownames(test) <- sample(rownames(ref.inbuilt), nrow(test))

    tmp <- tempfile()
    suppressWarnings(out <- runAnnotate(
        test,
        configureReferenceAnnotation("BlueprintEncodeData", "label.main"),
        output.dir=tmp,
        test.is.ensembl=TRUE,
        suppress.plots=TRUE,
        save.results=FALSE
    ))

    expect_s4_class(out$predictions[[1]], "DataFrame")
    expect_identical(nrow(out$predictions[[1]]), ncol(test))
    expect_gt(length(unlist(S4Vectors::metadata(out$predictions[[1]])$de.genes)), 0) # make sure we actually found some genes.
})

test_that("runAnnotate respects custom metadata", {
    tmp <- tempfile()
    out <- runAnnotate(
        test[,1:20], # cutting down the number of test cells to save time.
        list(
            configureReferenceAnnotation(ref, "label", ref.assay=1),
            configureReferenceAnnotation(ref, "label", ref.assay=1)
        ),
        metadata=list(custom=list(foo="bar")),
        output.dir=tmp,
        suppress.plots=TRUE
    )

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "pred-1", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$custom$foo, "bar")
    expect_match(meta$title, "reference `1`")

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "pred-2", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$custom$foo, "bar")
    expect_match(meta$title, "reference `2`")

    meta <- jsonlite::fromJSON(file.path(tmp, "results", "combined", "_metadata.json"), simplifyVector=FALSE)
    expect_identical(meta$custom$foo, "bar")
    expect_match(meta$title, "Combined")
})

test_that("runAnnotate can do a dry run", {
    tmp <- tempfile()
    output <- runAnnotate(
        test,
        configureReferenceAnnotation(ref, "label", ref.assay=1),
        output.dir=tmp,
        dry.run=TRUE
    )

    expect_null(output)
    fname <- file.path(tmp, "report.Rmd")
    expect_gt(length(readLines(fname)), 0L)
})
