.define_mito_rows_code <- function(se.name) {
    sprintf("is.mito <- local({
    rr <- SummarizedExperiment::rowRanges(%s)
    mito.names <- c('MT', 'M', 'chrM', 'chrMT')
    mito <- !is.na(S4Vectors::match(Seqinfo::seqnames(rr), mito.names))
    if (is(rr, 'GenomicRanges')) {
        S4Vectors::which(mito)
    } else {
        any(mito)
    }
})", se.name)
}

