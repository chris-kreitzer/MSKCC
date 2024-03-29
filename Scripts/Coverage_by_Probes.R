#' Function to calculate coverage from BAM file
#'
#' Takes a BAM file and an interval file as input and returns coverage for each
#' interval. Coverage should be then GC-normalized using the
#' \code{\link{correctCoverageBias}} function before determining purity and
#' ploidy with \code{\link{runAbsoluteCN}}. Uses the \code{scanBam} function
#' and applies low quality, duplicate reads as well as secondary alignment
#' filters.
#'
#'
#' @param bam.file Filename of a BAM file.
#' @param interval.file File specifying the intervals. Interval is expected in
#' first column in format CHR:START-END. 
#' @param output.file Optionally, write minimal coverage file. Can be read with
#' the \code{\link{readCoverageFile}} function.
#' @param index.file The bai index. This is expected without the .bai file
#' suffix, see \code{?scanBam}.
#' @param keep.duplicates Keep or remove duplicated reads.
#' @param chunks Split \code{interval.file} into specified number of chunks
#' to reduce memory usage.
#' @param ... Additional parameters passed to \code{ScanBamParam}.
#' @return Returns total and average coverage by intervals.
#' @author Markus Riester
#' @seealso \code{\link{preprocessIntervals}
#' \link{correctCoverageBias} \link{runAbsoluteCN}}
#' @examples
#'
#' bam.file <- system.file("extdata", "ex1.bam", package = "PureCN",
#'     mustWork = TRUE)
#' interval.file <- system.file("extdata", "ex1_intervals.txt",
#'     package = "PureCN", mustWork = TRUE)
#'
#' # Calculate raw coverage from BAM file. These need to be corrected for
#' # GC-bias using the correctCoverageBias function before determining purity
#' # and ploidy.
#' coverage <- calculateBamCoverageByInterval(bam.file = bam.file,
#'     interval.file = interval.file)
#'
#' @export calculateBamCoverageByInterval
#' @importFrom Rsamtools headerTabix ScanBamParam scanBamFlag
#'             scanBam scanFa scanFaIndex TabixFile
calculateBamCoverageByInterval <- function(bam.file, interval.file,
    output.file = NULL, index.file = bam.file, keep.duplicates = FALSE,
    chunks = 20, ...) {
    flog.debug("Used Ncells/Vcells: %s (Mb)",
        paste(suppressMessages(gc()[,2]), collapse = "/"))

    intervalGrAll <- readIntervalFile(interval.file, strict = FALSE,
        verbose = FALSE)
    # skip splitting with small panels
    if (is.null(chunks) || is.na(chunks) || chunks < 2 ||
        length(intervalGrAll) < chunks * 50) {
        chunkIds <- rep(1, length(intervalGrAll))
    } else {
        chunkIds <- cut(seq_along(intervalGrAll), chunks, labels = FALSE)
    }
    .calculateBamCoverageByInterval <- function(intervalGr, ...) {
        param <- ScanBamParam(what = c("pos", "qwidth", "flag"),
                    which = intervalGr,
                    flag = scanBamFlag(isUnmappedQuery = FALSE,
                                 isNotPassingQualityControls = FALSE,
                                 isSecondaryAlignment = FALSE,
                                 isDuplicate = NA
                         ),
                    ...
                 )

        xAll <- scanBam(bam.file, index = index.file, param = param)
        xDupFiltered <- .filterDuplicates(xAll)

        x <- xDupFiltered
        if (keep.duplicates) x <- xAll

        intervalGr$coverage <- vapply(seq_along(x), function(i)
            sum(coverage(IRanges(x[[i]][["pos"]], width = x[[i]][["qwidth"]]),
                shift = -start(intervalGr)[i], width = width(intervalGr)[i])), integer(1))

        intervalGr$average.coverage <- intervalGr$coverage / width(intervalGr)

        intervalGr$counts <- as.numeric(vapply(x, function(y) length(y$pos), integer(1)))
        intervalGr$duplication.rate <- 1 -
            vapply(xDupFiltered, function(y) length(y$pos), integer(1)) /
            vapply(xAll, function(y) length(y$pos), integer(1))

        flog.debug("Used Ncells/Vcells: %s (Mb)",
            paste(suppressMessages(gc()[,2]), collapse = "/"))
        
        x <- xAll <- xDupFiltered <- NULL    
        intervalGr
    }
    intervalGrAll <- split(intervalGrAll, chunkIds)

    intervalGr <- unlist(GRangesList(
        lapply(seq_along(intervalGrAll), function(i) {
            if (length(intervalGrAll) > 1) {
                flog.info("Processing %s (%i/%i)...", as.character(intervalGrAll[[i]][1]),
                    i, length(intervalGrAll))
            }
            .calculateBamCoverageByInterval(intervalGrAll[[i]], ...)
        }
    )))

    if (!is.null(output.file)) {
        .writeCoverage(intervalGr, output.file)
    }
    invisible(intervalGr)
}

.writeCoverage <- function(intervalGr, output.file) {
    tmp <- data.frame(
        Target = as.character(intervalGr),
        total_coverage = intervalGr$coverage,
        counts = intervalGr$counts,
        on_target = intervalGr$on.target,
        duplication_rate = intervalGr$duplication.rate
    )
    fwrite(tmp, file = output.file, row.names = FALSE, quote = FALSE,
        sep = " ", logical01 = TRUE, na = "NA")
}

.filterDuplicates <- function(x) {
    lapply(x, function(y) {
        idx <- y$flag < 1024
        lapply(y, function(z) z[idx])
    })
}