#' Commands for limma's mean-rank gene set test
#'
#' Generate R commands to run \code{\link[limma]{geneSetTest}} under various parametrizations.
#'
#' @inheritParams generateHypergeometricTestCommands
#' @param stat.name String containing the variable name of the numeric vector containing the test statistic for ranking genes.
#' The referenced vector should be of length equal to the number of genes.
#' @param sign.name String containing the variable name for the numeric vector of signed effect sizes (typically log-fold changes) for each gene.
#' The referenced vector should be of length equal to the number of genes.
#' It can be \code{NULL} if the statistics referenced by \code{stat.name} is already signed or \code{alternative="mixed"}.
#' @param use.sqrt Boolean indicating whether to take the square root of the unsigned statistics when converting them to signed values.
#' @param alternative String specifying the alternative hypothesis. 
#' This should be one of:
#' \itemize{
#' \item \code{"mixed"}: test each set for enrichment of genes with the largest absolute values of the statistics. 
#' \item \code{"up"}: test each set for enrichment of genes with the largest statistics, typically up-regulated.
#' \item \code{"down"}: test each set for enrichment of genes with the smallest statistics, typically down-regulated.
#' \item \code{"either"}: test each set for enrichment of genes in either direction.
#' }
#' @param seed Seed for the random number generator.
#' If \code{NULL}, no seed setting is performed.
#' @param args Named list of arguments to pass to \code{\link[limma]{geneSetTest}}.
#'
#' @section Direction of the alternative hypothesis:
#' In general, it is expected that the statistics in \code{stat.name} are signed values from a differential expression analysis, e.g., t-statistics, z-scores.
#' Larger positive and negative values should correspond to stronger up- and down-regulation, respectively.
#' This ensures that \code{alternative="up"} or \code{alternative="down"} will have the expected behavior.
#'
#' If the input statistics are unsigned, they are converted to their signed equivalents by multiplying it with the sign of the variable referenced by \code{sign}.
#' For F-statistics or chi-squared statistics, we can set \code{use.sqrt=TRUE} to convert them to t-statistics and z-scores, respectively.
#' This mimics the behavior of a test that produces signed test statistics.
#'
#' For \code{alternative="either"}, each gene set is tested separately for enrichment of up- and down-regulated genes.
#' The reported \code{Direction} is defined based on the one with a lower p-value.
#' In contrast, \code{alternative="mixed"} does not consider the sign at all, and only cares about enrichment of significant genes of any direction in each set.
#'
#' @author Aaron Lun
#'
#' @return Character vector containing the commands required to run \code{geneSetTest}.
#' Upon evaluation, this produces a \link[S4Vectors]{DFrame} with one row per gene set and the following columns:
#' \itemize{
#' \item \code{NumGenes}, the number of genes in each set.
#' \item \code{PValue}, the p-value for enrichment of highly-ranked genes in this set.
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value.
#' \item \code{Direction}, whether the set is more enriched for up- or down-regulated genes.
#' Only reported for \code{alternative="either"}.
#' }
#'
#' @examples
#' cat(generateGeneSetTestCommands("sets", "stat", "sign"), sep="\n")
#' cat(generateGeneSetTestCommands("sets", "stat", "sign", alternative="either"), sep="\n")
#'
#' @export
#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateGeneSetTestCommands <- function(sets.name, stat.name, sign.name, alternative=c("mixed", "up", "down", "either"), seed=NULL, use.sqrt=FALSE, args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    stat <- <%= STAT %>
:BEGIN sign-setup
    sign <- <%= SIGN %>
:END

    sizes <- unname(lengths(sets))
    keep <- sizes != 0L
    nonempty.sets <- unname(sets)[keep]

:BEGIN re-sign
    # Re-introducing the sign.
    stat <- <%= SANITIZER %> * base::sign(sign)

:END
:BEGIN seed
    set.seed(<%= SEED %>)
:END
    all.p <- rep(NA_real_, length(sets))
    all.p[keep] <- vapply(
        nonempty.sets,
        limma::geneSetTest,
        statistics=stat,
        alternative=<%= ALTERNATIVE %>,
        type=<%= TYPE %>,
:BEGIN more-args
:END
        FUN.VALUE=0
    )

    out <- S4Vectors::DataFrame(row.names=names(sets), NumGenes=sizes, PValue=all.p)
    out$FDR <- p.adjust(out$PValue, method='BH')
:BEGIN more-alt
    direction <- vapply(
        sets,
        function(x) mean(stat[x]) > 0,
        FUN.VALUE=TRUE,
        USE.NAMES=FALSE
    )
    out$Direction <- ifelse(direction, \"up\", \"down\")
:END
    out
})"

    parsed <- parseRmdTemplate(template)

    alternative <- match.arg(alternative)
    replacements <- list(ALTERNATIVE = deparseToString(alternative), SETS=sets.name, STAT=stat.name, SIGN=sign.name)
    if (!is.null(seed)) {
        parsed$seed <- NULL
    } else {
        replacements$SEED <- deparseToString(seed)
    }

    if (alternative == "mixed") {
        # geneSetTest just takes the absolute value internally so the sign doesn't matter.
        parsed[["re-sign"]] <- NULL
        parsed[["sign-setup"]] <- NULL
        replacements$TYPE <- "\"f\""

    } else {
        if (is.null(sign.name)) {
            parsed[["re-sign"]] <- NULL
            parsed[["sign-setup"]] <- NULL
        } else {
            if (use.sqrt) {
                replacements$SANITIZER <- "sqrt(abs(stat))"
            } else {
                replacements$SANITIZER <- "abs(stat)"
            }
        }
        replacements$TYPE <- "\"t\""
    }

    if (alternative != "either") {
        parsed[["more-alt"]] <- NULL
    }

    if (length(args)) {
        parsed[["more-args"]] <- .format_more_args(args, comma.last=TRUE)
    }

    formatted <- replacePlaceholders(parsed, replacements)
    unname(unlist(formatted, use.names=FALSE))
}
