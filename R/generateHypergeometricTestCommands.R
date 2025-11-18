#' Commands for a hypergeometric test
#'
#' Generate R commands to run a hypergeometric test under various parameterizations.
#' 
#' @param sets.name String containing the variable name for the list of gene sets.
#' Each set should be an integer vector of indices into the variable named by \code{is.sig.name}.
#' @param is.sig.name String containing the variable name for the logical vector indicating which genes are significant.
#' The referenced vector should be of length equal to the number of genes.
#' @param sign.name String containing the variable name for the numeric vector of signed effect sizes (typically log-fold changes) for each gene.
#' The referenced vector should be of length equal to the number of genes.
#' It can be \code{NULL} for \code{alternative="mixed"}.
#' @param alternative String specifying the alternative hypothesis.
#' This should be one of:
#' \itemize{
#' \item \code{"mixed"} tests for enrichment of any significant genes in each set, regardless of their direction.
#' \item \code{"up"} tests for enrichment of up-regulated genes in each set.
#' \item \code{"down"} tests for enrichment of down-regulated genes in each set.
#' \item \code{"either"} tests for enrichment of either up- or down-regulated genes in each set.
#' }
#'
#' @return Character vector containing the commands required to perform a hypergeometric test.
#' Upon evaluation, this produces a \link[S4Vectors]{DFrame} with one row per gene set and the following columns:
#' \itemize{
#' \item \code{NumGenes}, the number of genes in each set. 
#' \item \code{NumSig}, the number of significant genes in each set, possibly after filtering for the desired sign.
#' \item \code{PValue}, the p-value for enrichment from the hypergeometric test.
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value.
#' \item \code{Direction}, whether the set is more enriched for up- or down-regulated genes.
#' Only reported for \code{alternative="either"}.
#' }
#'
#' @section Direction of the alternative hypothesis:
#' For \code{alternative="either"}, each gene set is tested separately for enrichment of up- and down-regulated genes.
#' The reported \code{Direction} is defined based on the one with a lower p-value.
#' In contrast, \code{alternative="mixed"} does not consider the sign at all, and only cares about enrichment of significant genes of any direction in each set.
#'
#' @author Aaron Lun
#' 
#' @examples
#' cat(generateHypergeometricTestCommands("sets", "is.sig", "sign"), sep="\n")
#' cat(generateHypergeometricTestCommands("sets", "is.sig", "sign", alternative="either"), sep="\n")
#'
#' @export
#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateHypergeometricTestCommands <- function(sets.name, is.sig.name, sign.name, alternative=c("mixed", "up", "down", "either")) {
    template <- "local({
    # Simple implementation of the hypergeometric test.
    run_hypergeometric_test <- function(sets, is.sig) {
        num.universe <- length(is.sig)
        num.sig <- sum(is.sig)
        num.set <- unname(lengths(sets))
        num.set.sig <- vapply(sets, function(x) sum(is.sig[x]), FUN.VALUE=0L, USE.NAMES=FALSE)

        # Make sure to include the probability mass of the observed count in the p-value.
        args <- list(m=num.sig, n=num.universe - num.sig, k=num.set)
        pval <- do.call(phyper, c(args, list(q=num.set.sig, lower.tail=FALSE))) + do.call(dhyper, c(args, list(x=num.set.sig)))
        pval[num.set == 0L] <- NA # make sure empty sets do not affect the FDR correction.

        out <- S4Vectors::DataFrame(row.names=names(sets), NumGenes=num.set, NumSig=num.set.sig, PValue=pval)
        out$FDR <- p.adjust(out$PValue, method='BH')
        out
    }

    sets <- <%= SETS %>
    is.sig <- <%= IS_SIG %>
:BEGIN sign-setup
    sign <- <%= SIGN %>
:END

:BEGIN up
    run_hypergeometric_test(sets, is.sig & sign > 0)
:END
:BEGIN down 
    run_hypergeometric_test(sets, is.sig & sign < 0)
:END
:BEGIN mixed
    run_hypergeometric_test(sets, is.sig)
:END
:BEGIN either
    up.stats <- run_hypergeometric_test(sets, is.sig & sign > 0)
    down.stats <- run_hypergeometric_test(sets, is.sig & sign < 0)
    up.stats$NumSig <- up.stats$NumSig + down.stats$NumSig
    direction <- up.stats$PValue < down.stats$PValue
    up.stats$PValue <- 2 * pmin(up.stats$PValue, down.stats$PValue)
    up.stats$FDR <- p.adjust(up.stats$PValue, method=\"BH\")
    up.stats$Direction <- ifelse(direction, \"up\", \"down\")
    up.stats
:END
})"

    parsed <- parseRmdTemplate(template)

    alternative <- match.arg(alternative)
    if (alternative != "mixed" && is.null(sign.name)) {
        stop("'sign' should be supplied when 'alternative=\"", alternative, "\"")
    }

    if (alternative == "up") {
        parsed$down <- NULL
        parsed$mixed <- NULL
        parsed$either <- NULL

    } else if (alternative == "down") {
        parsed$up <- NULL
        parsed$mixed <- NULL
        parsed$either <- NULL

    } else if (alternative == "either") {
        parsed$up <- NULL
        parsed$down <- NULL
        parsed$mixed <- NULL

    } else {
        parsed$up <- NULL
        parsed$down <- NULL
        parsed$either <- NULL
        parsed$`sign-setup` <- NULL
        parsed$`sign-filter` <- NULL
    }

    formatted <- replacePlaceholders(parsed, list(SETS=sets.name, IS_SIG=is.sig.name, SIGN=sign.name))
    unname(unlist(formatted, use.names=FALSE))
}
