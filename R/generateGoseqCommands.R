#' Commands for goseq 
#'
#' Generate R commands to run \code{\link[goseq]{goseq}} under various parameterizations.
#' 
#' @inheritParams generateHypergeometricTestCommands
#' @param bias.name String containing the variable name for the numeric vector of biases (typically lengths or abundances) for each gene.
#' The referenced vector should be of length equal to the number of genes.
#'
#' @return Character vector containing the commands required to run \code{\link[goseq]{goseq}}.
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
#' @inheritSection generateHypergeometricTestCommands Direction of the alternative hypothesis
#'
#' @author Aaron Lun
#' 
#' @examples
#' cat(generateGoseqCommands("sets", "is.sig", "sign", "logcpm"), sep="\n")
#' cat(generateGoseqCommands("sets", "is.sig", "sign", "logcpm", alternative="either"), sep="\n")
#'
#' @export
#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateGoseqCommands <- function(sets.name, is.sig.name, sign.name, bias.name, alternative=c("mixed", "up", "down", "either"), args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    is.sig <- <%= IS_SIG %>
    bias <- <%= BIAS %>
:BEGIN sign-setup
    sign <- <%= SIGN %>
:END

    gene.names <- as.character(seq_along(is.sig))
    names(is.sig) <- names(bias) <- gene.names
    named.sets <- lapply(sets, function(x) gene.names[x])

    set.names <- as.character(seq_along(named.sets))
    names(named.sets) <- set.names
    sizes <- unname(lengths(named.sets))
    keep <- sizes != 0L
    named.sets <- named.sets[keep]

    # Make sure the bias is non-negative.
    min.bias <- min(bias)
    if (min.bias < 0) {
        bias <- bias - min.bias
    }

    run_goseq <- function(is.sig) {
        pwf <- goseq::nullp(is.sig, bias.data=bias)
        raw <- goseq::goseq(
            pwf,
            gene2cat=named.sets,
            use_genes_without_cat=TRUE, 
:BEGIN more-args
:END
        )

        m <- match(set.names, raw$category)
        out <- DataFrame(
            row.names=names(sets),
            NumGenes=sizes,
            NumSig=raw$numDEInCat[m],
            PValue=raw$over_represented_pvalue[m],
            FDR=p.adjust(out$PValue, method=\"BH\")
        )
        out$NumSig[!keep] <- 0L
        out 
    }

:BEGIN mixed
    run_goseq(is.sig)
:END
:BEGIN up
    run_goseq(is.sig & sign > 0)
:END
:BEGIN down
    run_goseq(is.sig & sign < 0)
:END
:BEGIN either
    up.stats <- run_goseq(is.sig & sign > 0)
    up.stats <- run_goseq(is.sig & sign < 0)
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
    if (alternative != "mixed" && is.null(sign)) {
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

    formatted <- replacePlaceholders(parsed, list(SETS=sets.name, IS_SIG=is.sig.name, SIGN=sign.name, BIAS=bias.name))
    unname(unlist(formatted, use.names=FALSE))
}
