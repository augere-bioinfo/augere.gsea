#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
.generateGoseqCommands <- function(sets.name, is.sig.name, sign.name, bias.name, alternative = "mixed", args = list()) {
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
        out <- S4Vectors::DataFrame(
            row.names=names(sets),
            NumGenes=sizes,
            NumSig=raw$numDEInCat[m],
            PValue=raw$over_represented_pvalue[m]
        )
        out$NumSig[!keep] <- 0L
        out$FDR <- p.adjust(out$PValue, method=\"BH\")
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
    down.stats <- run_goseq(is.sig & sign < 0)
    up.stats$NumSig <- up.stats$NumSig + down.stats$NumSig
    direction <- up.stats$PValue < down.stats$PValue
    up.stats$PValue <- 2 * pmin(up.stats$PValue, down.stats$PValue)
    up.stats$FDR <- p.adjust(up.stats$PValue, method=\"BH\")
    up.stats$Direction <- ifelse(direction, \"up\", \"down\")
    up.stats
:END
})"

    parsed <- parseRmdTemplate(template)

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
