#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
.generateHypergeometricTestCommands <- function(sets.name, is.sig.name, sign.name, alternative = "mixed") {
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
