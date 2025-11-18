#' Commands for pre-ranked CAMERA
#'
#' Generate R commands to run \code{\link[limma]{cameraPR}} under various parametrizations.
#'
#' @inheritParams generateHypergeometricTestCommands
#' @param args Named list of arguments to pass to \code{\link[limma]{cameraPR}}.
#'
#' @inheritSection generateGeneSetTestCommands Direction of the alternative hypothesis
#'
#' @details
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
#' @return Character vector containing the commands required to run \code{cameraPR}.
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
#' cat(generateCameraPrCommands("sets", "stat", "sign"), sep="\n")
#' cat(generateCameraPrCommands("sets", "stat", "sign", alternative="either"), sep="\n")
#'
#' @export
#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateCameraPrCommands <- function(sets.name, stat.name, sign.name, alternative=c("mixed", "up", "down", "either"), use.sqrt=FALSE, args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    stat <- <%= STAT %>
:BEGIN sign-setup
    sign <- <%= SIGN %>
:END

    old.names <- names(sets)
    names(sets) <- as.character(seq_along(sets))

    sizes <- unname(lengths(sets))
    keep <- sizes != 0L
    nonempty.sets <- sets[keep]

:BEGIN re-sign
    # Re-introducing the sign.
    stat <- <%= SANITIZER %> * base::sign(sign)

:END
:BEGIN de-sign
    # Removing the sign to avoid warnings.
    stat <- abs(stat)

:END
    out <- limma::cameraPR(
        statistic=stat,
        index=nonempty.sets,
:BEGIN non-directional 
        use.ranks=TRUE,
        directional=FALSE,
:END
        sort=FALSE,
:BEGIN more-args
:END
    )

:BEGIN one-sided
    # Converting them to one-sided p-values.
    out$PValue <- out$PValue / 2
    out$PValue <- ifelse(out$Direction == <%= DIRECTION %>, out$PValue, 1 - out$PValue)
    out$Direction <- NULL

:END
    out <- S4Vectors::DataFrame(out)
    colnames(out)[colnames(out) == \"NGenes\"] <- \"NumGenes\"
    out$FDR <- p.adjust(out$PValue, method='BH')
    out <- out[match(names(sets), rownames(out)),]
:BEGIN report-direction 
    out$Direction <- tolower(out$Direction)
:END
    rownames(out) <- old.names
    out
})"

    parsed <- parseRmdTemplate(template)

    alternative <- match.arg(alternative)
    replacements <- list(SETS=sets.name, STAT=stat.name, SIGN=sign.name)

    if (alternative == "mixed") {
        parsed[["re-sign"]] <- NULL
        parsed[["sign-setup"]] <- NULL
        parsed[["one-sided"]] <- NULL

    } else {
        parsed[["de-sign"]] <- NULL
        parsed[["non-directional"]] <- NULL

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

        if (alternative == "either") {
            parsed[["one-sided"]] <- NULL
        } else {
            replacements$DIRECTION <- deparseToString(if (alternative == "up") "Up" else "Down")
        }
    }

    if (alternative != "either") {
        parsed[["report-direction"]] <- NULL
    }

    if (length(args)) {
        parsed[["more-args"]] <- .format_more_args(args, comma.last=FALSE)
    }

    formatted <- replacePlaceholders(parsed, replacements)
    unname(unlist(formatted, use.names=FALSE))
}
