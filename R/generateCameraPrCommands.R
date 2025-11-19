#' Commands for pre-ranked CAMERA
#'
#' Generate R commands to run \code{\link[limma]{cameraPR}} under various parametrizations.
#'
#' @inheritParams generateHypergeometricTestCommands
#' @param args Named list of arguments to pass to \code{\link[limma]{cameraPR}}.
#'
#' @inheritSection generateGeneSetTestCommands Direction of the alternative hypothesis
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
#' cat(generateCameraPrCommands("sets", "stat"), sep="\n")
#' cat(generateCameraPrCommands("sets", "stat", alternative="either"), sep="\n")
#'
#' @export
#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateCameraPrCommands <- function(sets.name, stat.name, alternative=c("mixed", "up", "down", "either"), args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    stat <- <%= STAT %>

    old.names <- names(sets)
    names(sets) <- as.character(seq_along(sets))

    sizes <- unname(lengths(sets))
    keep <- sizes != 0L
    nonempty.sets <- sets[keep]

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
    replacements <- list(SETS=sets.name, STAT=stat.name)

    if (alternative == "mixed") {
        parsed[["one-sided"]] <- NULL

    } else {
        parsed[["de-sign"]] <- NULL
        parsed[["non-directional"]] <- NULL
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
