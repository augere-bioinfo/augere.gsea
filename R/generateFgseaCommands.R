#' Commands for FGSEA
#'
#' Generate R commands to run \code{\link[fgsea]{fgsea}} under various parametrizations.
#'
#' @inheritParams generateGeneSetTestCommands
#' @param args Named list of arguments to pass to \code{\link[fgsea]{fgsea}}.
#'
#' @inheritSection generateGeneSetTestCommands Direction of the alternative hypothesis
#'
#' @author Aaron Lun
#'
#' @return Character vector containing the commands required to run \code{fgsea}.
#' Upon evaluation, this produces a \link[S4Vectors]{DFrame} with one row per gene set and the following columns:
#' \itemize{
#' \item \code{NumGenes}, the number of genes in each set.
#' \item \code{ES}, the enrichment score.
#' \item \code{NES}, the normalized enrichment score.
#' \item \code{PValue}, the p-value for enrichment of highly-ranked genes in this set.
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value.
#' \item \code{Direction}, whether the set is more enriched for up- or down-regulated genes.
#' Only reported for \code{alternative="either"}.
#' }
#'
#' @examples
#' cat(generateFgseaCommands("sets", "stat"), sep="\n")
#' cat(generateFgseaCommands("sets", "stat", alternative="either"), sep="\n")
#'
#' @export
#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateFgseaCommands <- function(sets.name, stat.name, alternative=c("mixed", "up", "down", "either"), seed=NULL, args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    stat <- <%= STAT %>

    gene.names <- as.character(seq_along(stat))
    named.sets <- relist(gene.names[unlist(sets)], sets)
    names(stat) <- gene.names
    set.names <- as.character(seq_along(named.sets))
    names(named.sets) <- set.names

    sizes <- unname(lengths(named.sets))
    keep <- sizes != 0L
    named.sets <- named.sets[keep]

:BEGIN de-sign
    stat <- abs(stat) # Removing the sign so that fgsea only uses the magnitude.
:END
:BEGIN seed
    set.seed(<%= SEED %>)
:END
    raw <- fgsea::fgsea(
        named.sets,
        stat,
        scoreType=<%= SCORE_TYPE %>,
:BEGIN more-args
:END
    )

    out <- S4Vectors::DataFrame(
         row.names=raw$pathway,
         ES=raw$ES, 
         NES=raw$NES, 
         PValue=raw$pval, 
         FDR=raw$padj
    )

    out <- out[match(set.names, rownames(out)),,drop=FALSE]
    out <- cbind(S4Vectors::DataFrame(NumGenes=sizes), out)
    rownames(out) <- names(sets)
:BEGIN direction
    out$Direction <- ifelse(out$NES > 0, \"up\", \"down\")
:END
    out 
})"

    parsed <- parseRmdTemplate(template)

    alternative <- match.arg(alternative)
    replacements <- list(SETS=sets.name, STAT=stat.name)

    if (alternative == "mixed") {
        # FYI check out fgsea:::preparePathwaysAndStats where they just take
        # the absolute value for the calculations.
        replacements$SCORE_TYPE <- "\"pos\""
        parsed[["direction"]] <- NULL

    } else {
        replacements$SCORE_TYPE <- deparseToString(switch(alternative, up = "pos", down = "neg", either = "std"))
        parsed[["de-sign"]] <- NULL
        if (alternative != "either") {
            parsed[["direction"]] <- NULL
        }
    }

    if (length(args)) {
        parsed[["more-args"]] <- .format_more_args(args, comma.last=FALSE)
    }
    if (!is.null(seed)) {
        replacements$SEED <- seed
    } else {
        parsed[["seed"]] <- NULL
    }

    formatted <- replacePlaceholders(parsed, replacements)
    unname(unlist(formatted, use.names=FALSE))
}
