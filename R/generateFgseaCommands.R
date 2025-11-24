#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges CharacterList 
.generateFgseaCommands <- function(sets.name, stat.name, seed, leading.edge, gene.names.name, alternative=c("mixed", "up", "down", "either"), args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    stat <- <%= STAT %>

    # Creating new names, just in case the existing names aren't unique.
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
    set.seed(<%= SEED %>)
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
:BEGIN leading-edge

    leading.genes <- <%= GENE_NAMES %>[match(unlist(raw$leadingEdge), gene.names)]
    leading.by <- factor(rep(raw$pathway, lengths(raw$leadingEdge)), set.names)
    leading.edge <- S4Vectors::splitAsList(leading.genes, leading.by)
    stopifnot(identical(names(leading.edge), set.names))
    out$LeadingEdge <- unname(leading.edge)
:END

    out 
})"

    parsed <- parseRmdTemplate(template)

    alternative <- match.arg(alternative)
    replacements <- list(SETS=sets.name, STAT=stat.name, SEED=deparseToString(seed))

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

    if (!leading.edge) {
        # I don't think much of the concept of the "leading edge", but here it is.
        # Honestly, if I wanted to know which genes contribute to enrichment, I'd prefer to rank all DEGs in the set and work my way down.
        # This ensures that I won't miss important genes that are significantly DE but not included in the leading edge because they're not at the very top of the ranking.
        # If I really needed an explicit subset, I'd just take all DE genes in the set; this is intuitive and has some statistical basis w.r.t. error control.
        parsed[["leading-edge"]] <- NULL
    } else {
        replacements$GENE_NAMES <- gene.names.name
    }

    if (length(args)) {
        parsed[["more-args"]] <- .format_more_args(args, comma.last=FALSE)
    }

    formatted <- replacePlaceholders(parsed, replacements)
    unname(unlist(formatted, use.names=FALSE))
}
