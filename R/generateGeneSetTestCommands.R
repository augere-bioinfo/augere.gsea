#' @import augere.core
#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
generateGeneSetTestCommands <- function(sets.name, stat.name, alternative=c("mixed", "up", "down", "either"), seed=NULL, args=list()) {
    template <- "local({
    sets <- <%= SETS %>
    stat <- <%= STAT %>

    sizes <- unname(lengths(sets))
    keep <- sizes != 0L
    nonempty.sets <- unname(sets)[keep]

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
    replacements <- list(ALTERNATIVE = deparseToString(alternative), SETS=sets.name, STAT=stat.name)
    if (!is.null(seed)) {
        parsed$seed <- NULL
    } else {
        replacements$SEED <- deparseToString(seed)
    }

    if (alternative == "mixed") {
        # geneSetTest just takes the absolute value internally so any sign will automatically be ignored. 
        replacements$TYPE <- "\"f\""
    } else {
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
