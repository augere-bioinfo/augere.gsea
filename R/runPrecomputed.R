#' GSEA on precomputed statistics
#'
#' Run gene set enrichment analyses on precomputed statistics, usually from differential expression analyses.
#' This uses competitive gene set tests where the enrichment of \dQuote{interesting} genes within the set must be greater than that outside of the set.
#' 
#' @param x A data frame or \link[S4Vectors]{DataFrame} object containing various test statistics (columns) for genes (rows).
#' Rows should be named with the same gene identifiers used in \code{sets}.
#' Rows may contain \code{NA} values, which will be removed prior to analysis.
#'
#' Rows should \emph{not} be filtered to only the significant genes, this will be handled by \code{signif.threshold}.
#' Users should \strong{not} pass in a table containing only the significant results. 
#' 
#' Alternatively, an object returned by \link[augere.core]{wrapInput} that refers to a data frame or DataFrame.
#' @param sets A list or \link[IRanges]{CharacterList} of character vectors.
#' Each vector represents a gene set and contains the identifiers for genes in that set.
#' Identifiers should be consistent with the row names of \code{x}.
#' The list may also be named with the names of the sets.
#' @param methods Character vector specifying the gene set testing methods to use.
#' This can be any number of the following:
#' \itemize{
#' \item \code{"hypergeometric"}, a hypergeometric test for enrichment of significant genes in each set.
#' The set of significant genes is determined by \code{signif.field} and \code{signif.threshold}.
#' This is equivalent to a one-sided Fisher's exact test.
#' \item \code{"goseq"}, which uses \code{\link[goseq]{goseq}} from the \pkg{goseq} package.
#' This is much like the hypergeometric test but accounting for gene-specific biases in detection.
#' \item \code{"geneSetTest"}, which uses \code{\link[limma]{geneSetTest}} from the \pkg{limma} package.
#' This ranks genes according to some statistic (see \code{rank.field}) and tests for differences in the ranking of genes in and outside of the set.
#' \item \code{"fgsea"}, which uses \code{\link[fgsea]{fgsea}} from the \pkg{fgsea} package.
#' This also tests for differences in ranking but is faster than \code{geneSetTest}.
#' \item \code{"cameraPR"}, which uses \code{\link[limma]{cameraPR}} from the \pkg{limma} package.
#' This also tests for differences in ranking while accounting for variance inflation due to correlations between genes.
#' }
#' @param alternative String specifying the alternative hypothesis.
#' This should be one of:
#' \itemize{
#' \item \code{"mixed"} tests for enrichment of any significant genes in each set, regardless of their direction.
#' \item \code{"up"} tests for enrichment of up-regulated genes in each set.
#' \item \code{"down"} tests for enrichment of down-regulated genes in each set.
#' \item \code{"either"} tests for enrichment of either up- or down-regulated genes in each set.
#' Specifically, each gene set is tested separately for enrichment of up- and down-regulated genes, and the results are combined into a single p-value.
#' This differs from \code{"mixed"}, which does not consider the sign at all and only cares about enrichment of significant genes of any direction in each set.
#' }
#' @param signif.field String specifying the column of \code{x} to define significant genes, typically the adjusted p-value.
#' Only used for \code{method="hypergeometric"} and \code{"goseq"}.
#' @param signif.threshold Number specifying the upper threshold for significance to apply to the statistics from \code{signif.field}.
#' All genes with lower statistics are considered to be significant.
#' Only used for \code{method="hypergeometric"} and \code{"goseq"}.
#' @param rank.field String specifying the column of \code{x} containing test statistics for ranking, e.g., t-statistics, Z-scores.
#' This is generally expected to be signed such that the values with the largest magnitude are most significant and the sign represents some kind of directionality.
#' It may be unsigned if \code{sign.field} is provided, in which case the signs are used to convert the test statistics back to signed values;
#' or if \code{alternative="mixed"}, in which case the signs are ignored by each method.
#' Only used for \code{method="geneSetTest"}, \code{"fgsea"} and \code{"cameraPR"}.
#' @param sign.field String specifying the column of \code{x} containing a signed effect size, e.g., the log-fold change.
#' This is used to restore the sign to unsigned test statistics in \code{rank.field}.
#' It will also be used to define up- or down-regulated genes in \code{method="hypergeometric"} and \code{"goseq"} when \code{alternative} is not \code{"mixed"}.
#' @param rank.sqrt Boolean indicating whether to compute the square root of the \code{rank.field} statistic before restoring the sign with \code{sign.field}. 
#' For example, the F-statistic will be converted to a t-statistic while the likelihoiod ratio will be converted to a Z-score.
#' Only used for \code{method="geneSetTest"}, \code{"fgsea"} and \code{"cameraPR"}.
#' @param goseq.bias String specifying the column of \code{x} containing the per-gene detection bias.
#' Defaults to the gene abundance as larger counts provide greater power for detecting differential expression,
#' but if this is not available, the exonic gene length can also be used.
#' Only used for \code{method="goseq"}.
#' @param goseq.args Named list of additional arguments to pass to \code{\link[goseq]{goseq}}.
#' @param geneSetTest.args Named list of additional arguments to pass to \code{\link[limma]{geneSetTest}}.
#' @param cameraPR.args Named list of additional arguments to pass to \code{\link[limma]{cameraPR}}.
#' @param fgsea.leading.edge Boolean indicating whether the \dQuote{leading edge} should be stored from \code{\link[fgsea]{fgsea}}. 
#' @param fgsea.args Named list of additional arguments to pass to \code{\link[fgsea]{fgsea}}.
#' @param output.dir String containing the path to an output directory in which to write the Rmarkdown file and save results.
#' @param metadata Named list of additional metadata to store with each result.
#' @param author Character vector containg the names of the authors.
#' If \code{NULL}, defaults to the current user.
#' @param annotation Character vector specifying the columns of \code{\link[S4Vectors]{mcols}(sets)} to store in each result DataFrames.
#' @param dry.run Logical scalar indicating whether a dry run should be performed,
#' This generates the Rmarkdown report in \code{output.dir} but does not execute the analysis.
#' @param save.results Boolean indicating whether the results should be saved to file.
#'
#' @details
#' Some of the methods involve randomization, so for full reproducibility, users should call \code{\link{set.seed}} before running \code{runPrecomputed}.
#'
#' Note that, even if the user does not call \code{\link{set.seed}},
#' \code{runPrecomputed} will automatically insert \code{\link{set.seed}} statements into the Rmarkdown report 
#' prior to any GSEA functions that involve randomization.
#' Each \code{set.seed} call has a hard-coded seed to ensure that future compilation of the generated report will give the same result.
#' However, different calls to \code{runPrecomputed} will use different (randomly selected) seeds to avoid systematic biases.
#' Thus, if full reproducibility of \code{runPrecomputed} is required, users should set the seed themselves before calling \code{runPrecomputed}.
#'
#' @return
#' A Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir}. 
#' This contains all commmands used to reproduce the analysis. 
#'
#' If \code{dry.run=FALSE}, a list of \link[S4Vectors]{DataFrame}s is returned where each DataFrame contains the enrichment results for a \code{method}.
#' Each row corresponds to a gene set in \code{sets} and is named accordingly.
#' All DataFrames are guaranteed to have (at least) the following fields:
#' \itemize{
#' \item \code{NumGenes}, the number of genes in the set (after removing genes that were not tested).
#' \item \code{PValue}, the p-value for enrichment in each gene set.
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value. 
#' \item (if \code{alternative="either"}) \code{Direction}, whether the set is more enriched for \code{"up"}- or \code{"down"}-regulated genes.
#' }
#' Specific methods will have additional fields:
#' \itemize{
#' \item For \code{"fgsea"}, \code{ES} contains the enrichment score and \code{NES} contains the normalized enrichment score.
#' If \code{fgsea.leading.edge=TRUE}, \code{LeadingEdge} will contain a \link[IRanges]{CharacterList} with the names of genes in the \dQuote{leading edge} for each gene set.
#' \item For \code{"goseq"} and \code{"hypergeometric"}, \code{NumSig} contains the number of significant genes in each set.
#' (If \code{alternative} is \code{"up"} or \code{"down"}, this is filtered for significant genes of the desired sign.)
#' }
#' 
#' If \code{save.results=TRUE}, the results are saved in a \code{results} subdirectory of \code{output.dir}.
#'
#' If \code{dry.run=FALSE}, only the report is created, and \code{NULL} is returned.
#' 
#' @author Aaron Lun
#' 
#' @examples
#' all.genes <- sprintf("gene-%i", seq_len(1000))
#' 
#' library(S4Vectors)
#' tab <- DataFrame(
#'     row.names=all.genes,
#'     AveExpr=rnorm(length(all.genes)), 
#'     t=rnorm(length(all.genes)), 
#'     PValue=runif(length(all.genes)),
#'     FDR=runif(length(all.genes))
#' )
#' 
#' sets <- list(
#'     A=sample(all.genes, 20),
#'     B=sample(all.genes, 300),
#'     C=sample(all.genes, 10),
#'     D=sample(all.genes, 500)
#' )
#'
#' output.dir <- tempfile()
#' results <- runPrecomputed(
#'     tab,
#'     sets,
#'     signif.field="FDR",
#'     signif.threshold=0.05,
#'     rank.field="t",
#'     output=output.dir
#' )
#' results$hypergeometric
#' results$fgsea
#'
#' list.files(output.dir, recursive=TRUE)
#' 
#' @export
#' @import augere.core
runPrecomputed <- function(
    x,
    sets,
    signif.field, 
    signif.threshold, 
    rank.field,
    methods = c("hypergeometric", "goseq", "fgsea", "cameraPR"),
    alternative = c("mixed", "up", "down", "either"),
    rank.sqrt = FALSE,
    sign.field = NULL,
    goseq.bias = "AveExpr",
    goseq.args = list(),
    fgsea.leading.edge = FALSE,
    fgsea.args = list(),
    geneSetTest.args = list(),
    cameraPR.args = list(),
    metadata = NULL,
    annotation = NULL,
    author = NULL,
    output.dir = "precomputed", 
    dry.run = FALSE, 
    save.results = TRUE
) {
    restore.cache <- resetInputCache()
    on.exit(restore.cache(), after=FALSE, add=TRUE)

    if (is.null(author)) {
        author <- Sys.info()[["user"]]
    }

    dir.create(output.dir, showWarnings=FALSE, recursive=TRUE)
    fname <- file.path(output.dir, "report.Rmd")

    template <- system.file("templates", "precomputed.Rmd", package="augere.gsea")
    lines <- readLines(template)
    parsed <- parseRmdTemplate(lines)

    parsed[["setup-tab"]] <- processInputCommands(x, name="tab")
    parsed[["setup-sets"]] <- processInputCommands(sets, name="sets")

    alternative <- match.arg(alternative)
    replacements <- list(
        AUTHOR=paste(sprintf("  - %s", author), collapse="\n"),
        RAW_AUTHOR=deparseToString(as.list(author)),
        SIGN_FIELD=deparseToString(sign.field),
        ALTERNATIVE=deparseToString(alternative)
    )

    fields.in.use <- character(0)
    if (any(c("hypergeometric", "goseq") %in% methods)) {
        replacements$SIGNIF_FIELD <- deparseToString(signif.field)
        replacements$SIGNIF_THRESHOLD <- deparseToString(signif.threshold)
        fields.in.use <- c(fields.in.use, signif.field)
    } else {
        parsed[["setup-signif"]] <- NULL
    }

    if (any(c("fgsea", "cameraPR", "romer") %in% methods)) {
        replacements$RANK_FIELD <- deparseToString(rank.field)
        fields.in.use <- c(fields.in.use, rank.field)
    } else {
        parsed[["setup-rank"]] <- NULL
    }

    # We use the sign for the rank-based methods if it's provided, otherwise we assume that the rank stat is already sufficiently signed.
    use.sign <- any(c("fgsea", "cameraPR", "romer") %in% methods) && !is.null(sign.field)
    if (any(c("hypergeometric", "goseq") %in% methods) && alternative != "mixed") {
        # For hypergeometric/goseq, we actually need the sign if we're not testing the mixed alternative hypothesis.
        if (is.null(sign.field)) {
            stop("'sign.field' should be supplied when 'alternative=\"", alternative, "\"")
        }
        use.sign <- TRUE
    }
    if (use.sign) {
        replacements$SIGN_FIELD <- deparseToString(sign.field)
        if (rank.sqrt) {
            replacements$SANITIZER <- "sqrt(abs(rank.stat))"
        } else {
            parsed[["sqrt-txt"]] <- NULL
            replacements$SANITIZER <- "abs(rank.stat)"
        }
        fields.in.use <- c(fields.in.use, sign.field)
    } else {
        parsed[["re-sign"]] <- NULL
        parsed[["setup-sign"]] <- NULL
    }

    replacements$ALL_STAT_NAMES <- deparseToString(fields.in.use)

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        replacements$COMMON_METADATA <- deparseToString(metadata)
    } else {
        parsed[["create-common-metadata"]] <- NULL
    }

    if (!is.null(annotation)) {
        replacements$ANNO_FIELDS <- deparseToString(annotation)
    } else {
        parsed[["common-annotation"]] <- NULL
    }

    methods <- match.arg(methods, c("hypergeometric", "goseq", "geneSetTest", "cameraPR", "fgsea"), several.ok=TRUE)
    save.chunk.names <- character()
    rseed <- function() sample(.Machine$integer.max, 1)

    process_common <- function(y) {
        if (is.null(annotation)) {
            y[["add-annotation"]] <- NULL
        }
        if (!merge.metadata) {
            y[["merge-metadata"]] <- NULL
        } else {
            y[["no-merge-metadata"]] <- NULL
        }
        y
    }

    if ("hypergeometric" %in% methods) {
        replacements$HYPER_CMDS <- paste(.generateHypergeometricTestCommands("indices", "is.sig", "sign.stat", alternative=alternative), collapse="\n")
        parsed$hypergeometric <- process_common(parsed$hypergeometric)
        save.name <- "save-hyper"
        save.chunk.names <- c(save.chunk.names, save.name)
        replacements$HYPER_SAVE_CHUNK_NAME <- save.name
    } else {
        parsed$hypergeometric <- NULL
    }

    if ("goseq" %in% methods) {
        replacements$GOSEQ_CMDS <- paste(.generateGoseqCommands("indices", "is.sig", "sign.stat", "bias.stat", alternative=alternative, args=goseq.args), collapse="\n")
        replacements$GOSEQ_BIAS_FIELD <- deparseToString(goseq.bias)
        parsed$goseq <- process_common(parsed$goseq)
        save.name <- "save-goseq"
        save.chunk.names <- c(save.chunk.names, save.name)
        replacements$GOSEQ_SAVE_CHUNK_NAME <- save.name
    } else {
        parsed$goseq <- NULL
    }

    if ("geneSetTest" %in% methods) { 
        replacements$GST_CMDS <- paste(.generateGeneSetTestCommands("indices", "rank.stat", seed=rseed(), alternative=alternative, args=geneSetTest.args), collapse="\n")
        parsed$geneSetTest <- process_common(parsed$geneSetTest)
        save.name <- "save-gst"
        save.chunk.names <- c(save.chunk.names, save.name)
        replacements$GST_SAVE_CHUNK_NAME <- save.name
    } else {
        parsed$geneSetTest <- NULL
    }

    if ("fgsea" %in% methods) {
        replacements$FGSEA_CMDS <- paste(
            .generateFgseaCommands(
                sets.name="indices",
                stat.name="rank.stat",
                seed=rseed(),
                leading.edge=fgsea.leading.edge,
                alternative=alternative,
                gene.names.name="rownames(tab)",
                args=fgsea.args
            ),
            collapse="\n"
        )
        parsed$fgsea <- process_common(parsed$fgsea)
        save.name <- "save-fgsea"
        save.chunk.names <- c(save.chunk.names, save.name)
        replacements$FGSEA_SAVE_CHUNK_NAME <- save.name
    } else {
        parsed$fgsea <- NULL
    }

    if ("cameraPR" %in% methods) {
        replacements$CAMERA_CMDS <- paste(.generateCameraPrCommands("indices", "rank.stat", alternative=alternative, args=cameraPR.args), collapse="\n")
        parsed$cameraPR <- process_common(parsed$cameraPR)
        save.name <- "save-cameraPR"
        save.chunk.names <- c(save.chunk.names, save.name)
        replacements$CAMERA_SAVE_CHUNK_NAME <- save.name
    } else {
        parsed$cameraPR <- NULL
    }

    contents <- replacePlaceholders(parsed, replacements)
    writeRmd(contents, fname)

    if (dry.run) {
        return(NULL)
    }

    if (save.results) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- c("save-directory", save.chunk.names)
    }

    env <- new.env()
    compileReport(fname, env=env, skip.chunks=skip.chunks)
    env$all.results
}
