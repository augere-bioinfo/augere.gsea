# library(testthat); library(augere.gsea); source("test-runPrecomputed.R") 

set.seed(99999)
all.genes <- sprintf("gene-%i", seq_len(1000))

library(S4Vectors)
tab <- DataFrame(
    row.names=all.genes,
    AveExpr=rnorm(length(all.genes)), 
    t=rnorm(length(all.genes)), 
    PValue=runif(length(all.genes)),
    FDR=runif(length(all.genes))
)

sets <- list(
    A=sample(all.genes, 20),
    B=sample(all.genes, 300),
    C=sample(all.genes, 10),
    D=sample(all.genes, 500)
)

# Vanilla unsigned run to get started.
all.methods <- c("hypergeometric", "goseq", "geneSetTest", "fgsea", "cameraPR")
out <- tempfile()
stats <- runPrecomputed(tab, sets, methods=all.methods, output.dir=out, seed=42)

test_that("runPrecomputed works as expected", {
    expect_identical(names(stats), all.methods)
    expect_true(file.exists(file.path(out, "report.Rmd")))

    for (meth in all.methods) {
        df <- stats[[meth]]
        expect_s4_class(df, "DFrame")
        expect_identical(rownames(df), names(sets))
        expect_type(df$NumGenes, "integer")
        expect_type(df$FDR, "double")

        res <- augere.core::readResult(file.path(out, "results", meth))
        expect_identical(res$metadata$precomputed_gene_set_enrichment$method, meth)
        expect_identical(rownames(res$x), names(sets))
        expect_equal(res$x$PValue, df$PValue)
    }
})

test_that("runPrecomputed works with the various sign options", {
    tab2 <- tab
    tab2$t.abs <- abs(tab2$t)
    tab2$LogFC <- sign(tab2$t) * runif(nrow(tab2))

    # Adding some signedness.
    out2 <- tempfile()
    abs.stats <- runPrecomputed(tab2, sets, methods=all.methods, rank.field="t.abs", sign.field="LogFC", output.dir=out2, save.results=FALSE, seed=42)
    for (meth in all.methods) {
        expect_equal(stats[[meth]], abs.stats[[meth]])
    }

    # With some square rooting.
    out2 <- tempfile()
    tab2$F <- tab2$t^2
    sqrt.stats <- runPrecomputed(tab2, sets, methods=all.methods, rank.field="F", sign.field="LogFC", rank.sqrt=TRUE, output.dir=out2, save.results=FALSE, seed=42)
    for (meth in all.methods) {
        expect_equal(stats[[meth]], sqrt.stats[[meth]])
    }
})

test_that("runPrecomputed works with different alternatives", {
    out2 <- tempfile()
    expect_error(runPrecomputed(tab, sets, alternative="up", output.dir=out2, save.results=FALSE, seed=42), "'sign.field' should be supplied")

    out2 <- tempfile()
    rank.only <- runPrecomputed(tab, sets, alternative="up", methods=c("fgsea", "cameraPR"), output.dir=out2, seed=42)
    expect_identical(names(rank.only), c("fgsea", "cameraPR"))
    for (n in names(rank.only)) {
        res <- augere.core::readResult(file.path(out2, "results", n))
        expect_identical(res$metadata$precomputed_gene_set_enrichment$alternative, "up")
    }

    out2 <- tempfile()
    hyp.only <- runPrecomputed(tab, sets, alternative="down", sign.field="t", methods=c("hypergeometric", "goseq"), output.dir=out2, seed=42)
    expect_identical(names(hyp.only), c("hypergeometric", "goseq"))
    for (n in names(hyp.only)) {
        res <- augere.core::readResult(file.path(out2, "results", n))
        expect_identical(res$metadata$precomputed_gene_set_enrichment$alternative, "down")
    }
})

test_that("runPrecomputed works with only a subset of methods", {
    # No need for rank.
    simple.out <- tempfile()
    simple.methods <- c("hypergeometric", "goseq")
    simple.stats <- runPrecomputed(tab, sets, methods=simple.methods, rank.field=NULL, output.dir=simple.out, save.results=FALSE, seed=42)
    expect_identical(names(simple.stats), simple.methods)
    for (meth in simple.methods) {
        expect_equal(stats[[meth]], simple.stats[[meth]])
    }

    # Now rank-only.
    ro.out <- tempfile()
    ro.methods <- setdiff(all.methods, simple.methods)
    ro.stats <- runPrecomputed(tab, sets, methods=ro.methods, signif.field=NULL, output.dir=ro.out, save.results=FALSE, seed=42)
    expect_identical(names(ro.stats), ro.methods)
    for (meth in ro.methods) { # fortunately the seed is still consistent.
        expect_equal(stats[[meth]], ro.stats[[meth]])
    }
})

test_that("runPrecomputed works with the different output options", {
    # No saving.
    nosave.out <- tempfile()
    nosave.stats <- runPrecomputed(tab, sets, methods=all.methods, save.results=FALSE, output.dir=nosave.out, seed=42)
    expect_equal(stats, nosave.stats)
    expect_true(file.exists(file.path(nosave.out, "report.Rmd")))
    expect_false(file.exists(file.path(nosave.out, "results")))

    # Dry run only.
    dry.out <- tempfile()
    dry.stats <- runPrecomputed(tab, sets, methods=all.methods, save.results=FALSE, dry.run=TRUE, output.dir=dry.out, seed=42)
    expect_null(dry.stats)
    expect_true(file.exists(file.path(dry.out, "report.Rmd")))
    expect_false(file.exists(file.path(dry.out, "results")))
})

test_that("runPrecomputed respects some extra metadata", {
    # Metadata to write to disk. 
    common.out <- tempfile()
    common.stats <- runPrecomputed(tab, sets, methods=all.methods, output.dir=common.out, metadata=list(foo=1, bar=TRUE), seed=69)
    for (meth in all.methods) {
        res <- augere.core::readResult(file.path(common.out, "results", meth))
        expect_identical(res$metadata$precomputed_gene_set_enrichment$method, meth)
        expect_identical(res$metadata$foo, 1L)
        expect_true(res$metadata$bar)
    }

    # Some extra annotations.
    out <- tempfile()
    lsets <- List(sets)
    mcols(lsets)$whee <- runif(length(lsets))
    stats <- runPrecomputed(tab, lsets, methods=all.methods, output.dir=out, annotation="whee", save.results=FALSE, seed=69)
    for (meth in all.methods) {
        df <- stats[[meth]]
        expect_type(df$NumGenes, "integer")
        expect_type(df$FDR, "double")
        expect_identical(df$whee, mcols(lsets)$whee)
    }
})
