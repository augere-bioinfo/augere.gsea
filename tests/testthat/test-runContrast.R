# library(testthat); library(augere.gsea); source("test-runContrast.R")

ngenes <- 1000
set.seed(8888)
means <- runif(ngenes, 0, 10)
groups <- rep(LETTERS[1:3], each=4)

library(SummarizedExperiment)
counts <- matrix(rnbinom(ngenes * length(groups), mu=means, size=10), nrow=ngenes)
cd <- DataFrame(group=groups, block=rep(1:4, 3), age=seq_along(groups))
se <- SummarizedExperiment(list(counts=counts), colData=cd)
rownames(se) <- sprintf("gene-%i", seq_len(ngenes))

sets <- list(
    A = sample(rownames(se), 122),
    B = sample(rownames(se), 22),
    C = sample(rownames(se), 152),
    D = sample(rownames(se), 38),
    E = sample(rownames(se), 55)
)

output.default <- tempfile()
vanilla <- runContrast(se, sets, group="group", comparison=c("A", "C"), seed=100, output=output.default)

output.all <- tempfile()
all.methods <- c("mroast", "fry", "camera", "romer")
vanilla.all <- runContrast(se, sets, group="group", comparison=c("A", "C"), seed=100, save.results=FALSE, output=output.all, methods=all.methods)

test_that("runContrast works as expected", {
    expect_identical(names(vanilla), c("fry", "camera"))
    lines.def <- readLines(file.path(output.default, "report.Rmd"))
    expect_true(any(grepl("filterByExpr\\(.*group=", lines.def)))
    expect_true(any(grepl("voomWithQualityWeights", lines.def)))
    expect_true(any(grepl("robust=TRUE", lines.def)))
    expect_false(any(grepl("trend=FALSE", lines.def)))

    for (x in vanilla) {
        expect_s4_class(x, "DataFrame")
        expect_identical(rownames(x), names(sets))
        expect_type(x$NumGenes, "integer")
        expect_type(x$PValue, "double")
        expect_type(x$FDR, "double")
        expect_type(x$Direction, "character")
    }

    for (n in names(vanilla)) {
        roundtrip <- augere.core::readResult(file.path(output.default, "results", n))
        expect_identical(roundtrip$x, vanilla[[n]])
        expect_identical(n, roundtrip$metadata$differential_gene_set_test$method)
        expect_identical(roundtrip$metadata$differential_gene_set_test$contrast$left, list("A"))
        expect_identical(roundtrip$metadata$differential_gene_set_test$contrast$right, list("C"))
    }

    output.cov <- tempfile()
    covariates <- runContrast(se, sets, groups=NULL, comparison="age", covariate="age", seed=100, save.results=FALSE, output=output.cov)
    lines.cov <- readLines(file.path(output.cov, "report.Rmd"))
    expect_true(any(grepl("filterByExpr\\(.*design=", lines.cov)))
    expect_identical(names(vanilla), names(covariates))
    expect_false(identical(vanilla, covariates))

    output.block <- tempfile()
    blocked <- runContrast(se, sets, group="group", block="block", comparison=c("A", "C"), seed=100, save.results=FALSE, output=output.block)
    expect_identical(names(vanilla), names(blocked))
    expect_false(identical(vanilla, blocked))
})

test_that("runContrast works with subsetting", {
    output.nosub <- tempfile()
    nosub <- runContrast(se, sets, group="group", comparison=c("A", "C"), subset.groups=FALSE, seed=100, save.results=FALSE, output=output.nosub)
    expect_identical(names(nosub), names(vanilla))
    expect_false(identical(vanilla, nosub))

    output.subgroup <- tempfile()
    subgroup <- runContrast(
        se,
        sets,
        group="group",
        comparison=c("A", "C"),
        subset.groups=FALSE,
        subset.factor="group",
        subset.levels=c("A", "C"),
        seed=100,
        output=output.subgroup
    )
    expect_identical(vanilla, subgroup)

    roundtrip <- augere.core::readResult(file.path(output.subgroup, "results", "fry"))
    expect_identical(roundtrip$metadata$differential_gene_set_test$subset, "group IN (A,C)")
})

test_that("runContrast works for custom design and contrasts", {
    output.custom <- tempfile()
    custom <- runContrast(se, sets, design=~0 + group, contrast="groupB - groupA", seed=100, save.results=FALSE, output=output.custom)

    lines.custom <- readLines(file.path(output.custom, "report.Rmd"))
    expect_true(any(grepl("filterByExpr\\(.*design=", lines.custom)))

    # Comparing it to a simple setup.
    output.comp <- tempfile()
    comp <- runContrast(se, sets, group="group", comparison=c("B", "A"), subset.group=FALSE, seed=100, save.results=FALSE, output=output.comp)
    expect_identical(custom, comp)
})

test_that("runContrast works for duplicateCorrelation", {
    output.dupcor <- tempfile()
    dupcor <- runContrast(se, sets, group="group", dc.block="block", comparison=c("A", "C"), seed=100, save.results=FALSE, output=output.dupcor)

    lines.dupcor <- readLines(file.path(output.dupcor, "report.Rmd"))
    expect_true(any(grepl("duplicateCorrelation\\(", lines.dupcor)))

    expect_identical(names(vanilla), names(dupcor))
    expect_false(identical(vanilla, dupcor))
})

test_that("runContrast works with other limma-related options", {
    output.other <- tempfile()
    other <- runContrast(se, sets, group="group", comparison=c("A", "C"), quality=FALSE, trend=TRUE, robust=FALSE, seed=100, save.results=FALSE, output=output.other)

    lines.other <- readLines(file.path(output.other, "report.Rmd"))
    expect_true(any(grepl("trend=TRUE", lines.other)))
    expect_true(any(grepl("voom\\(", lines.other)))
    expect_false(any(grepl("robust=TRUE", lines.other)))

    expect_identical(names(vanilla), names(other))
    expect_false(identical(vanilla, other))
})

test_that("runContrast works with all or one methods", {
    output.one <- tempfile()
    one <- runContrast(se, sets, group="group", comparison=c("A", "C"), methods="camera", save.results=FALSE, output=output.one)
    expect_identical(names(one), "camera")

    all.methods <- c("mroast", "fry", "camera", "romer")
    output.all <- tempfile()
    all <- runContrast(se, sets, group="group", comparison=c("A", "C"), methods=all.methods, output=output.all)
    expect_identical(names(all), all.methods)
    expect_identical(all[names(vanilla)], vanilla)

    for (n in names(all)) {
        x <- all[[n]]
        expect_s4_class(x, "DataFrame")
        expect_identical(rownames(x), names(sets))
        expect_type(x$NumGenes, "integer")
        expect_type(x$PValue, "double")
        expect_type(x$FDR, "double")
        expect_type(x$Direction, "character")

        roundtrip <- augere.core::readResult(file.path(output.all, "results", n))
        expect_identical(roundtrip$x, all[[n]])
        expect_identical(n, roundtrip$metadata$differential_gene_set_test$method)
        expect_identical(roundtrip$metadata$differential_gene_set_test$contrast$left, list("A"))
        expect_identical(roundtrip$metadata$differential_gene_set_test$contrast$right, list("C"))
    }
})

test_that("runContrast works with disabled outputs", {
    output.disabled <- tempfile()
    disabled <- runContrast(se, sets, group="group", comparison=c("A", "C"), seed=100, save.results=FALSE, output=output.disabled)
    expect_false(file.exists(file.path(output.disabled, "results")))
    expect_true(file.exists(file.path(output.disabled, "report.Rmd")))
    expect_identical(names(disabled), names(vanilla))

    output.dry <- tempfile()
    dry <- runContrast(se, sets, group="group", comparison=c("A", "C"), seed=100, dry.run=TRUE, output=output.dry)
    expect_null(dry)
    expect_false(file.exists(file.path(output.disabled, "results")))
    expect_true(file.exists(file.path(output.disabled, "report.Rmd")))
})

test_that("runContrast works with customized annotation and metadata", {
    sets2 <- S4Vectors::List(sets)
    thing <- runif(length(sets))
    S4Vectors::mcols(sets2)$whee <- thing

    output.extra <- tempfile()
    extra <- runContrast(se, sets2, group="group", comparison=c("A", "C"), metadata=list(foo=1, bar=TRUE), annotation="whee", output=output.extra)
    expect_identical(names(extra), names(vanilla))

    for (n in names(extra)) {
        expect_identical(extra[[n]]$whee, thing)
        expect_type(extra[[n]]$FDR, "double")

        roundtrip <- augere.core::readResult(file.path(output.extra, "results", n))
        expect_identical(roundtrip$metadata$foo, 1L)
        expect_true(roundtrip$metadata$bar)
    }
})
