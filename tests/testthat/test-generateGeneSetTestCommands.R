# library(testthat); library(augere.gsea); source("test-generateGeneSetTestCommands.R")

set.seed(9999)

test_that("generateGeneSetTestCommands works as expected", {
    ngenes <- 1000
    sets <- list(
        A=1:20,
        B=50:100,
        C=200:300
    )

    cmds <- generateGeneSetTestCommands(sets.name="FOO", stat.name="BAR")

    env <- new.env()
    env$FOO <- sets
    env$BAR <- runif(ngenes, -1, 1)
    env$BAR[sets$B] <- rep(c(-10, 10), length.out=length(sets$B))

    result <- eval(parse(text=cmds), envir=env)
    expect_s4_class(result, "DFrame")
    expect_identical(colnames(result), c("NumGenes", "PValue", "FDR"))
    expect_identical(rownames(result), names(sets))

    expect_false(anyNA(result$FDR))
    expect_equal(result$NumGenes, c(20, 51, 101))

    expect_gt(result$PValue[1], 0.01)
    expect_lt(result$PValue[2], 0.01)
    expect_gt(result$PValue[3], 0.01)

    # p-values follow the enrichment trend correctly.
    env <- new.env()
    env$FOO <- sets
    env$BAR <- runif(ngenes, -1, 1)
    env$BAR[1:10] <- rep(c(-10, 10), length.out=10)
    env$BAR[71:80] <- rep(c(-10, 10), length.out=10)
    env$BAR[211:220] <- rep(c(-10, 10), length.out=10)

    result <- eval(parse(text=cmds), envir=env)
    expect_lt(result$PValue[1], result$PValue[2])
    expect_lt(result$PValue[2], result$PValue[3])

    # Further testing, this time with sets of the same size but different numbers of genes.
    env <- new.env()
    env$FOO <- list(alpha=1:50, bravo=101:150, charlie=201:250)
    env$BAR <- runif(ngenes, -1, 1)
    env$BAR[1:10] <- rep(c(-10, 10), length.out=10)
    env$BAR[111:130] <- rep(c(-10, 10), length.out=20)
    env$BAR[211:250] <- rep(c(-10, 10), length.out=40)

    result <- eval(parse(text=cmds), envir=env)
    expect_equal(result$NumGenes, rep(50, nrow(result)))
    expect_gt(result$PValue[1], result$PValue[2])
    expect_gt(result$PValue[2], result$PValue[3])
})

test_that("generateGeneSetTestCommands works with other alternative hypotheses", {
    ngenes <- 1000
    sets <- list(
        A=1:20,
        B=451:470,
        C=911:930
    )

    # Set A is only up, set B is only down, set C is 50:50.
    # Only the first half is significant.
    env <- new.env()
    env$FOO <- sets
    env$BAR <- runif(ngenes, -1, 1)
    env$BAR[1:10] <- 1
    env$BAR[451:460] <- -1
    env$BAR[911:920] <- rep(c(1, -1), length.out=10)

    cmds <- generateGeneSetTestCommands(sets.name="FOO", stat.name="BAR", alternative="up")
    result <- eval(parse(text=cmds), envir=env)
    expect_lt(result$PValue[1], result$PValue[3])
    expect_lt(result$PValue[3], result$PValue[2])

    cmds <- generateGeneSetTestCommands(sets.name="FOO", stat.name="BAR", alternative="down")
    result <- eval(parse(text=cmds), envir=env)
    expect_lt(result$PValue[2], result$PValue[3])
    expect_lt(result$PValue[3], result$PValue[1])

    cmds <- generateGeneSetTestCommands(sets.name="FOO", stat.name="BAR", alternative="either", seed=42)
    either.result <- eval(parse(text=cmds), envir=env)
    expect_equal(either.result$Direction[1], "up")
    expect_equal(either.result$Direction[2], "down")
    expect_lt(either.result$PValue[1], either.result$PValue[3])
    expect_lt(either.result$PValue[2], either.result$PValue[3])

    # Mixed should give us lower p-values than 'either' for the third set where the signs are mixed.
    cmds <- generateGeneSetTestCommands(sets.name="FOO", stat.name="BAR", alternative="mixed", seed=42)
    mixed.result <- eval(parse(text=cmds), envir=env)
    expect_lt(mixed.result$PValue[3], either.result$PValue[3])
})

test_that("generateGeneSetTestCommands works with empty sets", {
    ngenes <- 1000
    sets <- list(
        A=1:10,
        B=integer(0),
        C=990:1000
    )

    env <- new.env()
    env$FOO <- sets
    env$BAR <- runif(ngenes, -1, 1)

    cmds <- generateGeneSetTestCommands(sets.name="FOO", stat.name="BAR")
    result <- eval(parse(text=cmds), envir=env)
    expect_identical(which(is.na(result$PValue)), 2L)
})
