# library(testthat); library(augere.gsea); source("test-generateHypergeometricTestCommands.R")

test_that("generateHypergeometricTestCommands works as expected", {
    ngenes <- 1000
    sets <- list(
        A=1:20,
        B=50:100,
        C=200:300
    )

    cmds <- augere.gsea:::.generateHypergeometricTestCommands(sets.name="FOO", is.sig.name="BAR", sign.name=NULL)
    expect_false(any(grepl("sign <-", cmds)))

    env <- new.env()
    env$FOO <- sets
    env$BAR <- logical(ngenes)
    env$BAR[sets$B] <- TRUE

    result <- eval(parse(text=cmds), envir=env)
    expect_s4_class(result, "DFrame")
    expect_identical(colnames(result), c("NumGenes", "NumSig", "PValue", "FDR"))
    expect_identical(rownames(result), names(sets))

    expect_false(anyNA(result$FDR))
    expect_equal(result$NumGenes, c(20, 51, 101))
    expect_equal(result$NumSig, c(0, 51, 0))

    expect_equal(result$PValue[1], 1)
    expect_lt(result$PValue[2], 0.1)
    expect_equal(result$PValue[3], 1)

    # p-values follow the enrichment trend correctly.
    env <- new.env()
    env$FOO <- sets
    env$BAR <- logical(ngenes)
    env$BAR[10:20] <- TRUE
    env$BAR[70:80] <- TRUE
    env$BAR[210:220] <- TRUE

    result <- eval(parse(text=cmds), envir=env)
    expect_equal(result$NumSig, rep(11, nrow(result)))
    expect_lt(result$PValue[1], result$PValue[2])
    expect_lt(result$PValue[2], result$PValue[3])

    # Further testing, this time with sets of the same size but different numbers of genes.
    env <- new.env()
    env$FOO <- list(alpha=1:50, bravo=101:150, charlie=201:250)
    env$BAR <- logical(ngenes)
    env$BAR[c(1:10, 111:130, 201:240)] <- TRUE

    result <- eval(parse(text=cmds), envir=env)
    expect_equal(result$NumGenes, rep(50, nrow(result)))
    expect_equal(result$NumSig, c(10, 20, 40))
    expect_gt(result$PValue[1], result$PValue[2])
    expect_gt(result$PValue[2], result$PValue[3])
})

test_that("generateHypergeometricTestCommands works with other alternative hypotheses", {
    ngenes <- 1000
    sets <- list(
        A=1:20,
        B=451:470,
        C=911:930
    )

    expect_error(augere.gsea:::.generateHypergeometricTestCommands(sets.name="FOO", is.sig.name="BAR", sign.name=NULL, alternative="up"), "'sign'")

    # Set A is only up, set B is only down, set C is 50:50.
    # Only the first half is significant.
    signs <- rep(c(-1, 1), length.out=ngenes)
    signs[sets$A] <- 1
    signs[sets$B] <- -1

    env <- new.env()
    env$FOO <- sets
    env$BAR <- logical(ngenes)
    env$BAR[1:10] <- TRUE
    env$BAR[451:460] <- TRUE
    env$BAR[911:920] <- TRUE
    env$WHEE <- signs

    cmds <- augere.gsea:::.generateHypergeometricTestCommands(sets.name="FOO", is.sig.name="BAR", sign.name="WHEE", alternative="up")
    expect_true(any(grepl("sign <-", cmds)))
    result <- eval(parse(text=cmds), envir=env)
    expect_lt(result$PValue[1], result$PValue[3])
    expect_equal(result$PValue[2], 1)
    expect_lt(result$PValue[3], 0.1)
    expect_equal(result$NumSig, c(10, 0, 5))

    cmds <- augere.gsea:::.generateHypergeometricTestCommands(sets.name="FOO", is.sig.name="BAR", sign.name="WHEE", alternative="down")
    result <- eval(parse(text=cmds), envir=env)
    expect_equal(result$PValue[1], 1)
    expect_lt(result$PValue[2], result$PValue[3])
    expect_lt(result$PValue[3], 0.1)
    expect_equal(result$NumSig, c(0, 10, 5))

    cmds <- augere.gsea:::.generateHypergeometricTestCommands(sets.name="FOO", is.sig.name="BAR", sign.name="WHEE", alternative="either")
    result <- eval(parse(text=cmds), envir=env)
    expect_equal(result$Direction[1], "up")
    expect_equal(result$Direction[2], "down")
    expect_equal(result$NumSig, rep(10, nrow(result)))
    expect_identical(result$PValue[1], result$PValue[2])
    expect_lt(result$PValue[2], result$PValue[3])
})

test_that("generateHypergeometricTestCommands works with empty sets", {
    ngenes <- 1000
    sets <- list(
        A=1:10,
        B=integer(0),
        C=990:1000
    )

    env <- new.env()
    env$FOO <- sets
    env$BAR <- rbinom(ngenes, 1, 0.5) == 1

    cmds <- augere.gsea:::.generateHypergeometricTestCommands(sets.name="FOO", is.sig.name="BAR", sign.name=NULL)
    result <- eval(parse(text=cmds), envir=env)
    expect_identical(which(is.na(result$PValue)), 2L)
})
