# library(testthat); library(augere.gsea); source("test-generateFgseaCommands.R")

set.seed(9999)

test_that("generateFgseaCommands works as expected", {
    ngenes <- 1000
    sets <- list(
        A=1:20,
        B=50:100,
        C=200:300
    )

    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL)
    expect_false(any(grepl("sign <-", cmds)))

    env <- new.env()
    env$FOO <- sets
    env$BAR <- runif(ngenes, -1, 1)
    env$BAR[sets$B] <- rep(c(-10, 10), length.out=length(sets$B))

    suppressWarnings(result <- eval(parse(text=cmds), envir=env))
    expect_s4_class(result, "DFrame")
    expect_identical(colnames(result), c("NumGenes", "ES", "NES", "PValue", "FDR"))
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

    suppressWarnings(result <- eval(parse(text=cmds), envir=env))
    expect_lt(result$PValue[1], result$PValue[2])
    expect_lt(result$PValue[2], result$PValue[3])

    # Further testing, this time with sets of the same size but different numbers of genes.
    env <- new.env()
    env$FOO <- list(alpha=1:50, bravo=101:150, charlie=201:250)
    env$BAR <- runif(ngenes, -1, 1)
    env$BAR[1:10] <- rep(c(-10, 10), length.out=10)
    env$BAR[111:130] <- rep(c(-10, 10), length.out=20)
    env$BAR[211:250] <- rep(c(-10, 10), length.out=40)

    suppressWarnings(result <- eval(parse(text=cmds), envir=env))
    expect_equal(result$NumGenes, rep(50, nrow(result)))
    expect_gt(result$PValue[1], result$PValue[2])
    expect_gt(result$PValue[2], result$PValue[3])
})

test_that("generateFgseaCommands works with other alternative hypotheses", {
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

    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="up")
    suppressWarnings(result <- eval(parse(text=cmds), envir=env))
    expect_lt(result$PValue[1], result$PValue[3])
    expect_lt(result$PValue[3], result$PValue[2])

    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="down")
    suppressWarnings(result <- eval(parse(text=cmds), envir=env))
    expect_lt(result$PValue[2], result$PValue[3])
    expect_lt(result$PValue[3], result$PValue[1])

    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="either", seed=42)
    suppressWarnings(either.result <- eval(parse(text=cmds), envir=env))
    expect_equal(either.result$Direction[1], "up")
    expect_equal(either.result$Direction[2], "down")
    expect_lt(either.result$PValue[1], either.result$PValue[3])
    expect_lt(either.result$PValue[2], either.result$PValue[3])

    # Mixed should give us lower p-values than 'either' for the third set where the signs are mixed.
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="mixed", seed=42)
    suppressWarnings(mixed.result <- eval(parse(text=cmds), envir=env))
    expect_lt(mixed.result$PValue[3], either.result$PValue[3])
})

test_that("generateFgseaCommands works with the sign", {
    ngenes <- 1000
    sets <- list(
        A=1:20,
        B=50:100,
        C=200:300
    )
    stats <- runif(ngenes, -1, 1)

    env <- new.env()
    env$FOO <- sets
    env$BAR <- stats
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="up", seed=69)
    up.result <- eval(parse(text=cmds), envir=env)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="down", seed=69)
    down.result <- eval(parse(text=cmds), envir=env)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL, alternative="either", seed=69)
    either.result <- eval(parse(text=cmds), envir=env)

    # Checking that we get the same result after re-introducing the sign.
    env <- new.env()
    env$FOO <- sets
    env$BAR <- abs(stats)
    env$WHEE <- sign(stats)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name="WHEE", alternative="up", seed=69)
    up.result2 <- eval(parse(text=cmds), envir=env)
    expect_identical(up.result, up.result2)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name="WHEE", alternative="down", seed=69)
    down.result2 <- eval(parse(text=cmds), envir=env)
    expect_identical(down.result, down.result2)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name="WHEE", alternative="either", seed=69)
    either.result2 <- eval(parse(text=cmds), envir=env)
    expect_identical(either.result, either.result2)

    # Checking that we get the same result with a square root.
    env <- new.env()
    env$FOO <- sets
    env$BAR <- stats^2
    env$WHEE <- sign(stats)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name="WHEE", alternative="up", use.sqrt=TRUE, seed=69)
    up.result2 <- eval(parse(text=cmds), envir=env)
    expect_identical(up.result, up.result2)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name="WHEE", alternative="down", use.sqrt=TRUE, seed=69)
    down.result2 <- eval(parse(text=cmds), envir=env)
    expect_identical(down.result, down.result2)
    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name="WHEE", alternative="either", use.sqrt=TRUE, seed=69)
    either.result2 <- eval(parse(text=cmds), envir=env)
    expect_identical(either.result, either.result2)
})

test_that("generateFgseaCommands works with empty sets", {
    ngenes <- 1000
    sets <- list(
        A=1:10,
        B=integer(0),
        C=990:1000
    )

    env <- new.env()
    env$FOO <- sets
    env$BAR <- runif(ngenes, -1, 1)

    cmds <- generateFgseaCommands(sets.name="FOO", stat.name="BAR", sign.name=NULL)
    result <- eval(parse(text=cmds), envir=env)
    expect_identical(which(is.na(result$PValue)), 2L)
})
