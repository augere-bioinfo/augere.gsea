.format_more_args <- function(indent, comma.last=TRUE) {
    deparsed.args <- character()
    for (n in names(args)) {
        deparsed.args <- c(deparsed.args, sprintf("%s=%s", n, deparseToString(args[[n]])))
    }

    indented <- sprintf("%s%s", strrep(" ", indent), deparsed.args)
    if (comma.last) {
        indented <- sprintf("%s,", indented)
    } else {
        indented[-1] <- sprintf("%s,", indented[-1])
    }
    indented
}
