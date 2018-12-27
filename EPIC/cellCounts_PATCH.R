ns.orig <- get(".harmonizeDataFrames", envir = asNamespace("minfi"))
.harmonizeDataFrames.quickfix <- function (x, y) 
{
  stopifnot(is(x, "DataFrame"))
  stopifnot(is(y, "DataFrame"))
  x.only <- setdiff(names(x), names(y))
  y.only <- setdiff(names(y), names(x))
  if (length(x.only) > 0) {
    df.add <- x[1, x.only]
    is.na(df.add[1, ]) <- TRUE
    y <- cbind(y, df.add)
  }
  if (length(y.only) > 0) {
    df.add <- data.frame(matrix(ncol = length(y.only), nrow = dim(x)[1]))
    names(df.add) <- y.only
    x <- cbind(x, df.add)
  }
  list(x = x, y = y[, names(x)])
}
environment(.harmonizeDataFrames.quickfix) <- environment(ns.orig)
attributes(.harmonizeDataFrames.quickfix) <- attributes(ns.orig)
assignInNamespace(".harmonizeDataFrames", .harmonizeDataFrames.quickfix, ns="minfi")