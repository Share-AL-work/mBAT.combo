#' Define global variables
#'
#' @param df
#' @param columns
#'
#' @return
#' @export


columns_check <- function(df, columns) {
  df_name <- deparse(substitute(df))
  df_columns <- names(df)
  diff_columns <- setdiff(columns, df_columns)
  if (length(diff_columns) > 0) {
    stop("'", df_name, "' doesn't have column(s): ",
      paste(diff_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
}


#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
df_check <- function(df) {
  df_name <- deparse(substitute(df))
  if (!is.data.frame(df)) {
    stop("'", df_name, "' is not a data frame.",
      call. = FALSE
    )
  }
}


#' Title
#'
#' @param x
#' @param arg
#'
#' @return
#' @export
#'
#' @examples
nonneg_num_check <- function(x, arg) {
  if (!is.numeric(x) || is.na(x) || x < 0L || length(x) != 1L) {
    stop("'", arg, "' must be a non-negative number of length 1.",
      call. = FALSE
    )
  }
}

#' Title
#'
#' @param x
#' @param left
#' @param right
#' @param arg
#'
#' @return
#' @export
#'
#' @examples
num_between_check <- function(x, left, right, arg) {
  if (!is.numeric(x) || is.na(x) || x < left || x > right) {
    stop("'", arg, "' must be a number of length 1 between ",
      left, " and ", right, ".",
      call. = FALSE
    )
  }
}

#' Title
#'
#' @param x
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
integer_vector_check <- function(x, tol = .Machine$double.eps) {
  ## this is an empirical solution
  is_wholenumber <- function(x, tol = tol) {
    ## function from is.integer help page
    abs(x - round(x)) < tol && !is.na(x)
  }
  tryCatch(all(sapply(x, is_wholenumber, tol = tol)),
    error = function(e) invisible(FALSE)
  )
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
pretty_num <- function(x, ...) {
  prettyNum(x, big.mark = ",", scientific = FALSE, ...)
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
if_dup_ind <- function(x) {
  which(duplicated(x) | duplicated(x, fromLast = TRUE))
}
