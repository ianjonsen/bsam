#' Deprecated functions. 
#' 
#' \code{fitSSM}, \code{diagSSM}, and \code{plotSSM}, have been deprecated. Instead use
#' \code{fit_ssm}, \code{diag_ssm} and \code{map_ssm}.
#' @rdname bsam-deprecated
#' @param ... ignored
#' @export
#' @rdname bsam-deprecated 
fitSSM <- function(...) {
  .Deprecated("fit_ssm")
}

#' @export
#' @rdname bsam-deprecated
diagSSM <- function(...) {
  .Deprecated("diag_ssm")
}

#' @export
#' @rdname bsam-deprecated
plotSSM <- function(...) {
  .Deprecated("map_ssm")
}