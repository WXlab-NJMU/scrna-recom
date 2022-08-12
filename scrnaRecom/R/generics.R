#' 36L: This is the title.
#'
#' @description
#' This is the description.
#'
#' @details
#' These are further details.
#'
#' @param x desc
#' @param y desc
#' @return Returns ddf
#'
#' @rdname SeuratIntegration
#' @concept integration
#' @export SeuratIntegration
#'
#' @examples
#' add_numbers(1, 2) ## returns 3
#'
#' ## don't run this in calls to 'example(add_numbers)'
#' \dontrun{
#'    SeuratIntegration(2, 3)
#' }
#'
SeuratIntegration <- function(x, ...) {
  UseMethod(generic = 'SeuratIntegration', object = x)
}

