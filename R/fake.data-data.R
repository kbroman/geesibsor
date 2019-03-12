#' Fake data
#'
#' Fake data for illustrating and testing the R/geesibsor package.
#'
#' There are data on 696 individuals in 270 sibships.  The first column is the
#' phenotype.  The next is an indicator of family assignment.  The last 6
#' columns are covariates.
#'
#' @format A numeric matrix with 696 rows and 8 columns.
#'
#' @seealso \code{\link{gee.sibs.or}}
#'
#' @source These are completely fake data.  Email Karl W. Broman ,
#' \email{broman@wisc.edu}, for further information.
#'
#' @keywords datasets
#'
#' @examples
#' data(fake.data)
#' y <- fake.data[,1]
#' id <- fake.data[,2]
#' x <- cbind(1, fake.data[,-(1:2)])
#' gee.output <- gee.sibs.or(y, x, id)
#'
"fake.data"
