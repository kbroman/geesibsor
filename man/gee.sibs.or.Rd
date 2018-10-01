\name{gee.sibs.or}
\alias{gee.sibs.or}

\title{Estimate odds ratio for siblings adjusting for covariates}

\description{
  Use GEE to estimate a common sib-sib odds ratio for a binary
  phenotype, adjusting for covariates with a logistic model.
}

\usage{
gee.sibs.or(y, x, id, beta=NULL, gamma=0.5, give.se=TRUE, return.intercept=FALSE,
            maxit=1000, tol=1e-5, eta.tol=1e-9, trace=FALSE, debug=FALSE,
            method=c("gee2","gee1","identity"))}
\arguments{
 \item{y}{A vector of binary phenotypes (coded 0/1).}
 \item{x}{A numeric matrix of covariates, including an intercept term.}
 \item{id}{A vector of integers indicating family/group assignment.}
 \item{beta}{Optional starting values for the covariate coefficients; if
   \code{length(beta) = ncol(x)+1}, the last value is assumed to be a
   starting value for \code{gamma}.}
 \item{gamma}{Optional starting value for the log odds ratio.  If
   \code{gamma=0}, it is assumed to be 0 and is not estimated.}
 \item{give.se}{If true, calculate estimated standard errors.}
 \item{return.intercept}{If TRUE, the estimate for the intercept
   coefficient is included in the output.}
 \item{maxit}{Maximum number of iterations.  If 0, some debugging
   information is returned.}
 \item{tol}{Tolerance value for determining convergence in Newton's
   method to solve the GEE.}
 \item{eta.tol}{Tolerance value for determining convergence in Newton's
   method to calculate the eta parameters.}
 \item{trace}{Indicates whether to display tracing information; large
   indicators result in more verbose output.}
 \item{debug}{Indicates whether to display special debugging-related
   information.}
 \item{method}{Indicates whether to use GEE2, GEE1, or the identity
   matrix as the working covariance matrix.}
}

\value{
  If \code{give.se=TRUE}, the output is a matrix with four columns: the parameter
  estimates, estimated SEs, \eqn{(est/SE)^2}, and P-values.  If
  \code{give.se=FALSE}, a vector with only the parameter estimates is
  given.

  The first item is the sib-sib log odds ratio; the rest are the
  covariates' coefficients; if \code{return.intercept=FALSE}, the
  intercept is excluded from the output.
}

\details{
  We assume a set of randomly ascertained sibships with measurements on
  a binary phenotype and a set of covariates.  We use a logistic model
  to contact the expected phenotype and the covariates and assume a
  constant log odds ratio for sibs' phenotypes.  Parameters are
  estimated using generalized estimating equations (GEE).

  The idea is described in Liange and Beaty (1991).  All of the details
  appear in Liang et al. (1992) and Qaqish (1990).

  Key bits are coded in C; but there's still a lot done directly in R.
  Eventually I'd like to move more of the code to C, for the sake of speed.
}

\references{
  Liang, K.-Y. and Beaty, T. H. (1991) Measuring familial aggregation by
  using odds-ratio regression models.  \emph{Genetic Epidemiology},
  \bold{8}, 361--370.

  Liang, K.-Y., Zeger, S. L. and Qaqish, B. (1992) Multivariate
  regression analyses for categorical data (with discussion).
  \emph{J. Roy. Statist. Soc.} B, \bold{1}, 3--40.

  Qaqish, B. F. (1990) Multivariate regression models using generalized
  estimating equations.  Ph. D. Thesis, Department of Biostatistics,
  Johns Hopkins University, Baltimore, Maryland.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.data)
y <- fake.data[,1]
id <- fake.data[,2]
x <- cbind(1, fake.data[,-(1:2)])
gee.output <- gee.sibs.or(y, x, id, trace=TRUE)
}

\seealso{ \code{\link{fake.data}} }


\keyword{models}
