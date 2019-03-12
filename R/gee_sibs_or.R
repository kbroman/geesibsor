# gee.sibs.or
#
#' Estimate odds ratio for siblings adjusting for covariates
#'
#' Use GEE to estimate a common sib-sib odds ratio for a binary phenotype,
#' adjusting for covariates with a logistic model.
#'
#' We assume a set of randomly ascertained sibships with measurements on a
#' binary phenotype and a set of covariates.  We use a logistic model to
#' contact the expected phenotype and the covariates and assume a constant log
#' odds ratio for sibs' phenotypes.  Parameters are estimated using generalized
#' estimating equations (GEE).
#'
#' The idea is described in Liange and Beaty (1991).  All of the details appear
#' in Liang et al. (1992) and Qaqish (1990).
#'
#' Key bits are coded in C; but there's still a lot done directly in R.
#' Eventually I'd like to move more of the code to C, for the sake of speed.
#'
#' @param y A vector of binary phenotypes (coded 0/1).
#' @param x A numeric matrix of covariates, including an intercept term.
#' @param id A vector of integers indicating family/group assignment.
#' @param beta Optional starting values for the covariate coefficients; if
#' `length(beta) = ncol(x)+1`, the last value is assumed to be a starting
#' value for `gamma`.
#' @param gamma Optional starting value for the log odds ratio.  If
#' `gamma=0`, it is assumed to be 0 and is not estimated.
#' @param give.se If true, calculate estimated standard errors.
#' @param return.intercept If TRUE, the estimate for the intercept coefficient
#' is included in the output.
#' @param maxit Maximum number of iterations.  If 0, some debugging information
#' is returned.
#' @param tol Tolerance value for determining convergence in Newton's method to
#' solve the GEE.
#' @param eta.tol Tolerance value for determining convergence in Newton's
#' method to calculate the eta parameters.
#' @param trace Indicates whether to display tracing information; large
#' indicators result in more verbose output.
#' @param debug Indicates whether to display special debugging-related
#' information.
#' @param method Indicates whether to use GEE2, GEE1, or the identity matrix as
#' the working covariance matrix.
#'
#' @return If `give.se=TRUE`, the output is a matrix with four columns:
#' the parameter estimates, estimated SEs, \eqn{(est/SE)^2}, and P-values.  If
#' `give.se=FALSE`, a vector with only the parameter estimates is given.
#'
#' The first item is the sib-sib log odds ratio; the rest are the covariates'
#' coefficients; if `return.intercept=FALSE`, the intercept is excluded
#' from the output.
#'
#' @author Karl W Broman, \email{broman@@wisc.edu}
#'
#' @seealso [fake.data()]
#'
#' @references Liang, K.-Y. and Beaty, T. H. (1991) Measuring familial
#' aggregation by using odds-ratio regression models.  *Genetic
#' Epidemiology*, **8**, 361--370.
#'
#' Liang, K.-Y., Zeger, S. L. and Qaqish, B. (1992) Multivariate regression
#' analyses for categorical data (with discussion).  *J. Roy. Statist.
#' Soc.* B, **1**, 3--40.
#'
#' Qaqish, B. F. (1990) Multivariate regression models using generalized
#' estimating equations.  Ph. D. Thesis, Department of Biostatistics, Johns
#' Hopkins University, Baltimore, Maryland.
#'
#' @keywords models
#'
#' @examples
#' data(fake.data)
#' y <- fake.data[,1]
#' id <- fake.data[,2]
#' x <- cbind(1, fake.data[,-(1:2)])
#' gee.output <- gee.sibs.or(y, x, id, trace=TRUE)
#'
#' @importFrom stats glm binomial pchisq
#' @export gee.sibs.or
gee.sibs.or <-
function(y, x, id, beta=NULL, gamma=0.5, give.se=TRUE, return.intercept=FALSE,
         maxit=1000, tol=1e-5, eta.tol=1e-9, trace=FALSE, debug=FALSE,
         method=c("gee2","gee1","identity"))
{
  # check arguments
  if(!is.null(beta) && length(beta) != ncol(x)) {
    if(length(beta) == ncol(x)+1) {
      gamma <- beta[length(beta)]
      beta <- beta[-length(beta)]
    }
    else stop("length(beta) != ncol(x)")
  }
  if(length(y) != nrow(x)) stop("length(y) != nrow(x)")
  if(length(y) != length(id)) stop("length(y) != length(id)")
  if(any(y!=0 & y !=1)) stop("y should have values 0 and 1")
  method <- match.arg(method)

  if(!maxit) eta.maxit <- 1000
  else eta.maxit <- maxit

  # split id into family categories
  wh <- split(1:length(id), id)
  n.fam <- length(wh)
  n.covar <- ncol(x)

  # starting values
  if(is.null(beta))
    beta <- glm(y ~ -1 + x, family=binomial(link=logit))$coef

  names(beta) <- names(gamma) <- NULL
  if(trace) cat(0,beta,gamma,"\n")

  if(!gamma) est.gamma <- FALSE
  else est.gamma <- TRUE

  flag <- 0
  prev.val <- Inf

  result <- vector("list",n.fam)
  for(i in 1:maxit) {
    for(j in 1:n.fam) {
      result[[j]] <- gee.sibs.or.bits(y[wh[[j]]], x[wh[[j]],,drop=FALSE],
                                      beta, gamma, eta.maxit, eta.tol, trace)

      if(any(is.na(unlist(result[[j]]))))
        warning("gee.sibs.or.bits returned NAs", j, "\n")

      len <- length(wh[[j]])
      if(!est.gamma) { # assume gamma=0 (don't estimate it)
        result[[j]]$matC <- result[[j]]$matC[1:n.covar,1:len,drop=FALSE]
        result[[j]]$matBinv <- result[[j]]$matBinv[1:len,1:len,drop=FALSE]
        result[[j]]$vecA <- result[[j]]$vecA[1:len]
      }

      if(method=="gee2") temp <- result[[j]]$matC %*% solve(result[[j]]$matBinv)
      else if(method=="gee1") {
        result[[j]]$matC[1:n.covar,-(1:len)] <-
          result[[j]]$matC[-(1:n.covar),1:len] <- 0
        result[[j]]$matBinv[1:len,-(1:len)] <-
          result[[j]]$matBinv[-(1:len),1:len] <- 0
        temp <- result[[j]]$matC %*% solve(result[[j]]$matBinv)
      }
      else temp <- result[[j]]$matC # <--- using identity matrix as working covar mat

      if(j==1) {
        first.mat <- temp %*% t(result[[j]]$matC)
        temp2 <- temp %*% result[[j]]$vecA
        second.mat <- temp2
        third.mat <- temp2 %*% t(temp2)
      }
      else {
        first.mat <- first.mat + temp %*% t(result[[j]]$matC)
        temp2 <- temp %*% result[[j]]$vecA
        second.mat <- second.mat + temp2
        third.mat <- third.mat + temp2 %*% t(temp2)
      }
    }
    if(!maxit) return(first.mat,second.mat)

    step <- solve(first.mat, second.mat)

    cur.val <- max(abs(second.mat))
#    if(cur.val > prev.val) {
#      step <- step/32
#      cat("  ** step/32 **\n")
#    }

    if(est.gamma) {
      newbeta <- beta + step[-length(step)]
#      if(cur.val <= prev.val)
        newgamma <- gamma + step[length(step)]
#      else
#        newgamma <- gamma + step[length(step)]
    }
    else {
      newbeta <- beta + step
      newgamma <- gamma
    }

    prev.val <- cur.val

    if(all(abs(beta-newbeta)<tol) & abs(gamma-newgamma)<tol) {
      flag <- 1
      break
    }
    beta <- newbeta; gamma <- newgamma

    if(debug) {
      z <- as.numeric(second.mat)
      u <- rev(order(abs(z)))
      z <- z[u][1:min(c(length(z),3))]
      cat(i,z,gamma,"\n")
    }
    if(trace) cat(i,beta,gamma,"\n")
  }
  if(!flag) warning(" -Didn't converge.\n")

  if(trace>1) {
    cat("\n")
    cat(as.numeric(second.mat),"\n")
    cat("\n")
  }

  output <- c(gamma,beta)

  if(give.se) {
    first.mat <- solve(first.mat)
    var.mat <- first.mat %*% third.mat %*% t(first.mat)
    se <- sqrt(diag(var.mat))
    se <- c(se[length(se)],se[-length(se)]) # reorder so ln(OR) is first

    output <- cbind(output, se, (output/se)^2)
    output <- cbind(output, 1-pchisq(output[,3],1))

    colnames(output) <- c("est","se","Wald","P")
    rownames(output) <- c("lnOR",colnames(x))
    if(!return.intercept) output <- output[-2,,drop=FALSE]
  }
  else {
    names(output) <- c("lnOR",colnames(x))
    if(!return.intercept) output <- output[-2]
  }

  output
}


#' @useDynLib geesibsor, .registration=TRUE

gee.sibs.or.bits <-
function(y, x, beta, gamma, maxit=1000,
         eta.tol=1e-9, trace=FALSE)
{
  n.sibs <- length(y)
  n.covar <- ncol(x)
  if(nrow(x) != n.sibs) stop("nrow(x) != length(y)")
  if(ncol(x) != length(beta)) stop("ncol(x) != length(beta)")

  if(n.sibs > 1) n.vecA <- n.sibs + choose(n.sibs,2)
  else n.vecA <- n.sibs

  result <- .C("R_gee_sibs_or_bits",
               as.integer(n.sibs),
               as.integer(n.covar),
               as.integer(n.vecA),
               as.integer(y),
               as.double(t(x)), # input transpose of covariate matrix
               as.double(beta),
               as.double(gamma),
               matC = as.double(rep(0,(n.covar+1)*n.vecA)),
               matBinv = as.double(rep(0,n.vecA*n.vecA)),
               vecA = as.double(rep(0,n.vecA)),
               as.integer(maxit),
               as.double(eta.tol),
               as.integer(trace),
               PACKAGE="geesibsor")

  return(list(matC=matrix(result$matC, nrow=n.covar+1),
              matBinv=matrix(result$matBinv, nrow=n.vecA),
              vecA=result$vecA))
}
