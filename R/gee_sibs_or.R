######################################################################
# gee_sibs_or.R
#
# Karl W. Broman
#
# first written 17 Sept 2003
# last modified 16 Oct 2006
#
# The goal of this is to estimate the log odds ratio between siblings
# for some binary response, after controlling for covariates.  We
# seek to solve the GEE corresponding to logit{E(y|x)} = x beta
# and log OR(y1, y2 | x) = gamma.
#
# See the following two papers:
#
#   KY Liang and TH Beaty (1991) Measuring familial aggregation by
#   using odds-ratio regression models.  Genet Epidemiol 8:361-370.
#
#   KY Liang, SL Zeger, B Qaqish (1992) Multivariate regression
#   analyses for categorical data (with discussion). JRSS B
#   54(1):3-40.
#
# The latter paper has all of the details.
#
######################################################################

######################################################################
# gee.sibs.or
#
# id should be a categorical variable indicating the groups to which
# the individuals belong
#
# y should be a 0/1 response
#
# x should be a matrix of numeric covariates (including the intercept)
#
# We then use generalized estimating equations (GEE) to fit a
# model with logit{E(y|x)} = x beta and ln OR(y1,y2|x) = gamma
# for all individuals within a group. (We're thinking of siblings.)
#
# INPUT:
#   y = binary outcome
#   x = matrix of covariates (including intercept)
#   id = vector indicating family/group status
#
#   beta = starting values for beta's (and, in length = ncol(x)+1, for gamma)
#   gamma = starting value for gamma (if =0, gamma is assumed to be
#           0 and is not estimated)
#
#   give.se = if TRUE, return estimated SEs (as an attribute)
#
#   maxit = maximum number of iterations for all iterative methods
#   tol = tolerance for convergence in GEE
#   eta.tol = tolerance for convergence in calculating the etas
#   trace = indicates whether to print tracing info (larger values gives
#           more verbose information)
#   debug = if TRUE, print even more stuff useful for my own debugging
#           purposes.
#
#   method indicates whether to use the full GEE, GEE1 (independence of
#          mu's and eta's) or the identity matrix as the working
#          covariance matrix
#
######################################################################
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

#  if(!is.loaded("R_gee_sibs_or_bits")) {
#    dyn.load("gee_orreg.so")
#    cat(" -Loaded gee_orreg.so\n")
#  }

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



# end of gee_sibs_or.R
