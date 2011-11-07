/**********************************************************************
 * 
 * gee_sibs_or.h
 *
 * Karl W. Broman
 *
 * first written 17 Sept 2003
 * last modified 23 Sept 2003
 *
 * The goal of this is to estimate the log odds ratio between siblings 
 * for some binary response, after controlling for covariates.  We 
 * seek to solve the GEE corresponding to logit{E(y|x)} = x beta
 * and log OR(y1, y2 | x) = gamma.
 *
 * See the following two papers:
 *
 *   KY Liang and TH Beaty (1991) Measuring familial aggregation by
 *   using odds-ratio regression models.  Genet Epidemiol 8:361-370.
 * 
 *   KY Liang, SL Zeger, B Qaqish (1992) Multivariate regression
 *   analyses for categorical data (with discussion). JRSS B
 *   54(1):3-40.
 *
 * The latter paper has all of the details.
 *
 **********************************************************************/

/**********************************************************************
 * 
 * gee_sibs_or_bits
 * 
 * This function is the main function for the essential bits of the 
 * GEE calculation for a single family.  
 *
 * INPUT:
 *     n_sibs = number of siblings in the family
 *     n_covar = number of covariates
 *     n_vecA = n_sibs + choose(n_sibs, 2) = n_sibs + n_sibs*(n_sibs-1)/2
 *     y = binary phenotype (vector of length n_sibs)
 *     x = covariates (matrix of dimension n_covar x n_sibs; note transpose of usual)
 *     beta = current estimates of coefficients (vector of length n_covar)
 *     gamma = current estimate of log OR for siblings
 *
 * OUTPUT:
 *     MatC = d(mu, eta) / dtheta  [a matrix of size (n_covar+1) x n_vecA]
 *     MatBinv = cov(y, w)         [a matrix of size n_vecA x n_vecA] 
 *     vecA = (y-mu, w-eta)        [a vector of length n_vecA]
 *
 *  Note that MatC and MatBinv are indexed as MatC[col][row]
 **********************************************************************/
void gee_sibs_or_bits(int n_sibs, int n_covar, int n_vecA, 
		      int *y, double **X, double *beta, double gamma,
		      double **MatC, double **MatBinv, double *vecA,
		      int maxit, double eta_tol, int trace);

/* wrapper for call from R */
void R_gee_sibs_or_bits(int *n_sibs, int *n_covar, int *n_vecA, 
			int *y, double *x, double *beta, double *gamma,
			double *matC, double *matBinv, double *vecA,
			int *maxit, double *eta_tol, int *trace);

/**********************************************************************
 * deta12_dmu: calculates derivative of eta12 = E(y1*y2) with respect 
 *             to mu1 = E(y1)
 **********************************************************************/
double deta12_dmu(double mu1, double mu2, double gamma);

/**********************************************************************
 * deta12_dgamma: calculates the derivate of eta12=E(y1*y2) with 
 *                respect to gamma = log OR(y1, y2)
 **********************************************************************/
double deta12_dgamma(double mu1, double mu2, double gamma);


/* end of gee_sibs_or.h */
