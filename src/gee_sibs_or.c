/**********************************************************************
 * 
 * gee_sibs_or.c
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>

#include "gee_sibs_or.h"
#include "calc_eta.h"
#include "util.h"

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
 *  Note that MatC and MatBinv are indexed as MatC[col][row] col=phe, row=covar
 **********************************************************************/
void gee_sibs_or_bits(int n_sibs, int n_covar, int n_vecA, 
		      int *y, double **X, double *beta, double gamma,
		      double **MatC, double **MatBinv, double *vecA,
		      int maxit, double eta_tol, int trace)
{
  int i, j, k, s, u, v, *w, n_pairs;
  int **pair_index, ***trip_index, ****quad_index;
  double *mu, *eta12=0, *eta123=0, *eta1234=0;

  n_pairs = n_vecA - n_sibs;

  /* Allocate space for the mu's, w's and eta's */
  mu = (double *)R_alloc(n_sibs, sizeof(double));
  w = (int *)R_alloc(n_pairs, sizeof(int)); 
  if(n_sibs > 1)
    eta12 = (double *)R_alloc(n_sibs * (n_sibs - 1) / 2, sizeof(double));
  if(n_sibs > 2)
    eta123 = (double *)R_alloc(n_sibs * (n_sibs - 1) * (n_sibs - 2) / (3*2), 
			       sizeof(double));
  if(n_sibs > 3)
    eta1234 = (double *)R_alloc(n_sibs * (n_sibs - 1) * (n_sibs - 2) * (n_sibs - 3) / 
				(4*3*2), sizeof(double));
  
  /* allocate space for index matrices: simplify indexing eta12, eta123, and eta1234 */
  if(n_sibs > 1)
    allocate_imatrix(&pair_index, n_sibs, n_sibs);
  if(n_sibs > 2)
    allocate_3d_iarray(&trip_index, n_sibs, n_sibs, n_sibs);
  if(n_sibs > 3)
    allocate_4d_iarray(&quad_index, n_sibs, n_sibs, n_sibs, n_sibs);

  /* calculate the fitted values */
  for(i=0; i<n_sibs; i++) {
    mu[i] = 0.0;
    for(j=0; j<n_covar; j++) mu[i] += beta[j] * X[i][j];
    mu[i] = exp(mu[i])/(1.0 + exp(mu[i]));
  }
      
  /* calculate the cross-products, eta12's and pair_index (if n_sibs > 1) */
  for(i=0, u=0; i<n_sibs-1; i++) {
    for(j=i+1; j<n_sibs; j++, u++) {
      w[u] = y[i]*y[j];
      eta12[u] = calc_eta12(mu[i], mu[j], gamma, eta_tol);
      pair_index[i][j] = pair_index[j][i] = u;
    }
  }

  /* calculate the eta123's and trip_index (if n_sibs > 2) */
  for(i=0, u=0; i<n_sibs-2; i++) {
    for(j=i+1; j<n_sibs-1; j++) {
      for(k=j+1; k<n_sibs; k++, u++) {
	eta123[u] = calc_eta123(mu[i], mu[j], mu[k],
				eta12[pair_index[i][j]], 
				eta12[pair_index[i][k]],
				eta12[pair_index[j][k]],
				maxit, eta_tol, trace);
	trip_index[i][j][k] = trip_index[i][k][j] = trip_index[j][i][k] = 
	  trip_index[j][k][i] = trip_index[k][i][j] = trip_index[k][j][i] = u;
      }
    }
  }

  /* calculate the eta1234's and quad_index (if n_sibs > 3) */
  for(i=0, u=0; i<n_sibs-3; i++) {
    for(j=i+1; j<n_sibs-2; j++) {
      for(k=j+1; k<n_sibs-1; k++) {
	for(s=k+1; s<n_sibs; s++, u++) {
	  eta1234[u] = calc_eta1234(mu[i], mu[j], mu[k], mu[s],
				    eta12[pair_index[i][j]], 
				    eta12[pair_index[i][k]],
				    eta12[pair_index[i][s]],
				    eta12[pair_index[j][k]], 
				    eta12[pair_index[j][s]],
				    eta12[pair_index[k][s]],
				    eta123[trip_index[i][j][k]],
				    eta123[trip_index[i][j][s]],
				    eta123[trip_index[i][k][s]],
				    eta123[trip_index[j][k][s]],
				    maxit, eta_tol, trace);
	  quad_index[i][j][k][s] = quad_index[i][j][s][k] = 
	    quad_index[i][k][j][s] = quad_index[i][k][s][j] = 
	    quad_index[i][s][j][k] = quad_index[i][s][k][j] = 
	  quad_index[j][i][k][s] = quad_index[j][i][s][k] = 
	    quad_index[j][k][i][s] = quad_index[j][k][s][i] = 
	    quad_index[j][s][i][k] = quad_index[j][s][k][i] = 
	  quad_index[k][j][i][s] = quad_index[k][j][s][i] = 
	    quad_index[k][i][j][s] = quad_index[k][i][s][j] = 
	    quad_index[k][s][j][i] = quad_index[k][s][i][j] = 
	  quad_index[s][j][k][i] = quad_index[s][j][i][k] = 
	    quad_index[s][k][j][i] = quad_index[s][k][i][j] = 
	    quad_index[s][i][j][k] = quad_index[s][i][k][j] = u;
	}
      }
    }
  }

  /* calculate the vector A */
  for(i=0; i<n_sibs; i++) 
    vecA[i] = (double)y[i] - mu[i];
  for(i=0; i<n_pairs; i++) 
    vecA[n_sibs+i] = (double)w[i] - eta12[i];


  /* Calculate the matrix matBinv = cov(y w) */
  for(i=0; i<n_sibs; i++) { /* first cov(y) */
    MatBinv[i][i] = mu[i]*(1.0-mu[i]);
    for(j=i+1; j<n_sibs; j++) 
      MatBinv[i][j] = MatBinv[j][i] = eta12[pair_index[i][j]] - mu[i]*mu[j];
  }

  for(i=0; i<n_sibs-1; i++) { /* now cov(w) */
    for(j=i+1; j<n_sibs; j++) {
      u = n_sibs + pair_index[i][j]; /* position in MatBinv for first pair */

      for(k=0; k<n_sibs-1; k++) {
	for(s=k+1; s<n_sibs; s++) {

	  v = n_sibs + pair_index[k][s]; /* position in MatBinv for second pair */
    

	  if(i==k && j==s) { /* same pair twice: get variance */
	    MatBinv[u][v] = eta12[pair_index[i][j]]*(1.0-eta12[pair_index[i][j]]);
	  }
	  else if(i==k) { /* pairs (i,j) and (i,s) */
	    MatBinv[u][v] = 
	      eta123[trip_index[i][j][s]] - eta12[pair_index[i][j]] * 
	      eta12[pair_index[i][s]];
	  }
	  else if(i==s) { /* pairs (i,j) and (k,i) */
	    MatBinv[u][v] = 
	      eta123[trip_index[i][j][k]] - eta12[pair_index[i][j]] * 
	      eta12[pair_index[k][i]];
	  }
	  else if(j==k) { /* pairs (i,j) and (j,s) */
	    MatBinv[u][v] = 
	      eta123[trip_index[i][j][s]] - eta12[pair_index[i][j]] *
	      eta12[pair_index[j][s]];
	  }
	  else if(j==s) { /* pairs (i,j) and (k,j) */
	    MatBinv[u][v] = 
	      eta123[trip_index[i][k][j]] - eta12[pair_index[i][j]] *
	      eta12[pair_index[k][j]];
	  }
	  else { /* (i,j) and (k,s) */
	    MatBinv[u][v] = 
	      eta1234[quad_index[i][j][k][s]] - eta12[pair_index[i][j]] *
	      eta12[pair_index[k][s]];
	  }
	}
      }
    }
  }

  for(i=0; i<n_sibs; i++) { /* finally, cov(y, w) */
    for(j=0; j<n_sibs-1; j++) {
      for(k=j+1; k<n_sibs; k++) {
	v = n_sibs + pair_index[j][k];
	
	if(i==j) { /* i and (i,k) */
	  MatBinv[i][v] = MatBinv[v][i] = 
	    eta12[pair_index[i][k]]*(1.0-mu[i]);
	}
	else if(i==k) { /* i and (j,i) */
	  MatBinv[i][v] = MatBinv[v][i] = 
	    eta12[pair_index[j][i]]*(1.0-mu[i]);
	}
	else { /* i and (j,k) */
	  MatBinv[i][v] = MatBinv[v][i] =
	    eta123[trip_index[i][j][k]] - mu[i]*eta12[pair_index[j][k]];
	}
      }
    }
  }


  /* The last bit: MatC = d(mu, eta) / dtheta */
  for(i=0; i<n_sibs; i++) { /* first derivative of fitted values rel to betas */
    for(k=0; k<n_covar; k++) /* note that MatC indexed as MatC[phe][covar] */
      MatC[i][k] = X[i][k]*mu[i]*(1.0-mu[i]);

    MatC[i][n_covar] = 0.0;
  }

  /* now derivative of the eta12's */
  for(i=0; i<n_sibs-1; i++) {
    for(j=i+1; j<n_sibs; j++) {
      u = n_sibs + pair_index[i][j]; 

      /*
      temp1 = 1.0 - (mu[i]+mu[j])*(1.0 - exp(gamma));
      temp2 = sqrt(temp1*temp1 - 4.0*(exp(gamma)-1.0)*exp(gamma)*mu[i]*mu[j]); 
      */

      for(k=0; k<n_covar; k++) {

	/*
	MatC[u][k] = (MatC[i][k] + MatC[j][k])/2.0 - 
	  (2.0*temp1*(exp(gamma)-1)*(MatC[i][k] + MatC[j][k]) - 
	   4.0*(exp(gamma)-1.0)*exp(gamma)*(mu[i]*MatC[j][k] + mu[j]*MatC[i][k])) /
	  (4.0 * (exp(gamma) - 1.0) * temp2); 
	*/

	MatC[u][k] = deta12_dmu(mu[i], mu[j], gamma) * MatC[i][k] +
	  deta12_dmu(mu[j], mu[i], gamma) * MatC[j][k];
      }

      /*
      MatC[u][n_covar] = -exp(gamma)/(2.0*(exp(gamma)-1.0)*(exp(gamma)-1.0)) -
	((2.0*(exp(gamma)-1.0)*(2.0*temp1*(mu[i]+mu[j])*exp(gamma)) -
	  4.0*(2.0*exp(2.0*gamma)-exp(gamma))*mu[i]*mu[j])/(2.0*temp2) - 
	 temp2*2.0*exp(gamma)) / (4.0*(exp(gamma)-1.0)*(exp(gamma)-1.0));
      */

      MatC[u][n_covar] = deta12_dgamma(mu[i], mu[j], gamma);
    }
  }
}


/* wrapper for call from R */
void R_gee_sibs_or_bits(int *n_sibs, int *n_covar, int *n_vecA, 
			int *y, double *x, double *beta, double *gamma,
			double *matC, double *matBinv, double *vecA,
			int *maxit, double *eta_tol, int *trace)
{
  double **X, **MatC, **MatBinv;

  /* reorganize x as a double-indexed matrix, X */
  reorg_dmatrix(x, &X, *n_covar, *n_sibs);

  /* reorganize matC as double-indexed matrix MatC */
  reorg_dmatrix(matC, &MatC, *n_covar+1, *n_vecA);

  /* reorganize matBinv as double-indexed matrix MatBinv */
  reorg_dmatrix(matBinv, &MatBinv, *n_vecA, *n_vecA);

  gee_sibs_or_bits(*n_sibs, *n_covar, *n_vecA, y, X, beta, *gamma,
		   MatC, MatBinv, vecA, *maxit, *eta_tol, *trace);
}


/**********************************************************************
 * deta12_dmu: calculates derivative of eta12 = E(y1*y2) with respect 
 *             to mu1 = E(y1)
 **********************************************************************/
double deta12_dmu(double mu1, double mu2, double gamma)
{
  double exp_gamma, qi_prime, q, ai_prime, a;

  exp_gamma = exp(gamma);

  a = 1.0 + (mu1 + mu2)*(exp_gamma - 1.0);
  q = a*a - 4.0 * exp_gamma * (exp_gamma - 1.0) * mu1 * mu2;

  ai_prime = exp_gamma - 1.0;
  qi_prime = 2.0*a*ai_prime - 4.0 * exp_gamma * (exp_gamma - 1.0) * mu2;

  return( (2.0*ai_prime - qi_prime/sqrt(q))/(4.0*(exp_gamma-1.0)) );
}

/**********************************************************************
 * deta12_dgamma: calculates the derivate of eta12=E(y1*y2) with 
 *                respect to gamma = log OR(y1, y2)
 **********************************************************************/
double deta12_dgamma(double mu1, double mu2, double gamma)
{
  double exp_gamma, qij_prime, q, aij_prime, a;

  exp_gamma = exp(gamma);

  a = 1.0 + (mu1 + mu2)*(exp_gamma - 1.0);
  q = a*a - 4.0 * exp_gamma * (exp_gamma - 1.0) * mu1 * mu2;

  aij_prime = (mu1 + mu2);
  qij_prime = 2.0*a*aij_prime - 4.0 * (2.0*exp_gamma - 1.0) * mu1 * mu2;

  return( exp_gamma * ((aij_prime - 0.5*qij_prime/sqrt(q))*(exp_gamma-1.0) - (a - sqrt(q))) /
	  (2.0*(exp_gamma-1.0)*(exp_gamma-1.0)) );
}


/* end of gee_sibs_or.c */
