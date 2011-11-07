/**********************************************************************
 * 
 * calc_eta.c
 *
 * Karl W. Broman
 *
 * first written 17 Sept 2003
 * last modified 19 Sept 2003
 *
 * These functions calculate E(y1*y2), E(y1*y2*y3) and E(y1*y2*y3*y4)
 *
 * See gee_sibs_or.c for details
 *
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <R_ext/PrtUtil.h>
#include "calc_eta.h"

/**********************************************************************
 * calc_eta12:     Calculate eta12 = E(y1*y2 | x)
 *
 * mu1 = E(y1 | x), mu2 = E(y2 | x)
 * gamma = log OR (y1, y2 | x)
 *
 * See equation (9) on pg 9 in Liang et al. (1992).
 **********************************************************************/
double calc_eta12(double mu1, double mu2, double gamma, double tol) 
{
  double b;

  if(fabs(gamma) < tol) return(mu1*mu2); /* looks independent */

  gamma = exp(gamma); /* convert to odds ratio */

  /* this is simply solving a quadratic equation */
  b = 1.0 - (mu1+mu2)*(1.0-gamma);
  return( (b - sqrt(b*b - 4.0*(gamma-1.0)*gamma*mu1*mu2)) / (2.0 * (gamma-1)) );
}

/* wrapper for call from R */
void R_calc_eta12(double *mu1, double *mu2, double *gamma, double *tol, double *result)
{
  *result = calc_eta12(*mu1, *mu2, *gamma, *tol);
}






/**********************************************************************
 * calc_eta123:     Calculate eta123 = E(y1*y2*y3 | x)
 *
 * mu1 = E(y1 | x), mu2 = E(y2 | x), mu3 = E(y3 | x) 
 * 
 * eta12 = E(y1*y2 | x) (similarly eta13 and eta23) calculated 
 * from mu1, mu2, mu3 and gamma, using calc_eta12().
 *
 * We assume that OR(y1,y2 | x, y3=1) = OR(y1,y2 | x, y3=0).
 *
 * We use Newton's method; maxit = max no. iterations,
 *     tol = tolerance for convergence
 **********************************************************************/
double calc_eta123(double mu1, double mu2, double mu3, 
		   double eta12, double eta13, double eta23,
		   int maxit, double tol, int trace)
{
  int i, flag=0;
  double eta123=0.0, old_eta123, f[2];

  /* starting value = independence */
  old_eta123 = mu1*mu2*mu3;

  flag = 0;
  for(i=0; i<maxit; i++) {
    /* the function (and its derivative) whose zero we seek */
    f_eta123(old_eta123, mu1, mu2, mu3, eta12, eta13, eta23, f);

    /* One step of Newton's method */
    eta123 = old_eta123 - f[0] / f[1];
      
    if(trace>2) 
      Rprintf("%4d %8.6f %8.6f %8.6f\n", i+1, eta123, f[0], f[1]);

    /* check for convergence */ 
    if(fabs(eta123 - old_eta123) < tol && fabs(f[0]) < tol) {
      flag = 1;
      break;
    }

    old_eta123 = eta123;
  }
  if(!flag) warning(" calc_eta123(): didn't converge\n");
  return(eta123);
}

/* wrapper for call from R */
void R_calc_eta123(double *mu1, double *mu2, double *mu3, 
		   double *eta12, double *eta13, double *eta23,
		   int *maxit, double *tol, int *trace, double *result)
{
  *result = calc_eta123(*mu1, *mu2, *mu3, *eta12, *eta13, *eta23,
			*maxit, *tol, *trace);
}

/* calculates polynomial f(x) for which solution of f(x)=0 gives eta123 */
/* also calculates its derivative f(x) is placed in f[0]; f'(x) in f[1] */
void f_eta123(double x, double mu1, double mu2, double mu3, 
	      double eta12, double eta13, double eta23, double *f)
{
  double a[8];
  
  /* the terms in the polynomial */
  a[0] = x;                   /* P(000) */
  a[1] = (mu1-eta12-eta13+x); /* P(100) */
  a[2] = (mu2-eta12-eta23+x); /* P(010) */
  a[3] = (mu3-eta13-eta23+x); /* P(001) */
  a[4] = (eta12-x);           /* P(110) */
  a[5] = (eta13-x);           /* P(101) */
  a[6] = (eta23-x);           /* P(011) */
  a[7] = (1.0-mu1-mu2-mu3+eta12+eta13+eta23-x); /* P(000) */

  /* the function */
  f[0] = a[0]*a[1]*a[2]*a[3] - a[4]*a[5]*a[6]*a[7];

  /* derivative */
  f[1] = a[1]*a[2]*a[3] + a[0]*a[2]*a[3] + a[0]*a[1]*a[3] + a[0]*a[1]*a[2] 
    + a[5]*a[6]*a[7] + a[4]*a[6]*a[7] + a[4]*a[5]*a[7] + a[4]*a[5]*a[6];
}
  





/**********************************************************************
 * calc_eta1234:     Calculate eta1234 = E(y1*y2*y3*y4 | x)
 *
 * mu1 = E(y1 | x), mu2 = E(y2 | x), mu3 = E(y3 | x), mu4 = E(y4 | x)
 * 
 * eta12 = E(y1*y2 | x) (similarly eta13, eta23, eta24, etc) calculated 
 * from mu1, mu2, mu3 and gamma, using calc_eta12().
 *
 * eta123 = E(y1*y2*y3 | x) (similarly for eta124, eta134, eta234)
 * calculated from mu1, mu2, mu3, gamma using calc_eta123().
 *
 * We assume that OR(y1,y2 | x, y3=1) = OR(y1,y2 | x, y3=0) 
 *
 *            and OR(y1,y2 | x, y3=y4=1) + OR(y1,y2 | x, y3=y4=0) =
 *                   OR(y1,y2 | x, y3=1,y4=0) + OR(y1,y2 | x, y3=0,y4=1)
 *
 * We use Newton's method; maxit = max no. iterations,
 *     tol = tolerance for convergence
 **********************************************************************/
double calc_eta1234(double mu1, double mu2, double mu3, double mu4,
		    double eta12, double eta13, double eta14, 
		    double eta23, double eta24, double eta34,
		    double eta123, double eta124, double eta134, double eta234,
		    int maxit, double tol, int trace)
{
  int i, flag=0;
  double eta1234=0.0, old_eta1234, f[2];

  /* starting value = independence */
  old_eta1234 = mu1*mu2*mu3*mu4;

  flag = 0;
  for(i=0; i<maxit; i++) {
    /* the function whose zero we seek */
    f_eta1234(old_eta1234, mu1, mu2, mu3, mu4, eta12, eta13, eta14,
	      eta23, eta24, eta34, eta123, eta124, eta134, eta234, f);

    /* One step of Newton's method */
    eta1234 = old_eta1234 - f[0] / f[1];
      
    if(trace>2) 
      Rprintf("%4d %8.6f %8.6f %8.6f\n", i+1, eta1234, f[0], f[1]);

    /* check for convergence */ 
    if(fabs(eta1234 - old_eta1234) < tol && fabs(f[0]) < tol) {
      flag = 1;
      break;
    }

    old_eta1234 = eta1234;
  }
  if(!flag) warning(" calc_eta1234(): didn't converge\n");
  return(eta1234);
}

/* wrapper for call from R */
void R_calc_eta1234(double *mu1, double *mu2, double *mu3, double *mu4,
		    double *eta12, double *eta13, double *eta14, 
		    double *eta23, double *eta24, double *eta34,
		    double *eta123, double *eta124, double *eta134, double *eta234,
		    int *maxit, double *tol, int *trace, double *result)
{
  *result = calc_eta1234(*mu1, *mu2, *mu3, *mu4, *eta12, *eta13, *eta14, *eta23,
			 *eta24, *eta34, *eta123, *eta124, *eta134, *eta234,
			 *maxit, *tol, *trace);
}

/* calculates polynomial f(x) for which solution of f(x)=0 gives eta1234 */
/* also calculates its derivative f(x) is placed in f[0]; f'(x) in f[1] */
double f_eta1234(double x, double mu1, double mu2, double mu3, double mu4,
		 double eta12, double eta13, double eta14, 
		 double eta23, double eta24, double eta34,
		 double eta123, double eta124, double eta134, double eta234,
		 double *f)
{
  double a[16], temp;
  int i, j;

  a[0] = x; /* P(0000) */

  a[1] = eta12 - eta123 - eta124 + x;  /* P(1100) */
  a[2] = eta13 - eta123 - eta134 + x;  /* P(1010) */
  a[3] = eta14 - eta124 - eta134 + x;  /* P(1001) */
  a[4] = eta23 - eta123 - eta234 + x;  /* P(0110) */
  a[5] = eta24 - eta124 - eta234 + x;  /* P(0101) */
  a[6] = eta34 - eta134 - eta234 + x;  /* P(0011) */

  a[7] = 1.0 - mu1 - mu2 - mu3 - mu4 + eta12 + eta13 + eta14 +
    eta23 + eta24 + eta34 - eta123 - eta124 - eta134 - eta234 + x; /* P(0000) */

  a[8]  = eta123 - x; /* P(1110) */
  a[9]  = eta124 - x; /* P(1101) */
  a[10] = eta134 - x; /* P(1011) */
  a[11] = eta234 - x; /* P(0111) */

  a[12] = mu1 - eta12 - eta13 - eta14 + eta123 + eta124 + eta134 - x; /* P(1000) */
  a[13] = mu2 - eta12 - eta23 - eta24 + eta123 + eta124 + eta234 - x; /* P(0100) */
  a[14] = mu3 - eta13 - eta23 - eta34 + eta123 + eta134 + eta234 - x; /* P(0010) */
  a[15] = mu4 - eta14 - eta24 - eta34 + eta124 + eta134 + eta234 - x; /* P(0001) */

  f[0] = a[0]*a[1]*a[2]*a[3]*a[4]*a[5]*a[6]*a[7] - 
    a[8]*a[9]*a[10]*a[11]*a[12]*a[13]*a[14]*a[15];

  f[1] = 0.0;
  /* Derivative: want sum of products like a[1]*...*a[7] + a[0]*a[2]*...*a[7] + ... */
  for(i=0; i<8; i++) {
    temp=1; 
    for(j=0; j<8; j++) {
      if(i != j) temp *= a[j];
    }
    f[1] += temp;
  }

  for(i=8; i<16; i++) {
    temp=1; 
    for(j=8; j<16; j++) {
      if(i != j) temp *= a[j];
    }
    f[1] += temp;
  }

  return(-1.0);
}

/* end of calc_eta.c */
