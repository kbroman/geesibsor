/**********************************************************************
 * 
 * calc_eta.h
 *
 * Karl W. Broman
 *
 * first written 17 Sept 2003
 * last modified 17 Sept 2003
 *
 * These functions calculate E(y1*y2), E(y1*y2*y3) and E(y1*y2*y3*y4)
 *
 * See gee_sibs_or.c for details
 *
 **********************************************************************/

/**********************************************************************
 * calc_eta12:     Calculate eta12 = E(y1*y2 | x)
 *
 * mu1 = E(y1 | x), mu2 = E(y2 | x)
 * gamma = log OR (y1, y2 | x)
 *
 * See equation (9) on pg 9 in Liang et al. (1992).
 **********************************************************************/
double calc_eta12(double mu1, double mu2, double gamma, double tol);

/* wrapper for call from R */
void R_calc_eta12(double *mu1, double *mu2, double *gamma, double *tol, double *result);






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
		   int maxit, double tol, int trace);

/* wrapper for call from R */
void R_calc_eta123(double *mu1, double *mu2, double *mu3, 
		   double *eta12, double *eta13, double *eta23,
		   int *maxit, double *tol, int *trace, double *result);

/* calculates polynomial f(x) for which solution of f(x)=0 gives eta123 */
/* also calculates its derivative f(x) is placed in f[0]; f'(x) in f[1] */
void f_eta123(double x, double mu1, double mu2, double mu3, 
	      double eta12, double eta13, double eta23, double *f);
  





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
		    int maxit, double tol, int trace);

/* wrapper for call from R */
void R_calc_eta1234(double *mu1, double *mu2, double *mu3, double *mu4,
		    double *eta12, double *eta13, double *eta14, 
		    double *eta23, double *eta24, double *eta34,
		    double *eta123, double *eta124, double *eta134, double *eta234,
		    int *maxit, double *tol, int *trace, double *result);

/* calculates polynomial f(x) for which solution of f(x)=0 gives eta1234 */
/* also calculates its derivative f(x) is placed in f[0]; f'(x) in f[1] */
double f_eta1234(double x, double mu1, double mu2, double mu3, double mu4,
		 double eta12, double eta13, double eta14, 
		 double eta23, double eta24, double eta34,
		 double eta123, double eta124, double eta134, double eta234,
		 double *f);

/* end of calc_eta.h */
