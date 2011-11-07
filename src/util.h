/**********************************************************************
 * 
 * util.h
 *
 * Karl W. Broman
 *
 * first written 17 Sept 2003
 * last modified 17 Sept 2003
 *
 * Utility functions for GEE/odds ratio regression
 *
 * See gee_sibs_or.c for the key functions
 *
 **********************************************************************/

/**********************************************************************
 * reorg_dmatrix
 *
 * Reorganize a singly-indexed matrix as a doubly-indexed matrix.
 * 
 * Afterwards, index x by X[col][row]
 **********************************************************************/
void reorg_dmatrix(double *x, double ***X, int n_row, int n_col);

/**********************************************************************
 * reorg_imatrix
 *
 * Reorganize a singly-indexed matrix as a doubly-indexed matrix.
 * 
 * Afterwards, index x by X[col][row]
 **********************************************************************/
void reorg_imatrix(int *x, int ***X, int n_row, int n_col);

/**********************************************************************
 * allocate_imatrix
 *
 * allocate matrix of ints; indexed as imat[row][col]
 **********************************************************************/
void allocate_imatrix(int ***imat, int n_row, int n_col);

/**********************************************************************
 * allocate_3d_iarray
 *
 * allocate 3-dimensional array of ints; indexed as imat[i1][i2][i3]
 **********************************************************************/
void allocate_3d_iarray(int ****imat, int dim1, int dim2, int dim3);

/**********************************************************************
 * allocate_4d_iarray
 *
 * allocate 4-dimensional array of ints; indexed as imat[i1][i2][i3][i4]
 **********************************************************************/
void allocate_4d_iarray(int *****imat, int dim1, int dim2, int dim3, int dim4);

/* end of util.h */

