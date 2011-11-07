/**********************************************************************
 * 
 * util.c
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

#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include "util.h"

/**********************************************************************
 * reorg_dmatrix
 *
 * Reorganize a singly-indexed matrix as a doubly-indexed matrix.
 * 
 * Afterwards, index x by X[col][row]
 **********************************************************************/
void reorg_dmatrix(double *x, double ***X, int n_row, int n_col)
{
  int i;
  
  *X = (double **)R_alloc(n_col, sizeof(double *));
  
  (*X)[0] = x;
  for(i=1; i<n_col; i++) 
    (*X)[i] = (*X)[i-1] + n_row;
}

/**********************************************************************
 * reorg_imatrix
 *
 * Reorganize a singly-indexed matrix as a doubly-indexed matrix.
 * 
 * Afterwards, index x by X[col][row]
 **********************************************************************/
void reorg_imatrix(int *x, int ***X, int n_row, int n_col)
{
  int i;
  
  *X = (int **)R_alloc(n_col, sizeof(int *));
  
  (*X)[0] = x;
  for(i=1; i<n_col; i++) 
    (*X)[i] = (*X)[i-1] + n_row;
}

/**********************************************************************
 * allocate_imatrix
 *
 * allocate matrix of ints; indexed as imat[row][col]
 **********************************************************************/
void allocate_imatrix(int ***imat, int n_row, int n_col)
{
  int *a;
  
  /* contingous space for the thing */
  a = (int *)R_alloc(n_row * n_col, sizeof(int));
  
  reorg_imatrix(a, imat, n_col, n_row);
}

/**********************************************************************
 * allocate_3d_iarray
 *
 * allocate 3-dimensional array of ints; indexed as imat[i1][i2][i3]
 **********************************************************************/
void allocate_3d_iarray(int ****imat, int dim1, int dim2, int dim3)
{
  int *a, **b;
  int i, j;
  
  /* contingous space for the thing */
  a = (int *)R_alloc(dim1 * dim2 * dim3, sizeof(int));
  
  *imat = (int ***)R_alloc(dim1, sizeof(int **));
  b = (int **)R_alloc(dim1*dim2, sizeof(int *));

  (*imat)[0] = b;
  for(i=1; i<dim1; i++) 
    (*imat)[i] = (*imat)[i-1] + dim2;
  
  for(i=0; i<dim1; i++) 
    for(j=0; j<dim2; j++)
      (*imat)[i][j] = a + i*dim2*dim3 + j*dim3;
}

/**********************************************************************
 * allocate_4d_iarray
 *
 * allocate 4-dimensional array of ints; indexed as imat[i1][i2][i3][i4]
 **********************************************************************/
void allocate_4d_iarray(int *****imat, int dim1, int dim2, int dim3, int dim4)
{
  int *a, **b, ***c;
  int i, j, k;
  
  /* contingous space for the thing */
  a = (int *)R_alloc(dim1 * dim2 * dim3 * dim4, sizeof(int));
  
  *imat = (int ****)R_alloc(dim1, sizeof(int ***));
  c = (int ***)R_alloc(dim1*dim2, sizeof(int **));
  b = (int **)R_alloc(dim1*dim2*dim3, sizeof(int *));

  (*imat)[0] = c;
  for(i=1; i<dim1; i++) 
    (*imat)[i] = (*imat)[i-1] + dim2;
  
  for(i=0; i<dim1; i++) 
    for(j=0; j<dim2; j++)
      (*imat)[i][j] = b + (i*dim2 + j)*dim3;

  for(i=0; i<dim1; i++) 
    for(j=0; j<dim2; j++)
      for(k=0; k<dim3; k++)
	(*imat)[i][j][k] = a + ((i*dim2+j)*dim3+k)*dim4;
}

/* end of util.c */

