/************************************************************************
Name:        matrix_fn.c
Author:      Weil Lai
Description: internal functions that calculate eigenvectors,
             eigenvalues, matrix cross product, and 
             removal of columns with zeros in sparse indicator matrix
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
************************************************************************/

#include "sigPathway.h"

/************************************************************************/
/*
  Modified version of function from array.c of R 2.2.0
*/
void matprod(double *x, int nrx, int ncx,
	     double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    int i;
    double one = 1.0, zero = 0.0;

    /* assume we are not dealing with matrices containing NA/NaNs */
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
      F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
		      x, &nrx, y, &nry, &zero, z, &nrx);
    }else /* zero-extent operations should return zeroes */
      for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}
/************************************************************************/


/************************************************************************/
/* 
  Based off a special case of modLa_rs in Lapack.c of R 2.2.0
*/
void eigen(double *xin, int *nn, double *values, double *vectors)
{
  int n, i, j, lwork, info = 0;
  char jobv[1], uplo[1], range[1];
  double *work, *rx, *rvalues, *z, tmp;

  int liwork, *iwork, itmp, m;
  double vl = 0.0, vu = 0.0, abstol = 0.0; 
  int il, iu, *isuppz;

  n = *nn;
  jobv[0] = 'V';
  uplo[0] = 'L';
  range[0] = 'A';

  rx = (double *) Calloc(n*n, double);
  Memcpy(rx, xin, n*n);

  rvalues = (double *) Calloc(n, double);
  z = (double *) Calloc(n*n, double);

  isuppz = (int *) Calloc(2*n, int);

  /* ask for optimal size of work arrays */
  lwork = -1; liwork = -1;
  F77_CALL(dsyevr)(jobv, range, uplo, &n, rx, &n,
		   &vl, &vu, &il, &iu, &abstol, &m, rvalues,
		   z, &n, isuppz,
		   &tmp, &lwork, &itmp, &liwork, &info);
  if (info != 0)
    Rprintf("Error code %d from Lapack routine 'dsyevr'", info);

  lwork = (int) tmp;
  liwork = itmp;
  
  work = (double *) Calloc(lwork, double);
  iwork = (int *) Calloc(liwork, int);
  F77_CALL(dsyevr)(jobv, range, uplo, &n, rx, &n,
		   &vl, &vu, &il, &iu, &abstol, &m, rvalues,
		   z, &n, isuppz,
		   work, &lwork, iwork, &liwork, &info);
  if (info != 0)
    Rprintf("Error code %d from Lapack routine 'dsyevr'", info);

  for(i = 0; i < n; i++)
    values[i] = rvalues[n-1-i];

  for(j = 0; j < n; j++)
    Memcpy(&vectors[n*(n-1-j)], &z[n*j], n);

  Free(rx);
  Free(rvalues);
  Free(z);
  Free(isuppz);
  Free(work);
  Free(iwork);
}
/************************************************************************/


/************************************************************************/
/*
  Calculate convariance of two vectors
*/
double covar(double *x, double *y, int n)
{
  int i;
  double m_x = 0.0, m_y = 0.0, ans = 0.0;

  for(i = 0; i < n; i++) {
    m_x += x[i];
    m_y += y[i];
  }
  m_x /= n;
  m_y /= n;

  for(i = 0; i < n; i++)
    ans += (x[i] - m_x)*(y[i] - m_y);

  return(ans/(n-1));
}
/************************************************************************/


/************************************************************************/
/*
  Calculate covariance matrix
*/
void covar_mat(double *x, int *nr, int *nc, int *fill, double *output)
{
  int i, j;

  /* calculate covariance and store results to lower triangle of output */
  for(i = 0; i < *nc; i++) {
    for(j = 0; j <= i; j++)
      output[i + *nc*j] = covar(&x[*nr*i], &x[*nr*j], *nr);
  }

  if( *fill == 1 ) {
    /* populate upper triangle of output with values from lower triangle */
    for(i = 0; i < *nc-1; i++) {
      for(j = i+1; j < *nc; j++)
	output[i + *nc*j] = output[j + *nc*i];
    }
  }
}
/************************************************************************/


/************************************************************************/
/*
  Remove probe sets not represented by any of the gene sets
  This pre-processing function speeds up the permutation steps in calc_NEk
*/
void remove_zero_cols(int *nprobesV, int *indexV, int *n_gs, int *n_ps,
		      int *indexV_new, int *ps_in_gs)
{
  int h, i, np;
  int *col_sumsV, *remap_indexV;

  /* calculate column sums */
  np = 0;
  for(i = 0; i < *n_gs; i++)
    np += nprobesV[i];

  col_sumsV = (int *) Calloc(*n_ps, int);
  for(i = 0; i < np; i++)
    col_sumsV[indexV[i]]++;

  /* identify which columns have sums > 0 and make remapping vector */
  remap_indexV = (int *) Calloc(*n_ps, int);
  h = 0;
  for(i = 0; i < *n_ps; i++) {
    if(col_sumsV[i] > 0) {
      ps_in_gs[i] = 1;
      remap_indexV[i] = h;
      h++;
    }
    else
      ps_in_gs[i] = 0;
  }

  /* regenerate new indices (because probe sets not associated with 
     gene sets of interest will be removed) */
  for(i = 0; i < np; i++)
    indexV_new[i] = remap_indexV[indexV[i]];

  Free(col_sumsV);
  Free(remap_indexV);
}
/************************************************************************/
