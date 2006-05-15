/************************************************************************
Name:        calc_weights.c
Author:      Weil Lai
Description: internal functions for calculating variable weights
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
************************************************************************/

#include "sigPathway.h"

/************************************************************************/
/*
  Used in estimating weights (called in calcWeight_common)
*/
void internal_weight(double *e_vectors, double *e_values, int *np,
		     double *qq0, double *lambda, int *only_weights,
		     double *weight, double *minweight)
{
  int i;
  double *tempV, *temp_weight;

  tempV = (double *) Calloc(*np, double);
  temp_weight = (double *) Calloc(*np, double);
  
  for(i = 0; i < *np; i++)
    tempV[i] = qq0[i]/(e_values[i]+*lambda);
  matprod(e_vectors, *np, *np, tempV, *np, 1, temp_weight);

  if(*only_weights == 1) {
    Memcpy(weight, temp_weight, *np);
  }
  else {
    *minweight = temp_weight[0];
    for(i = 1; i < *np; i++) {
      if(temp_weight[i] < *minweight)
	*minweight = temp_weight[i];
    }
  }
    
  Free(tempV);
  Free(temp_weight);
}
/************************************************************************/


/************************************************************************/
/* 
  Used in estimating weights
*/
void calcWeight_common(double *v, int np, int nused, int *verbose, double *w)
{
  int i, j, zero = 0, one = 1, temp_int;
  double *sdV, *rho;
  double *e_values, *e_vectors, *qq0, c0, lambda, lambdaL, lambdaU;
  double minweight, tmp;
  const double tolerance = 0.001;

  sdV = (double *) Calloc(np, double);
  for(i = 0; i < np; i++)
    sdV[i] = sqrt(v[i+np*i]);
  
  rho = (double *) Calloc(np * np, double);
  for(i = 0; i < np; i++) {
    for(j = 0; j <= i; j++) {
      temp_int = i + np*j;
      rho[temp_int] = v[temp_int] / sdV[i] / sdV[j];
    }
  }
  Free(sdV);
  
  e_values = (double *) Calloc(np, double);
  e_vectors = (double *) Calloc(np * np, double);
  eigen(rho, &np, e_values, e_vectors);
  Free(rho);
  
  if(*verbose == 1)
    Rprintf("\tDone with calculating eigenvectors and eigenvalues\n");
    
  /* make sure initial lambda is positive */
  c0 = (e_values[0] > tolerance) ? e_values[0] : tolerance;
  for(i = 1; i < nused; i++) {
    if(e_values[i] < c0 && e_values[i] > tolerance)
      c0 = e_values[i];
  }
  lambda = c0;
  if(*verbose == 1)
    Rprintf("\tInitial lambda = %f\n", lambda);
  
  qq0 = (double *) Calloc(np, double);
  for(j = 0; j < np; j++) {
    for(i = 0; i < np; i++)
      qq0[j] += e_vectors[i + np*j];
  }
  
  internal_weight(e_vectors, e_values, &np, qq0, &lambda,
		  &zero, &tmp, &minweight);
  
  if(*verbose == 1)
    Rprintf("\tInitial minweight = %f\n", minweight);
  
  lambdaL = 0.0;
  while(minweight <= 0) {
    lambdaL = lambda;
    lambda = 2*lambdaL;
    
    if(*verbose == 1) {
      Rprintf("\tlambdaL = %f\n", lambdaL);
      Rprintf("\tlambda = %f\n", lambda);
    }
    
    internal_weight(e_vectors, e_values, &np, qq0, &lambda,
		    &zero, &tmp, &minweight);
    
    if(*verbose == 1)
      Rprintf("\tNew minweight = %f\n", minweight);
  }

  /* find smallest lambda (within tolerance) that gives the smallest 
     non-negative minimum weight */
  lambdaU = lambda;
  if(lambdaU != c0) {
    while(lambdaU - lambdaL > tolerance && minweight <= 0) {
      lambda = 0.5*(lambdaL + lambdaU);
      
      internal_weight(e_vectors, e_values, &np, qq0, &lambda,
		      &zero, &tmp, &minweight);
      
      if(minweight <= 0)
	lambdaL = lambda;
      else
	lambdaU = lambda;
    }      
  }
  
  internal_weight(e_vectors, e_values, &np, qq0, &lambda,
		  &one, w, &tmp);
  
  Free(e_values);
  Free(e_vectors);
  Free(qq0);
}
/************************************************************************/


/************************************************************************/
/*
  Estimate weights for 2 groups
*/
void calcWeights2Groups(double *Y, int *n_ps, int *ncol, double *phenotype,
			int *indexV,
			int *nprobesV, int *n_gs, int *verbose,
			double *weightsV)
{
  int h, i, j, n0, n1, k1, k2, np, subset0, subset1, zero = 0, temp_int, nused;
  double *Y_subset0, *Y_subset1, *v0, *v1, *v;

  n0 = n1 = 0;
  for(i = 0; i < *ncol; i++) {
    if(phenotype[i] == 0.0)
      n0++;
    else
      n1++;
  }
  

  k1 = k2 = 0;

  for(h = 0; h < *n_gs; h++) {
    if(*verbose == 1)
      Rprintf("h = %d\n", h);

    np = nprobesV[h];
    Y_subset0 = (double *) Calloc(n0 * np, double);
    Y_subset1 = (double *) Calloc(n1 * np, double);

    for(i = 0; i < np; i++) {
      subset0 = subset1 = 0;
      for(j = 0; j < *ncol; j++) {
	if(phenotype[j] == 0.0) {
	  Y_subset0[subset0 + n0*i] = Y[indexV[k1] + *n_ps*j];
	  subset0++;
	}else {
	  Y_subset1[subset1 + n1*i] = Y[indexV[k1] + *n_ps*j];
	  subset1++;
	}
      }
      k1++;
    }

    if(*verbose == 1)
      Rprintf("\tDone subsetting Y\n");

    v0 = (double *) Calloc(np * np, double);
    v1 = (double *) Calloc(np * np, double);
    v = (double *) Calloc(np * np, double);

    /* Calculate lower triangle of covariance matrix */
    covar_mat(Y_subset0, &n0, &np, &zero, v0);
    covar_mat(Y_subset1, &n1, &np, &zero, v1);
    for(i = 0; i < np; i++) {
      for(j = 0; j <= i; j++) {
	temp_int = i + np*j;
	v[temp_int] = v0[temp_int]/(n0) + v1[temp_int]/(n1);
      }
    }
    Free(Y_subset0);
    Free(Y_subset1);
    Free(v0);
    Free(v1);

    nused = (np < (*ncol-2)) ? np : (*ncol-2);

    calcWeight_common(v, np, nused, verbose, &weightsV[k2]);
    k2 += np;
    Free(v);
  }

}
/************************************************************************/


/************************************************************************/
/*
  Estimate weights for 1 group
*/
void calcWeights1Group(FUNC_STAT stat_fn,
		       double *Y, int *n_ps, int *ncol, double *phenotype,
		       int *nsim,
		       int *indexV,
		       int *nprobesV, int *n_gs, int *verbose,
		       double *weightsV)
{
  int h, i, j, k1, k2, np, zero = 0, one = 1, nused;
  double *z_null_subset, *z_null, *v;

  z_null = (double *) Calloc(*nsim * (*n_ps), double);
  null_c(stat_fn, &zero, Y, n_ps, ncol, phenotype, &one, nsim, z_null);

  k1 = k2 = 0;
  for(h = 0; h < *n_gs; h++) {
    if(*verbose == 1)
      Rprintf("h = %d\n", h);

    np = nprobesV[h];
    z_null_subset = (double *) Calloc(*nsim * np, double);

    for(j = 0; j < np; j++, k1++) {
      for(i = 0; i < *nsim; i++)
	z_null_subset[i + *nsim * j] = z_null[i + *nsim * indexV[k1]];
    }

    v = (double *) Calloc(np * np, double);
    covar_mat(z_null_subset, nsim, &np, &zero, v);
    Free(z_null_subset);

    nused = (np < *ncol) ? np : *ncol;
    
    calcWeight_common(v, np, nused, verbose, &weightsV[k2]);
    k2 += np;
    Free(v);
  }

  Free(z_null);

}
/************************************************************************/


/************************************************************************/
void calcWeights(FUNC_STAT stat_fn,
		 double *Y, int *n_ps, int *ncol, double *phenotype,
		 int *n_gs, int *ngroups, int *nsim,
		 int *nprobesV, int *indexV, int *verbose, double *weightsV)
{
  if( *ngroups == 2 ) {
    calcWeights2Groups(Y, n_ps, ncol, phenotype, indexV, nprobesV,
		      n_gs, verbose, weightsV);
  }else {
    calcWeights1Group(stat_fn,
		      Y, n_ps, ncol, phenotype, nsim, indexV, nprobesV,
		      n_gs, verbose, weightsV);
  }
}
/************************************************************************/
