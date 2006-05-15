/************************************************************************
Name:        calc_gene_set_internal.c
Author:      Weil Lai
Description: contains internal functions to calculate normalized gene
             set statistics, p-values, and q-values
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
  - April 2006
    - Added more detailed comments
************************************************************************/

#include "sigPathway.h"

/************************************************************************/
/*
  Normalize the gene set (Tk, Ek, modified GSEA) statistics
  with permuted gene set statistics
  (1) For a gene set k, count the # t_set_null[,k] <= t_set[k]
  (2) If the count is > 0 and < nsim, then transform the gene set statistic
      with qnorm()
  (3) If the count is == 0 or == nsim, then transform the gene set statistic
      to a z-score
*/
void calc_internal(int *verbose, int *n_gs, int *nsim,
		   double *t_set, double *t_set_null, double *t_set_new)
{
  int i, j;
  double temp, mean, var;

  /* calculate the standardized test statistics */
  for( i = 0; i < *n_gs; i++ ) {
    temp = 0.0;
    for( j = 0; j < *nsim; j++ ) {
      if( t_set_null[*nsim*i + j] <= t_set[i] )
	temp += 1.0;
    }

    if( temp > 0 && temp < *nsim ) {
      t_set_new[i] = qnorm(temp / *nsim, 0.0, 1.0, 1, 0);
    }else {
      mean = 0.0;
      for( j = 0; j < *nsim; j++ )
	mean += t_set_null[*nsim*i + j];
      mean /= *nsim;
      
      var = 0.0;
      for( j = 0; j < *nsim; j++ ) {
	var += (t_set_null[*nsim*i + j] - mean) * 
	  (t_set_null[*nsim*i + j] - mean);
      }
      var /= (*nsim - 1);

      t_set_new[i] = (t_set[i] - mean) / sqrt(var);
    }
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating t_set_new\n");
}
/************************************************************************/


/************************************************************************/
/*
  Calculate the normalized null distribution, p-values, and q-values of
  each gene set
  (1) Normalize the null distribution of each gene set by ranking the values
      within each gene set's null distribution and then transforming those
      values with qnorm()
  (2) Calculate p-values with pchisq()
  (3) Estimate p_null (a.k.a. pi_0 in other literature) by "the left
      derivative of the smallest concave function greater than the
      empirical distribution function of" the p-values (Tian et al. 2005)
  (4) Calculate q-values with calcQVFast(), which uses computational
      shortcuts to implement the last equation in the appendix of Tian et al.
*/
void calc_internal2(int *verbose, int *n_gs, int *nsim,
		  double *t_set_null, double *t_set_new,
		  double *t_set_new_null, double *p_null,
		  double *pV, double *qV)
{
  int i, j, xylen;
  double *inV, *rankQV, *x, *y, *resV;

  /* calculate the null distributions of the standardized test statistics */
  inV = (double *) Calloc(*nsim, double);
  rankQV = (double *) Calloc(*nsim, double);

  for( j = 0; j < *n_gs; j++ ) {
    for( i = 0; i < *nsim; i++ )
      inV[i] = -t_set_null[i + *nsim*j];

    rank_avg2(nsim, inV, rankQV);
    for( i = 0; i < *nsim; i++ )
      t_set_new_null[i + *nsim*j] = qnorm(rankQV[i]/(*nsim+1), 0.0, 1.0, 1, 0);
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating t_set_new_null\n");

  /* calculate the proportion of nulls */
  xylen = *n_gs+1;
  x = (double *) Calloc(xylen, double);
  y = (double *) Calloc(xylen, double);
  for( i = 0; i < *n_gs; i++ ) {
    pV[i] = pchisq(t_set_new[i]*t_set_new[i], 1, 0, 0);
    x[i+1] = pV[i];
    y[i+1] = (i+1.0)/(*n_gs);
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating pV\n");

  resV = (double *) Calloc(*n_gs, double);
  R_rsort(x, xylen);
  maj(&xylen, x, y, resV);
  
  p_null[0] = 1.0;
  for( i = 0; i < *n_gs; i++ ) {
    if( resV[i] < p_null[0] && x[i+1] < 0.95 )
      p_null[0] = resV[i];
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating p_null\n");

  /* calculate the q-value for each gene set */

  /*
    double nrej, temp_sum, qtemp;
    for( i = 0; i < *n_gs; i++ ) {
      nrej = 0.0;
      for( j = 0; j < *n_gs; j++ ) {
        if( fabs(t_set_new[j]) >= fabs(t_set_new[i]) )
	  nrej++;
      }
      
      temp_sum = 0.0;
      for( j = 0; j < (*nsim)*(*n_gs); j++ ) {
        if( fabs(t_set_new_null[j]) >= fabs(t_set_new[i]) )
	  temp_sum++;
      }
  
      qtemp = p_null[0] * temp_sum / (*nsim) / nrej;
      if( qtemp > 1.0 )
        qtemp = 1.0;
  
      qV[i] = qtemp;
    }
  */

  calcQVFast(t_set_new, n_gs, p_null, t_set_new_null, nsim, qV);
    
  if( *verbose == 1 )
    Rprintf("Finished calculating qV\n");

  Free(inV);
  Free(rankQV);
  Free(x);
  Free(y);
  Free(resV);

  if( *verbose == 1 )
    Rprintf("Finished freeing temporary vectors\n");
}
/************************************************************************/


/************************************************************************/
/*
  Rank an input vector (based off do_rank() in sort.c of R 2.2.0)
*/
void rank_avg2(int *m, double *inV, double *rank)
{
  int i, j, k, *indexV;
  double *tempV;

  tempV = (double *) Calloc(*m, double);
  Memcpy(tempV, inV, *m);

  indexV = (int *) Calloc(*m, int);
  for(i = 0; i < *m; i++)
    indexV[i] = i;

  rsort_with_index(tempV, indexV, *m);
  i = 0;

  while(i < *m) {
    j = i;

    while( (j < *m-1) && (inV[indexV[j]] == inV[indexV[j+1]]) )
      j++;
    
    if(i != j) {
      for(k = i; k <= j; k++)
	rank[indexV[k]] = (i + j + 2) / 2.0;
    }else {
      rank[indexV[i]] = i + 1.0;
    }

    i = j + 1;
  }
  
  Free(indexV);
  Free(tempV);
}
/************************************************************************/


/************************************************************************/
/*
  Used in estimating p_null
*/
void maj(int *n, double *x, double *y, double *resV)
{
  int i,j;
  int start, end;
  int *group;

  double temp, maxslope;

  start = 0;
  group = (int *) Calloc(*n, int);

  while( start < *n-1 ) {

    j = 0;
    for( i = start+1; i < *n; i++ ) {
      if(x[i] > x[start]) {
	group[j] = i;
	j++;
      }
    }

    maxslope = (y[group[0]] - y[start]) / (x[group[0]] - x[start]);
    end = group[0];
    for( i = 1; i < j; i++ ) {
      temp = (y[group[i]] - y[start]) / (x[group[i]] - x[start]);
      if( temp >= maxslope ) {
	maxslope = temp;
	end = group[i];
      }
    }

    for( i = start; i < end; i++ )
      resV[i] = maxslope;

    start = end;
  }

  Free(group);
}
/************************************************************************/


/************************************************************************/
/*
  Calculate q-values

  This function estimates q-values by comparing the number of rejected
  hypotheses in the permuted test statistic to the number of rejected 
  hypotheses in the observed test statistic.  The appendix section of
  Tian et al. PNAS paper contains more details about this estimation.

  We sort the absolute values of the normalized test statistic vector (tsn)
  and the normalized null test statistic vector (tsn_null) to reduce 
  the number of redundant comparisons.
*/
void calcQVFast(double *tsn, int *ngs,
		double *p_null,
		double *tsn_null, int *nsim,
		double *qV)
{
  int i, j, i_original, len_tsn_null;
  double qtemp, *tsn_abs_sorted, *tsn_null_abs_sorted;
  int *nrejV, *idx_tsn_abs_sorted;

  len_tsn_null = *nsim * (*ngs);

  nrejV = (int *) Calloc(*ngs, int);
  tsn_abs_sorted = (double *) Calloc(*ngs, double);
  idx_tsn_abs_sorted = (int *) Calloc(*ngs, int);

  for(i = 0; i < *ngs; i++) {
    tsn_abs_sorted[i] = fabs(tsn[i]);
    idx_tsn_abs_sorted[i] = i;
  }

  rsort_with_index(tsn_abs_sorted, idx_tsn_abs_sorted, *ngs);

  i = 0;
  nrejV[idx_tsn_abs_sorted[*ngs-1]] = 1;
  while(i < *ngs-1) {
    if(tsn_abs_sorted[i+1] == tsn_abs_sorted[i]) {
      i_original = i;
      
      while(i < *ngs-1) {
	if(tsn_abs_sorted[i+1] != tsn_abs_sorted[i])
	  break;
	i++;
      }

      for(j = i_original; j <= i; j++)
	nrejV[idx_tsn_abs_sorted[j]] = *ngs - i_original;
    }else {
      nrejV[idx_tsn_abs_sorted[i]] = *ngs - i;
    }
    
    i++;
  }

  tsn_null_abs_sorted = (double *) Calloc(len_tsn_null, double);
  for(i = 0; i < len_tsn_null; i++)
    tsn_null_abs_sorted[i] = fabs(tsn_null[i]);

  R_qsort(tsn_null_abs_sorted, 1, len_tsn_null);

  i = 0;
  j = 0;
  while((i < *ngs) & (j < len_tsn_null)) {
    if(tsn_null_abs_sorted[j] >= tsn_abs_sorted[i]) {
      qtemp = *p_null * (len_tsn_null-j) / (*nsim) /
	nrejV[idx_tsn_abs_sorted[i]];
      if(qtemp > 1.0)
	qtemp = 1.0;
      qV[idx_tsn_abs_sorted[i]] = qtemp;

      i++;
    }else
      j++;
  }

  Free(nrejV);
  Free(tsn_abs_sorted);
  Free(tsn_null_abs_sorted);
  Free(idx_tsn_abs_sorted);
}
/************************************************************************/
