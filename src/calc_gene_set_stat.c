/************************************************************************
Name:        calc_gene_set_stat.c
Author:      Weil Lai
Description: contains gene set statistical functions (N-GSEA, NT_k, NE_k)
             used for sigPathway package
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
  - April 2006
    - Added more detailed comments
    - Added options for complete (instead of random) permutation
    - Modified functions to allow the use of other test statistics (e.g.,
      f-test from 1-way ANOVA)
************************************************************************/

#include "sigPathway.h"

/************************************************************************/
/*
  Calculate NES statistics (GSEA using the two sample t-statistic
                            as the association measure and with normalization)
  Note: This function is not as optimized as those seen for
        calc_NTk and calc_NEk.  Expect this function to use more time and 
        resources than calc_NTk and calc_NEk.  It is provided so that
        investigators can reproduce the results in Tian et al.
*/
void calc_GSEA(double *Y, int *n_ps, int *ncol, double *phenotype,
	       int *n_gs, int *nsim,
	       int *nprobesV, int *indexV,
	       int *urp, int *verbose,

	       double *es,
	       double *nes,
	       double *p_null, double *gs_pV, double *gs_qV)
{
  int *index_original, *index, *index_inverted, h, i, j, k;
  int zero = 0, two = 2;
  double *tV, *tmpV = NULL;
  double *pperm;
  double *c1, *c2, *s1, temp, smax, fabs_smax;
  double *es_null, *nes_null;

  tV = (double *) Calloc(*n_ps, double);
  
  t_R(Y, n_ps, ncol, phenotype, &two, &zero, tV, tmpV, tmpV);

  /* calculate es for each gene set */
  index_original = (int *) Calloc(*n_ps, int);
  index = (int *) Calloc(*n_ps, int);
  index_inverted = (int *) Calloc(*n_ps, int);
  
  for( i = 0; i < *n_ps; i++ )
    index_original[i] = i;
  Memcpy(index, index_original, *n_ps);
  revsort(tV, index, *n_ps);

  for( i = 0; i < *n_ps; i++ )
    index_inverted[index[i]] = i;

  c1 = (double *) Calloc(*n_gs, double);
  c2 = (double *) Calloc(*n_gs, double);

  s1 = (double *) Calloc(*n_ps, double);

  h = 0;
  for( i = 0; i < *n_gs; i++ ) {
    c1[i] = -sqrt(((double)nprobesV[i])/(*n_ps-nprobesV[i]));
    c2[i] = sqrt(((double)(*n_ps-nprobesV[i]))/nprobesV[i]);

    for( j = 0; j < *n_ps; j++ )
      s1[j] = c1[i];

    for( j = 0; j < nprobesV[i]; j++ ) {
      s1[index_inverted[indexV[h]]] = c2[i];
      h++;
    }

    temp = s1[0];
    smax = temp;
    fabs_smax = fabs(smax);

    for( j = 1; j < *n_ps; j++ ) {
      temp += s1[j];
      if( fabs(temp) > fabs_smax ) {
	smax = temp;
	fabs_smax = fabs(smax);
      }
    }

    es[i] = smax;
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating es\n");

  /*
    Calculate the null distributions of gene set statistics
      - urp == 0: complete permutation (by lexicographical order)
      - urp != 0: random permutation

    In both cases, we avoid permuted phenotypes that are equal to original
    phenotypes.
  */
  pperm = (double *) Calloc(*ncol, double);
  es_null = (double *) Calloc(*nsim * (*n_gs), double);

  if(*urp == 0) {
    Memcpy(pperm, phenotype, *ncol);
    R_rsort(pperm, *ncol);
  }

  i = 0;
  while(i < *nsim) {
    if(*urp != 0)
      sampleNR_double(phenotype, pperm, ncol);

    if( memcmp(pperm, phenotype, *ncol*sizeof(double)) != 0 ) {

      t_R(Y, n_ps, ncol, pperm, &two, &zero, tV, tmpV, tmpV);
      
      Memcpy(index, index_original, *n_ps);
      revsort(tV, index, *n_ps);
      for( j = 0; j < *n_ps; j++ )
	index_inverted[index[j]] = j;
      
      h = 0;
      for( j = 0; j < *n_gs; j++ ) {
	for( k = 0; k < *n_ps; k++ )
	  s1[k] = c1[j];
	
	for( k = 0; k < nprobesV[j]; k++ ) {
	  s1[index_inverted[indexV[h]]] = c2[j];
	  h++;
	}
	
	temp = s1[0];
	smax = temp;
	fabs_smax = fabs(smax);
	for( k = 1; k < *n_ps; k++ ) {
	  temp += s1[k];
	  if( fabs(temp) > fabs_smax ) {
	    smax = temp;
	    fabs_smax = fabs(smax);
	  }
	}
	
	es_null[i + *nsim*j] = smax;
      }
      
      i++;
    }

    if( *urp == 0 && i < *nsim )
      get_next_perm(ncol, pperm);
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating es_null\n");

  Free(tV);
  Free(index_original);
  Free(index);
  Free(index_inverted);
  Free(c1);
  Free(c2);
  Free(s1);
  Free(pperm);

  calc_internal(verbose, n_gs, nsim, es, es_null, nes);

  nes_null = (double *) Calloc(*nsim * (*n_gs), double);
  calc_internal2(verbose, n_gs, nsim, es_null, nes, nes_null,
		 p_null, gs_pV, gs_qV);
  Free(es_null);
  Free(nes_null);
}
/************************************************************************/


/************************************************************************/
/*
  Calculate NGS_k statistics (permute by row -- probes)

  NOTES:
  This function is a more generalized version of calc_NTk.
*/
void calc_NGSk(double *tV, int *n_ps, int *n_gs, int *nsim,
	       int *nprobesV, int *indexV, int *verbose,

	       double *t_set,
	       double *t_set_new,
	       double *p_null, double *gs_pV, double *gs_qV)
{
  int h, i, j, k;
  double temp, *tpV;
  double *t_set_null, *t_set_new_null;
  
  /* calculate the test statistics for each gene set */
  for(h = 0, i = 0; i < *n_gs; i++) {
    temp = 0.0;
    for(j = 0; j < nprobesV[i]; j++, h++)
      temp += tV[indexV[h]];
    t_set[i] = temp;
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating t_set\n");

  
  /* 
     calculate the null distributions of the test statistics
     - not considering using complete permutation because
       n_ps is usually between 1000 and 100000, which means that n_ps! is
       too large for complete permutation
  */
  tpV = (double *) Calloc(*n_ps, double);
  t_set_null = (double *) Calloc(*nsim * (*n_gs), double);

  i = 0;
  while(i < *nsim) {
    sampleNR_double(tV, tpV, n_ps);
    if( memcmp(tpV, tV, *n_ps*sizeof(double)) != 0 ) {
      for(h = 0, j = 0; j < *n_gs; j++) {
	temp = 0.0;
	for(k = 0; k < nprobesV[j]; k++, h++)
	  temp += tpV[indexV[h]];
	t_set_null[i + *nsim*j] = temp;
      }   

      i++;
    }
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating t_set_null\n");

  Free(tpV);    

  calc_internal(verbose, n_gs, nsim, t_set, t_set_null, t_set_new);

  t_set_new_null = (double *) Calloc(*nsim * (*n_gs), double);
  calc_internal2(verbose, n_gs, nsim, t_set_null, t_set_new, t_set_new_null,
		 p_null, gs_pV, gs_qV);
  Free(t_set_null);
  Free(t_set_new_null);
}
/************************************************************************/


/************************************************************************/
/* 
  Calculate NE_k statistics (permute by column -- samples)

  NOTES: 
         wType = 1: constant weights
         wType = 2: variable weights

         Calculating variable weights does take considerable time and
         memory to complete.

         This function uses a sparse indicator matrix to speed up
         calculations and reduce memory consumption.

	 'urp' stands for "use random permutation"
*/
void calc_NEk(FUNC_STAT stat_fn,
	      double *Y, int *n_ps, int *ncol, double *phenotype,
	      int *n_gs, int *nsim, 
	      int *nprobesV, int *indexV,
	      int *ngroups, int *wType,
	      int *urp, int *verbose,
	      
	      double *t_set,
	      double *t_set_new,
	      double *p_null, double *gs_pV, double *gs_qV)
{
  int h, i, j, k, zero = 0;
  double temp, *pperm, *tV, *tmpV = NULL;
  double *weightsV, *t_set_null, *t_set_new_null;

  /* calculate weights (if specified) */
  if( *wType == 2 ) {
    h = 0;
    for(i = 0; i < *n_gs; i++)
      h += nprobesV[i];
    weightsV = (double *) Calloc(h, double);
    calcWeights(stat_fn, Y, n_ps, ncol, phenotype, n_gs, ngroups, nsim, 
		nprobesV, indexV, &zero, weightsV);
    if( *verbose == 1 )
      Rprintf("Finished calculating variable weights\n");
  }else {
    weightsV = (double *) Calloc(1, double);
  }

  /* calculate probe (set) statistics */
  tV = (double *) Calloc(*n_ps, double);
  (*stat_fn)(Y, n_ps, ncol, phenotype, ngroups, &zero, tV, tmpV, tmpV);

  /* calculate gene set (pathway) statistics */
  if(*wType == 2) {
    for(h = 0, i = 0; i < *n_gs; i++) {
      temp = 0.0;
      for(j = 0; j < nprobesV[i]; j++, h++)
	temp += tV[indexV[h]] * weightsV[h];
      t_set[i] = temp;
    }
  }else {
    for(h = 0, i = 0; i < *n_gs; i++) {
      temp = 0.0;
      for(j = 0; j < nprobesV[i]; j++, h++)
	temp += tV[indexV[h]];
      t_set[i] = temp;
    }
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating t_set\n");


  /*
    Calculate the null distributions of gene set statistics
      - urp == 0: complete permutation (by lexicographical order)
      - urp != 0: random permutation

    In both cases, we avoid permuted phenotypes that are equal to original
    phenotypes.
  */
  pperm = (double *) Calloc(*ncol, double);
  t_set_null = (double *) Calloc(*nsim * (*n_gs), double);

  if(*urp == 0) {
    Memcpy(pperm, phenotype, *ncol);
    R_rsort(pperm, *ncol);
  }

  i = 0;
  while(i < *nsim) {
    if(*urp != 0)
      sampleNR_double(phenotype, pperm, ncol);
   
    if( memcmp(pperm, phenotype, *ncol*sizeof(double)) != 0 ) {

      (*stat_fn)(Y, n_ps, ncol, pperm, ngroups, &zero, tV, tmpV, tmpV);

      if(*wType == 2) {
	for(h = 0, j = 0; j < *n_gs; j++) {
	  temp = 0.0;
	  for(k = 0; k < nprobesV[j]; k++, h++)
	    temp += tV[indexV[h]] * weightsV[h];
	  t_set_null[i + *nsim*j] = temp;
	}
      }else {
	for(h = 0, j = 0; j < *n_gs; j++) {
	  temp = 0.0;
	  for(k = 0; k < nprobesV[j]; k++, h++)
	    temp += tV[indexV[h]];
	  t_set_null[i + *nsim*j] = temp;
	}
      }

      i++;
    }
    
    if( *urp == 0 && i < *nsim )
      get_next_perm(ncol, pperm);
  }

  if( *verbose == 1 )
    Rprintf("Finished calculating t_set_null\n");

  Free(weightsV);
  Free(tV);
  Free(pperm);

  calc_internal(verbose, n_gs, nsim, t_set, t_set_null, t_set_new);

  t_set_new_null = (double *) Calloc(*nsim * (*n_gs), double);
  calc_internal2(verbose, n_gs, nsim, t_set_null, t_set_new, t_set_new_null,
		 p_null, gs_pV, gs_qV);
  Free(t_set_null);
  Free(t_set_new_null);
}
/************************************************************************/
