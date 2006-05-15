/************************************************************************
Name:        sigPathway.c
Author:      Weil Lai
Description: wrapper functions to run pathway analysis
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
  - April 2006
    - Added option for the algorithms to use complete permutation
      when nsim is determined (elsewhere) to be greater than
      the total number of unique permutations
************************************************************************/

#include "sigPathway.h"

/************************************************************************/
/*
  Main wrapper function for the pathway analysis calculations

  I/O    Variable   R class   Length      Description
  ---------------------------------------------------------------------------
  INPUT  Y          numeric   n_ps*ncol   Expression matrix in vector,
                                          column-major format
  INPUT  n_ps       integer   1           # probe sets (also # rows in Y)
  INPUT  ncol       integer   1           # samples (also # cols in Y)
  INPUT  phenotype  numeric   ncol        Phenotype vector
  INPUT  n_gs       integer   1           # gene sets
  INPUT  nsim       integer   1           # max permutations
  INPUT  nprobesV   integer   n_gs        Sparse gene set indicator matrix
                                          (# probe sets in each gene set)
  INPUT  indexV     integer   n_gs        Sparse gene set indicator matrix
                                          (the column indices of Y for each
                                           "1" in the sparse matrix)
  INPUT  ngroups    integer   1           # groups
  INPUT  testType   character (varies)    Either "GSEA", "NGSk", "NTk", "NEk"
  INPUT  weightType character (varies)    Either "constant" or "variable"
  INPUT  urp        integer   1           != 0: always use random permutation
                                          == 0: if the number of unique
                                                phenotype permutations is less
                                                than nsim, use complete perm
  INPUT  verbose    integer   1           == 1: print messages to screen
                                          == 0: no messages
  OUTPUT t_set      numeric   n_gs        gene set statistics
  OUTPUT t_set_new  numeric   n_gs        normalized gene set statistics
  OUTPUT p_null     numeric   n_gs        null p-value for q-value estimation
  OUTPUT gs_pV      numeric   n_gs        p-values for gene sets
  OUTPUT gs_qV      numeric   n_gs        q-values for gene sets
  ---------------------------------------------------------------------------

  Notes:
    (1) 'testType' and 'weightType' are internally double pointers
        to character matrices.  The analyze_SP_C() uses
        strcmp(testType[0], "GSEA"), for example, to match the R character
        vector to the test conditions.
    (2) There is no parameter checking done in the C code.  All the parameter
        checking and variable type coercion are done in the R code.
    (3) The code does not respond to user interruptions (e.g., Control-Break).
    (4) The code does not free memory if the program breaks because of
        insufficient memory as indicated by Calloc().
    (5) To add your own test statistic (e.g., Wilcoxon), create a new function
        in calc_probe_stat.c and reference it in analyze_SP_C() with the
        following guidelines:
        (a) The function must have the FUNC_STAT prototype (see t_R())
        (b) The function prototype must be included in sigPathway.h
*/
void analyze_SP_C(double *Y, int *n_ps, int *ncol, double *phenotype,
		  int *n_gs, int *nsim,
		  int *nprobesV, int *indexV, int *ngroups,
		  char **testType, char **weightType, 
		  int *urp, int *verbose,

		  double *t_set, 
		  double *t_set_new, 
		  double *p_null, double *gs_pV, double *gs_qV)
{
  int adjust_t_set = 1;
  FUNC_STAT stat_fn;

  if(strcmp(testType[0], "GSEA") == 0) {
    if(*verbose == 1)
      Rprintf("Entering GSEA code...\n");

    calc_GSEA(Y, n_ps, ncol, phenotype, n_gs, nsim, nprobesV, indexV, 
	      urp, verbose,
	      t_set, t_set_new, p_null, gs_pV, gs_qV);
    /* Do not adjust t_set by the gene set size */
    adjust_t_set = 0;

  }else if(strcmp(testType[0], "NTk") == 0) {
    if(*verbose == 1)
      Rprintf("Entering NTk code...\n");
    
    if(*ngroups > 2) {
      if(*verbose == 1)
	Rprintf("Using f-statistics...\n");
      stat_fn = &f_R;
    }else if(*ngroups == 2) {
      if(*verbose == 1)
	Rprintf("Using 2-group t-statistics...\n");
      stat_fn = &t_R;
    }else {
      if(*verbose == 1)
	Rprintf("Using Pearson correlation coefficient and Fisher's z...\n");
      stat_fn = &z_R;
    }

    double *statV, *tmpV;
    statV = (double *) Calloc(*n_ps, double);
    tmpV = NULL;

    int zero = 0;
    (*stat_fn)(Y, n_ps, ncol, phenotype, ngroups, &zero, statV, tmpV, tmpV);

    calc_NGSk(statV, n_ps, n_gs, nsim, nprobesV, indexV, verbose,
	      t_set, t_set_new, p_null, gs_pV, gs_qV);

    Free(statV);

  }else if(strcmp(testType[0], "NEk") == 0) {
    if(*verbose == 1)
      Rprintf("Entering NEk code...\n");

    /*
      To speed up NEk calculations, we remove probe sets (from the expression
      matrix) that are not represented in any of the gene sets.  This requires
      us to subset the expression matrix and adjust the indices of
      the sparse gene set matrix.
    */
    int len_indexV, h, i, j;
    int *indexV_new, *ps_in_gs, nps2;

    len_indexV = 0;
    for(i = 0; i < *n_gs; i++)
      len_indexV += nprobesV[i];
    
    indexV_new = (int *) Calloc(len_indexV, int);
    ps_in_gs = (int *) Calloc(*n_ps, int);

    remove_zero_cols(nprobesV, indexV, n_gs, n_ps, indexV_new, ps_in_gs);

    nps2 = 0;
    for(i = 0; i < *n_ps; i++)
      nps2 += ps_in_gs[i];

    double *Y2;
    Y2 = (double *) Calloc(nps2 * (*ncol), double);

    h = 0;
    for(i = 0; i < *n_ps; i++) {
      if(ps_in_gs[i] == 1) {
	for(j = 0; j < *ncol; j++) {
	  Y2[h + nps2*j] = Y[i + *n_ps*j];
	}
	h++;
      }
    }

    Free(ps_in_gs);

    if(*verbose == 1)
      Rprintf("Finished subsetting expression matrix...\n");

    int wType = (strcmp(weightType[0], "variable") == 0) ? 2:1;
    if(wType == 2) {
      if(*ngroups <= 2) {
	if(*verbose == 1)
	  Rprintf("Calculating variable weights...\n");
      }else {
	if(*verbose == 1) {
	  Rprintf("Variable weights currently not implemented for f-statistics...\n");
	  Rprintf("Not calculating weights...\n");
	}
	wType = 1;
      }
    }


    if(*ngroups > 2) {
      if(*verbose == 1)
	Rprintf("Using f-statistics...\n");
      stat_fn = &f_R;
    }else if(*ngroups == 2) {
      if(*verbose == 1)
	Rprintf("Using 2-group t-statistics...\n");
      stat_fn = &t_R;
    }else {
      if(*verbose == 1)
	Rprintf("Using Pearson correlation coefficient and Fisher's z...\n");
      stat_fn = &z_R;
    }


    calc_NEk(stat_fn,
	     Y2, &nps2, ncol, phenotype, n_gs, nsim, nprobesV, indexV_new,
	     ngroups, &wType, 
	     urp, verbose,
	     t_set, t_set_new, p_null, gs_pV, gs_qV);
    
    Free(indexV_new);
    Free(Y2);

  }else if(strcmp(testType[0], "NGSk") == 0) {
    if(*verbose == 1)
      Rprintf("Entering NGSk code with user-supplied statistics...\n");
    calc_NGSk(Y, n_ps, n_gs, nsim, nprobesV, indexV, verbose,
	      t_set, t_set_new, p_null, gs_pV, gs_qV);

  }else {
    error("'%s' is not a valid test type", testType[0]);
  } 


  /* 
     Adjust t_set so that it has the MEAN of the test statistic for each
     gene set studied.  The internal code was written to do comparisons with
     the SUM of the test statistic to calculate t_set_new and downstream
     results.

     Exception: GSEA does not need this adjustment because it has a different
                definition for gene set scores
  */
  if( adjust_t_set == 1 ) {
    int i;
    for(i = 0; i < *n_gs; i++)
      t_set[i] /= nprobesV[i];
  }

  if(*verbose == 1)
    Rprintf("Finished running %s code...\n", testType[0]);

}
/************************************************************************/
