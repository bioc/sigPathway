/************************************************************************
Name:        calc_probe_stat.c
Author:      Weil Lai
Description: contains probe specific statistical functions
             used for sigPathway package
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
    - Made the one and two sample cases have the same structure in their
      function declarations
    - Changed the way phenotypes are counted so that the phenotype no longer
      needs to be ordered
  - April 2006
    - Removed t_null_R(), f_null_R(), and z_null_R() in favor of the
      modular function null_c()
    - Added null_R() so that the user can call null_c() from R with .C()
  - July 14, 2006
    - By popular request, switched the ordering of the means subtracted
      in the 2 sample t-statistic calculations to match the conventions and 
      results users get when replicating the calculations with multtest 
      (and possibly other Bioconductor packages)
************************************************************************/

#include "sigPathway.h"


/************************************************************************/
/*
  Calculate null distributions of sample statistics
*/
void null_c(FUNC_STAT stat_fn, int *all_ptype,
	    double *Y, int *nrow, int *ncol,
	    double *phenotype, int *ngroups, int *nsim, double *statnV)
{
  int i, j, zero = 0;
  double *pperm, *statV, *tmpV = NULL;

  pperm = (double *) Calloc(*ncol, double);
  statV = (double *) Calloc(*nrow, double);

  i = 0;
  while(i < *nsim) {
    sampleNR_double(phenotype, pperm, ncol);

    /*
      all_ptype == 0: skip if permuted phenotype is equal to 
                      non-permuted phenotype
      all_ptype != 0: do all permutations, even if permuted phenotype is
                      equal to non-permuted phenotype
    */
    if( *all_ptype == 0 ) {
      if( memcmp(pperm, phenotype, *ncol*sizeof(double)) != 0 ) {
	(*stat_fn)(Y, nrow, ncol, pperm, ngroups, &zero, statV, tmpV, tmpV);
	/* copy permuted statistics to statnV, which has 'nsim' rows and 
	   'nrow' columns */
	for(j = 0; j < *nrow; j++)
	  statnV[i + *nsim * j] = statV[j];
	i++;
      }
    }else {
      (*stat_fn)(Y, nrow, ncol, pperm, ngroups, &zero, statV, tmpV, tmpV);
      /* copy permuted statistics to statnV, which has 'nsim' rows and 
	 'nrow' columns */
      for(j = 0; j < *nrow; j++)
	statnV[i + *nsim * j] = statV[j];
      i++;
    }
  }

  Free(pperm);
  Free(statV);
}
/************************************************************************/


/************************************************************************/
/*
  R wrapper for null_c
 */
void null_R(char **statType, int *all_ptype,
	    double *Y, int *n_ps, int *ncol,
	    double *phenotype, int *ngroups, int *nsim, double *statnV)
{
  FUNC_STAT stat_fn;
  if( strcmp(statType[0], "F") == 0 ) {
    stat_fn = &f_R;
  }else if( strcmp(statType[0], "T") == 0 ) {
    stat_fn = &t_R;
  }else if( strcmp(statType[0], "Z") == 0 ) {
    stat_fn = &z_R;
  }else {
    error("'%s' is not a supported test type", statType[0]);
    return;
  }

  null_c(stat_fn, all_ptype, Y, n_ps, ncol, phenotype, ngroups, nsim, statnV);
}
/************************************************************************/


/************************************************************************/
/*
  Calculate two sample t-statistics with unequal variances
*/
void t_R(double *Y, int *nrow, int *ncol, 
	 double *phenotype, int *ngroups, int *calc_pV,
	 double *tV, double *dfV, double *pV)
{
  int i,j,n0,n1;
  double mean[2], sse[2];

  n0 = n1 = 0;
  for( i = 0; i < *ncol; i++ ) {
    if(phenotype[i] == 0.0)
      n0++;
    else
      n1++;
  }

  for( i = 0; i < *nrow; i++ ) {
    mean[0] = mean[1] = sse[0] = sse[1] = 0.0;

    for( j = 0; j < *ncol; j++ ) {
      if(phenotype[j] == 0.0)
	mean[0] += Y[i + *nrow*j];
      else
	mean[1] += Y[i + *nrow*j];
    }

    mean[0] /= n0;
    mean[1] /= n1;

    for( j = 0; j < *ncol; j++ ) {
      if(phenotype[j] == 0.0)
	sse[0] += (Y[i + *nrow*j] - mean[0]) * (Y[i + *nrow*j] - mean[0]);
      else
	sse[1] += (Y[i + *nrow*j] - mean[1]) * (Y[i + *nrow*j] - mean[1]);
    }

    /* WL: switched the ordering of the means */
    tV[i] = (mean[1]-mean[0]) / sqrt(sse[0]/(n0)/(n0-1) + sse[1]/(n1)/(n1-1));

    if(*calc_pV == 1) {
      dfV[i] = (sse[0]/(n0)/(n0-1) + sse[1]/(n1)/(n1-1)) * 
	(sse[0]/(n0)/(n0-1) + sse[1]/(n1)/(n1-1)) /
	(sse[0]*sse[0]/(n0)/(n0)/(n0-1)/(n0-1)/(n0-1) + 
	 sse[1]*sse[1]/(n1)/(n1)/(n1-1)/(n1-1)/(n1-1));
      
      pV[i] = 2*pt(fabs(tV[i]), dfV[i], 0, 0);
    }
  }
}
/************************************************************************/


/************************************************************************/
/*
  Calculate one sample z-statistics after converting the association between
  the expression values and phenotypes (expressed as Pearson correlation
  coefficients) to Fisher's z: 0.5*log((1+rho)/(1-rho))/sqrt(n-3)
*/
void z_R(double *Y, int *nrow, int *ncol, 
	 double *phenotype, int *ngroups, int *calc_pV,
	 double *zV, double *cV, double *pV)
{
  int i,j;
  double mean_x, mean_y, sse_x, sse_y, sse_xy, rho;

  for(i = 0; i < *nrow; i++) {
    mean_x = mean_y = sse_x = sse_y = sse_xy = 0.0;

    for(j = 0; j < *ncol; j++) {
      mean_x += phenotype[j];
      mean_y += Y[i + *nrow*j];
    }
    mean_x /= *ncol;
    mean_y /= *ncol;

    for(j = 0; j < *ncol; j++) {
      sse_x += (phenotype[j] - mean_x)*(phenotype[j] - mean_x);
      sse_y += (Y[i + *nrow*j] - mean_y)*(Y[i + *nrow*j] - mean_y);
      sse_xy += (phenotype[j] - mean_x)*(Y[i + *nrow*j] - mean_y);
    }
    
    rho = sse_xy/sqrt(sse_x*sse_y);

    zV[i] = 0.5*log((1+rho)/(1-rho))*sqrt(*ncol-3);

    if(*calc_pV == 1) {
      cV[i] = rho;
      
      pV[i] = pchisq(zV[i]*zV[i], 1, 0, 0);
    }
  }

}
/************************************************************************/


/************************************************************************/
/*
  Calculate multi-sample f-statistics
*/
void f_R(double *Y, int *nrow, int *ncol,
	 double *phenotype, int *ngroups, int *calc_pV,
	 double *fV, double *tmpV, double *pV)
{
  int i, j, *ptype, *ni, df_among, df_within;
  double temp, Yi_mean, *Yig_mean, ss_among, ss_within;

  df_among = *ngroups - 1;
  df_within = *ncol - *ngroups;

  ptype = (int *) Calloc(*ncol, int);

  ni = (int *) Calloc(*ngroups, int);
  Yig_mean = (double *) Calloc(*ngroups, double);

  /* count number of samples per group */
  for(i = 0; i < *ncol; i++) {
    ptype[i] = ((int) phenotype[i]);
    ni[ptype[i]]++;
  }

  for(i = 0; i < *nrow; i++) {
    /* calculate mean */
    Yi_mean = 0.0;
    memset(Yig_mean, 0, *ngroups*sizeof(double));

    for(j = 0; j < *ncol; j++) {
      temp = Y[i + *nrow * j];
      Yi_mean += temp;
      Yig_mean[ptype[j]] += temp;
    }
    Yi_mean /= *ncol;

    for(j = 0; j < *ngroups; j++)
      Yig_mean[j] /= ni[j];

    /* calculate ss_among */
    ss_among = 0.0;
    for(j = 0; j < *ngroups; j++)
      ss_among += ni[j]*(Yig_mean[j]-Yi_mean)*(Yig_mean[j]-Yi_mean);

    /* calculate ss_within */
    ss_within = 0.0;
    for(j = 0; j < *ncol; j++) {
      temp = (Y[i + *nrow*j]-Yig_mean[ptype[j]]);
      ss_within += temp*temp;
    }

    /* calculate f-statistic and p-value */
    fV[i] = (ss_among / df_among) / (ss_within / df_within);

    if(*calc_pV == 1)
      pV[i] = pf(fV[i], df_among, df_within, 0, 0);
  }


  Free(ptype);
  Free(ni);
  Free(Yig_mean);
}
/************************************************************************/
