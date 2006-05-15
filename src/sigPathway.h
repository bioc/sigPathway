/************************************************************************
Name:        sigPathway.c
Author:      Weil Lai
Description: C header file for all C functions used in sigPathway package
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
  - April 2006
    - Changes made to function headers as required by changes made in the
      c code
************************************************************************/

#ifndef SIGPATHWAY_H
#define SIGPATHWAY_H


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>


typedef void (*FUNC_STAT)(double *, int *, int *, 
			  double *, int *, int *,
			  double *, double *, double *);

/* calc_gene_set_internal.c */
void calc_internal(int *verbose, int *n_gs, int *nsim,
		   double *t_set, double *t_set_null, double *t_set_new);
void calc_internal2(int *verbose, int *n_gs, int *nsim,
		    double *t_set_null, double *t_set_new,
		    double *t_set_new_null, double *p_null,
		    double *pV, double *qV);
void rank_avg2(int *m, double *inV, double *rank);
void maj(int *n, double *x, double *y, double *resV);
void calcQVFast(double *tsn, int *ngs,
		double *p_null,
		double *tsn_null, int *nsim,
		double *qV);

/* calc_gene_set_stat.c */
void calc_GSEA(double *Y, int *n_ps, int *ncol, double *phenotype,
	       int *n_gs, int *nsim,
	       int *nprobesV, int *indexV,
	       int *urp, int *verbose,

	       double *es,
	       double *nes,
	       double *p_null, double *gs_pV, double *gs_qV);

void calc_NGSk(double *tV, int *n_ps, int *n_gs, int *nsim,
	       int *nprobesV, int *indexV, int *verbose,

	       double *t_set,
	       double *t_set_new,
	       double *p_null, double *gs_pV, double *gs_qV);

void calc_NEk(FUNC_STAT stat_fn,
	      double *Y, int *n_ps, int *ncol, double *phenotype,
	      int *n_gs, int *nsim, 
	      int *nprobesV, int *indexV,
	      int *ngroups, int *wType,
	      int *urp, int *verbose,
	      
	      double *t_set,
	      double *t_set_new,
	      double *p_null, double *gs_pV, double *gs_qV);

/* calc_probe_stat.c */
void null_c(FUNC_STAT stat_fn, int *all_ptype,
	    double *Y, int *nrow, int *ncol,
	    double *phenotype, int *ngroups, int *nsim, double *statnV);
void null_R(char **statType, int *all_ptype,
	    double *Y, int *n_ps, int *ncol,
	    double *phenotype, int *ngroups, int *nsim, double *statnV);
void t_R(double *Y, int *nrow, int *ncol, 
	 double *phenotype, int *ngroups, int *calc_pV,
	 double *tV, double *dfV, double *pV);
void z_R(double *Y, int *nrow, int *ncol,
	 double *phenotype, int *ngroups, int *calc_pV,
	 double *zV, double *cV, double *pV);
void f_R(double *Y, int *nrow, int *ncol,
	 double *phenotype, int *ngroups, int *calc_pV,
	 double *fV, double *tmpV, double *pV);

/* calc_weights.c */
void internal_weight(double *e_vectors, double *e_values, int *np,
		     double *qq0, double *lambda, int *only_weights,
		     double *weight, double *minweight);
void calcWeight_common(double *v, int np, int nused, int *verbose, double *w);
void calcWeights2Groups(double *Y, int *n_ps, int *ncol, double *phenotype,
			int *indexV,
			int *nprobesV, int *n_gs, int *verbose,
			double *weightsV);
void calcWeights1Group(FUNC_STAT stat_fn,
		       double *Y, int *n_ps, int *ncol, double *phenotype,
		       int *nsim,
		       int *indexV,
		       int *nprobesV, int *n_gs, int *verbose,
		       double *weightsV);
void calcWeights(FUNC_STAT stat_fn,
		 double *Y, int *n_ps, int *ncol, double *phenotype,
		 int *n_gs, int *ngroups, int *nsim,
		 int *nprobesV, int *indexV, int *verbose, double *weightsV);

/* matrix_fn.c */
void matprod(double *x, int nrx, int ncx,
	     double *y, int nry, int ncy, double *z);
void eigen(double *xin, int *nn, double *values, double *vectors);
double covar(double *x, double *y, int n);
void covar_mat(double *x, int *nr, int *nc, int *fill, double *output);
void remove_zero_cols(int *nprobesV, int *indexV, int *n_gs, int *n_ps,
		      int *indexV_new, int *ps_in_gs);

/* shuffle.c */
void sampleNR_int(int *inV, int *outV, int *m);
void sampleNR_double(double *inV, double *outV, int *m);

/* permute.c */
void get_next_perm(int *n, double *pperm);
void count_perm(int *n, double *phenotype, int *ngroups, double *ncperm);
double factorial(int N);


#endif
