/************************************************************************
Name:        permute.c
Author:      Weil Lai
Description: Contains functions to permute phenotypes in lexicographic
             order and count the total number of unique permutations
Change Logs:
  - April 13, 2006
    - Initial version
************************************************************************/

#include "sigPathway.h"

/************************************************************************/
/*
  Based on the pseudocode at
  http://www.cut-the-knot.org/do_you_know/AllPerm.shtml

  The website's author (Alexander Bogomolny) referenced this algorithm
  from page 71 of

  E.W. Dijkstra, _A Discipline of Programming_, Prentice-Hall, 1997

  Take a look at calc_NEk() in calc_gene_set_stat.c for an example on how
  get_next_perm() is used
*/
void get_next_perm(int *n, double *pperm)
{
  int i, j;
  double temp_double;

  i = *n-1;
  j = *n;

  while(pperm[i-1] >= pperm[i])
    i--;

  while(pperm[j-1] <= pperm[i-1])
    j--;

  temp_double = pperm[i-1];
  pperm[i-1] = pperm[j-1];
  pperm[j-1] = temp_double;

  i++;
  j = *n;

  while(i < j) {
    temp_double = pperm[i-1];
    pperm[i-1] = pperm[j-1];
    pperm[j-1] = temp_double;
    
    i++;
    j--;
  }

}
/************************************************************************/


/************************************************************************/
/*
  Computes the number of unique permutations based on phenotype and the
  number of groups.

  e.g., 

  Phenotype               ngroups            ncperm
  -----------------------------------------------------------------------
  c(0.5,0.7,0.8,0.9)         1 (continuous)  4! = 24
  c(0.5,0.7,0.8,0.8)         1 (continuous)  4!/2! = 12
  c(0,0,0,1,1,1)             2               6!/3!/3! = 20
  c(0,1,2,0,1,2,0,1,2)       3               9!/3!/3!/3! = 1680

*/
void count_perm(int *n, double *phenotype, int *ngroups, double *ncperm)
{
  int i, j, *nV;
  double numerator, temp;

  nV = (int *) Calloc(*ngroups, int);
  for(i = 0; i < *n; i++)
    nV[(int) phenotype[i]]++;

  temp = numerator = 1.0;
  for(i = 0; i < *ngroups; i++) {
    for(j = 1; j <= nV[i]; j++) {
      temp *= numerator / j;
      numerator += 1.0;
    }
  }

  *ncperm = temp;
  Free(nV);
}
/************************************************************************/


/************************************************************************/
/* Currently unused in the program but kept for potential future use */
double factorial(int N)
{
  int i;
  double ans = 1.0;
  for(i = 1; i <= N; i++)
    ans *= (double) i;
  return(ans);
}
/************************************************************************/
