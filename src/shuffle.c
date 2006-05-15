/************************************************************************
Name:        shuffle.c
Author:      Weil Lai
Description: contains sampling without replacement functions
             used for sigPathway package
Change Logs:
  - March 30, 2006
    - Initial version created to modularize sigPathway.c v. 1.1-3
      for easier maintenance
    - Changed m to *m so that the user may call the below functions in R
      with .C()
************************************************************************/

#include "sigPathway.h"

/*
  Random sampling of integers without replacement
  Use set.seed(...) in R to set initial seed
*/
void sampleNR_int(int *inV, int *outV, int *m)
{
  int i, j, temp;
  Memcpy(outV, inV, *m);
  GetRNGstate();
  for( i = 0; i < *m-1; i++ ) {
    j = i + (int)(*m-i)*unif_rand();
    temp = outV[j];
    outV[j] = outV[i];
    outV[i] = temp;
  }
  PutRNGstate();
}


/*
  Random sampling of doubles without replacement
  Use set.seed(...) in R to set initial seed
*/
void sampleNR_double(double *inV, double *outV, int *m)
{
  int i, j;
  double temp;
  Memcpy(outV, inV, *m);
  GetRNGstate();
  for( i = 0; i < *m-1; i++ ) {
    j = i + (int)(*m-i)*unif_rand();
    temp = outV[j];
    outV[j] = outV[i];
    outV[i] = temp;
  }
  PutRNGstate();
}
