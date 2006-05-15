#########################################################################
## Name:        estimateNumPerm.R
## Author:      Weil Lai
## Description: Estimates the number of permutations needed for sample
##              (i.e., column) permutation
##
## Change Log:
##   - April 13, 2006
##     - initial version
#########################################################################

#########################################################################
## estimateNumPerm() counts the total number of unique permutations of
## the phenotype vector
#########################################################################

estimateNumPerm <- function(phenotype, ngroups)  {
  .checkInputs.phenotype.ngroups(phenotype, ngroups)
  
  if(ngroups == 1)
    ngroups <- length(unique(phenotype))

  res <- .C("count_perm",
            n = as.integer(length(phenotype)),
            phenotype = as.numeric(as.integer(factor(phenotype)) - 1),
            ngroups = as.integer(ngroups),
            ncperm = numeric(1),
            PACKAGE = "sigPathway"
            )[[4]]
  return(res)
}
