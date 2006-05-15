#########################################################################
## Name:        checkInputs.R
## Author:      Weil Lai
## Description: Internal R functions which check the users' parameter
##              inputs
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## The below functions are used internally to check the parameter inputs
##
#########################################################################

.checkInputs.tab.phenotype.ngroups <- function(tab, phenotype, ngroups)  {
  if(any(is.na(tab)) || !all(is.finite(as.matrix(tab))))
    stop("'tab' cannot contain undefined values\n")
  .checkInputs.phenotype.ngroups(phenotype, ngroups)
}

########################################################################
.checkInputs.phenotype.ngroups <- function(phenotype, ngroups)  {
  if(any(is.na(phenotype)))
    stop("'phenotype' cannot contain undefined values\n")
  if(is.na(ngroups) || !is.finite(ngroups) || !(ngroups >= 1))
    stop("'ngroups' must be greater than or equal to 1.\n")
  if(is.factor(phenotype))
    stop("'phenotype' cannot be a vector of factors\n")
  
  if(ngroups >= 2)  {
    temp <- levels(factor(phenotype))
    if(length(temp) != ngroups)  {
      stop("Mismatch in 'ngroups' (", ngroups, ") and the number of groups (",
           length(temp), ") in 'phenotype'\n")
    }
  }
  ##else if(ngroups == 2)  {
  ##  temp <- levels(factor(phenotype))
  ##  if(length(temp) != 2 || any(is.na(match(temp, c("0","1")))))
  ##    stop("'phenotype' can only consist of 0 and 1s\n")
  ##}
}

########################################################################
.checkInputs.gsList <- function(gsList)  {
  if(length(gsList$nprobesV) == 0 || length(gsList$indexV) == 0)
    stop("Invalid input given for parameter 'gsList'\n")
}

########################################################################
.checkInputs.nsim <- function(nsim)  {
  if(is.na(nsim) || !is.finite(nsim) || nsim <= 0)
    stop("'nsim' must be a positive integer.\n")
}

########################################################################
.checkInputs.annotpkg <- function(annotpkg)  {
  if(require(annotpkg, character.only = TRUE))  {
    return(annotpkg)
  }else  {
    cat("Warning: package", annotpkg, "not found.  Setting 'annotpkg' to NULL\n")
    return(NULL)
  }
}

########################################################################
