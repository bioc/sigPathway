#########################################################################
## Name:        calcProbeStat.R
## Author:      Weil Lai
## Description: Calculates the t, z, or f statistics for each probe (set)
##              based on the number of groups in the users' data
##
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## calcTStatFast() lets the user calculate the f, t, or z statistics
## quickly along with their corresponding p-values
##
## calcTNullFast() lets the user calculate the permuted distribution
## of the f, t, or z statistics
#########################################################################

calcTStatFast <- function(tab, phenotype, ngroups = 2)
{
  .checkInputs.tab.phenotype.ngroups(tab, phenotype, ngroups)

  if( ngroups > 2 )
    useCFunc <- "f_R"
  else if( ngroups == 2 )
    useCFunc <- "t_R"
  else
    useCFunc <- "z_R"

  if(ngroups > 1)
    ptype <- as.numeric(factor(phenotype)) - 1
  else
    ptype <- as.numeric(phenotype)

  res <- .C(useCFunc,
            Y = as.numeric(as.matrix(tab)),
            nrow = as.integer(nrow(tab)),
            ncol = as.integer(ncol(tab)),
            phenotype = as.numeric(ptype),
            ngroups = as.integer(ngroups),
            calcpV = as.integer(1),
            statV = numeric(nrow(tab)),
            tmpV = numeric(nrow(tab)),
            pV = numeric(nrow(tab)),
            PACKAGE = "sigPathway"
            )[7:9]

  if( ngroups > 2 )
    return(list(pval = res$pV, tstat = res$statV))
  else if( ngroups == 2 )
    return(list(pval = res$pV, tstat = res$statV, df = res$tmpV))
  else
    return(list(pval = res$pV, tstat = res$statV, rho = res$tmpV))
}

########################################################################

calcTNullFast <- function(tab, phenotype, nsim, ngroups = 2,
                          allphenotypes = FALSE)
{
  .checkInputs.tab.phenotype.ngroups(tab, phenotype, ngroups)

  if( ngroups > 2 )
    statType <- "F"
  else if( ngroups == 2 )
    statType <- "T"
  else
    statType <- "Z"

  if(ngroups > 1)
    ptype <- as.numeric(factor(phenotype)) - 1
  else
    ptype <- as.numeric(phenotype)

  res <- .C("null_R",
            statType = as.character(statType),
            allphenotypes = as.integer(allphenotypes == TRUE),
            Y = as.numeric(as.matrix(tab)),
            nrow = as.integer(nrow(tab)),
            ncol = as.integer(ncol(tab)),
            phenotype = as.numeric(ptype),
            ngroups = as.integer(ngroups),
            nsim = as.integer(nsim),
            nullV = numeric(nsim*nrow(tab)),
            PACKAGE = "sigPathway"
           )[[9]]
  dim(res) <- c(nsim, nrow(tab))
  
  return(res)
  
}

########################################################################
