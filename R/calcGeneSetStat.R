#########################################################################
## Name:        calcGeneSetStat.R
## Author:      Weil Lai
## Description: Calculates the gene set statistics (wrapper functions
##              for the C code, which does most of the work)
##
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## .calcSP() is the internal R wrapper function which sends parameters
## to analyze_SP_C() of sigPathway.c
##
## The following R functions are wrapper functions of .calcSP():
##   calculate.GSEA(), calculate.NGSk(), calculate.NTk(), calculate.NEk()
#########################################################################

.calcSP <- function(tab, phenotype, gsList, nsim, testType, weightType,
                    ngroups, verbose, alwaysUseRandomPerm)
{
  .checkInputs.tab.phenotype.ngroups(tab, phenotype, ngroups)
  .checkInputs.gsList(gsList)
  .checkInputs.nsim(nsim)

  urp <- 1

  if( (testType == "NEk" | testType == "GSEA") &
       alwaysUseRandomPerm == FALSE )  {
    ## check to see if the number of unique permutations is less than
    ## the number of user-specified permutations
    np <- estimateNumPerm(phenotype, ngroups)
    if(nsim > np)  {
      urp <- 0
      nsim <- np - 1
      temp.text <-
        c("'nsim' is greater than the number of unique permutations\n",
          paste("Changing 'nsim' to ", nsim,
                ", excluding the unpermuted case\n", sep = "")
          )
      
      cat(temp.text)
    }
  }

  if(ngroups > 1)
    ptype <- as.numeric(factor(phenotype)) - 1
  else
    ptype <- as.numeric(phenotype)
  
  ngs <- length(gsList$nprobesV)
  res <- .C("analyze_SP_C",
            Y = as.numeric(as.matrix(tab)),
            nps = as.integer(nrow(tab)),
            ncol = as.integer(ncol(tab)),
            phenotype = as.numeric(ptype),
            ngs = as.integer(ngs),
            nsim = as.integer(nsim),
            nprobesV = as.integer(gsList$nprobesV),
            indexV = as.integer(gsList$indexV),
            ngroups = as.integer(ngroups),
            testType = as.character(testType),
            weightType = as.character(weightType),
            urp = as.integer(urp),
            verbose = as.integer(verbose == TRUE),
            
            t.set = numeric(ngs),
            t.set.new = numeric(ngs),
            p.null = numeric(1),
            p.value = numeric(ngs),
            q.value = numeric(ngs),
            PACKAGE = "sigPathway"
            )[c(5:6, 14:18)]
  return(res)
}

######################################################################

calculate.GSEA <- function(tab, phenotype, gsList, nsim = 1000,
                           verbose = FALSE, alwaysUseRandomPerm = FALSE)
{
  return(.calcSP(tab, phenotype, gsList, nsim,
                 "GSEA", "constant", 2, verbose, alwaysUseRandomPerm)
         )
}  

######################################################################

calculate.NGSk <- function(statV, gsList, nsim = 1000, verbose = FALSE,
                           alwaysUseRandomPerm = FALSE)  {
  
  if(any(is.na(statV)) || !all(is.finite(statV)))
    stop("The input vector 'statV' should only contain finite values\n")

  statV2 <- matrix(statV, ncol = 1)
  return(.calcSP(statV2, c(0,1), gsList, nsim,
                 "NGSk", "constant", 2, verbose, alwaysUseRandomPerm)
         )
}

######################################################################

calculate.NTk <- function(tab, phenotype, gsList, nsim = 1000, ngroups = 2,
                          verbose = FALSE, alwaysUseRandomPerm = FALSE)
{
  return(.calcSP(tab, phenotype, gsList, nsim,
                 "NTk", "constant", ngroups, verbose, alwaysUseRandomPerm)
         )
}

######################################################################

calculate.NEk <- function(tab, phenotype, gsList, nsim = 1000,
                          weightType = c("constant", "variable"), ngroups = 2,
                          verbose = FALSE, alwaysUseRandomPerm = FALSE)
{
  weightType <- match.arg(weightType)

  return(.calcSP(tab, phenotype, gsList, nsim,
                 "NEk", weightType, ngroups, verbose, alwaysUseRandomPerm)
         )
}

########################################################################
