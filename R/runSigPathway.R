#########################################################################
## Name:        runSigPathway.R
## Author:      Weil Lai
## Description: convenient all-in-one wrapper function for beginning users
##              to analyze their microarray data with sigPathway
##
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## runSigPathway() is a large wrapper function for selectGeneSets(),
## calculate.NTk(), calculate.NEk(), rankPathways(), calcTStatFast(),
## and getPathwayStatistics().  Its purpose is to let the user run
## pathway analysis without the fuss of typing out function arguments
## for multiple intermediate functions()
##
#########################################################################

runSigPathway <- function(G, minNPS = 20, maxNPS = 500,
                          tab, phenotype, nsim = 1000,
                          weightType = c("constant", "variable"), ngroups = 2,
                          npath = 25, verbose = FALSE, allpathways = FALSE,
                          annotpkg = NULL, alwaysUseRandomPerm = FALSE)
{
  weightType <- match.arg(weightType)
  if(!is.null(annotpkg))
    annotpkg <- .checkInputs.annotpkg(annotpkg)

  cat("Selecting the gene sets\n")
  gsList <- selectGeneSets(G, rownames(tab), minNPS, maxNPS)

  cat("Calculating NTk statistics for each selected gene set\n")
  list.NTk <- calculate.NTk(tab, phenotype, gsList, nsim, ngroups,
                            verbose, alwaysUseRandomPerm)

  if(weightType == "constant")  {
    cat("Calculating NEk statistics for each selected gene set\n")
    methodNames <- c("NTk", "NEk")
  }else  {
    cat("Calculating weighted NEk statistics for each selected gene set\n")
    methodNames <- c("NTk", "NEk*")
  }
  
  list.NEk <- calculate.NEk(tab, phenotype, gsList, nsim, weightType,
                            ngroups, verbose, alwaysUseRandomPerm)
  
  if(!allpathways)  {
    cat("Summarizing the top", npath, "pathways from each statistic\n",
        sep = " ")
  }else  {
    cat("Summarizing the top", npath, "pathways\n", sep = " ")
  }
  
  df.pathways <- rankPathways(list.NTk, list.NEk, G, tab, phenotype, gsList,
                              ngroups, methodNames, npath, allpathways)

  sList <- calcTStatFast(tab, phenotype, ngroups)
  list.gPS <- getPathwayStatistics(tab, phenotype, G, df.pathways[,"IndexG"],
                                   ngroups, sList, FALSE, annotpkg)

  cat("Done!  Use the writeSigPathway() function to write results to HTML\n")

  list.parameters <-
    list(nprobes = nrow(tab), nsamples = ncol(tab),
         phenotype = phenotype, ngroups = ngroups,
         minNPS = minNPS, maxNPS = maxNPS, ngs = list.NTk$ngs,
         nsim.NTk = list.NTk$nsim, nsim.NEk = list.NEk$nsim,
         weightType = weightType,
         annotpkg = annotpkg, npath = npath, allpathways = allpathways)
         
  return(list(
              gsList = gsList,
              list.NTk = list.NTk,
              list.NEk = list.NEk,
              df.pathways = df.pathways,
              list.gPS = list.gPS,
              parameters = list.parameters
              )
         )
}

########################################################################
