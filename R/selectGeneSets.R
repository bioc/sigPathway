#########################################################################
## Name:        selectGeneSets.R
## Author:      Weil Lai
## Description: Selects pathways to analyze based on the probes and the
##              min/max pathway size used
##
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## selectGeneSets() identifies the gene sets that are represented in
## the user's expression matrix and encodes the results into a sparse
## matrix.
#########################################################################

selectGeneSets <- function(G, probeID, minNPS = 20, maxNPS = 500)
{
  if(class(G) != "list" || length(G) == 0 || is.null(G[[1]]$probes))
    stop("Invalid input given for parameter 'G'\n")
  if(class(minNPS) != "numeric" || is.na(minNPS) || minNPS <= 0)
    stop("Invalid input given for parameter 'minNPS'\n")
  if(class(maxNPS) != "numeric" || is.na(maxNPS) || maxNPS <= 0)
    stop("Invalid input given for parameter 'maxNPS'\n")
  if(class(probeID) != "character" || length(probeID) < minNPS)
    stop("Invalid input given for parameter 'probeID'\n")

  if(maxNPS <= minNPS)
    stop("minNPS needs to be less than maxNPS\n")
  
  nset <- length(G)
  
  set.sizeAP <- set.sizeCP <- integer(nset)
  for (i in 1:nset) {
    set.sizeAP[i] <- length(G[[i]]$probes)
    set.sizeCP[i] <- sum(!is.na(match(G[[i]]$probes, probeID)))
  }

  ind.pass <- which(set.sizeCP >= minNPS & set.sizeAP < maxNPS)
  ngs <- length(ind.pass)

  if(ngs == 0)  {
    stop("No gene sets found with the given parameters.\nPlease check that 'probeID' contains a character vector of probe IDs,\nthat minNPS is less than maxNPS,\nand that the values for minNPS and maxNPS are not too strigent.")
  }
  
  indexV <- NULL
  nprobesV <- integer(ngs)
  for(i in 1:ngs)  {
    temp.ind <- na.omit(match(G[[ind.pass[i]]]$probes, probeID)-1)
    indexV <- c(indexV, temp.ind)
    nprobesV[i] <- length(temp.ind)
  }

  return(list(nprobesV = nprobesV,
              indexV = indexV,
              indGused = ind.pass))
}

########################################################################
