#########################################################################
## Name:        getPathwayStat.R
## Author:      Weil Lai
## Description: gives the user the ability to look at the mean, stdev,
##              and test statistics of individual probes for each gene set
##
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## getPathwayStatistics() and getPathwayStatistics.NGSk() lets the user
## extract probe set statistics and annotations associated with each
## selected gene set
#########################################################################

getPathwayStatistics <- function(tab, phenotype, G, index, ngroups = 2,
                                 statList = NULL, keepUnknownProbes = FALSE,
                                 annotpkg = NULL)
{
  if(!is.null(annotpkg))
    annotpkg <- .checkInputs.annotpkg(annotpkg)
  
  if(class(G) != "list" || length(G) == 0 || is.null(G[[1]]$probes))
    stop("Invalid input given for parameter 'G'\n")

  .checkInputs.tab.phenotype.ngroups(tab, phenotype, ngroups)

  if(any(is.na(index)) || !all(is.finite(index)))
    stop("The input vector 'index' cannot contain undefined values\n")
  if(any(index <= 0) || any(index > length(G)))
    stop("'index' needs refer to valid indices in 'G'\n")
  
  if(is.null(statList))
    statList <- calcTStatFast(tab, phenotype, ngroups)

  resList <- list()
  
  for(i in 1:length(index))  {

    psid <- G[[index[i]]]$probes
    ind.psid <- match(psid, rownames(tab))
    psid.nona <- psid[!is.na(ind.psid)]
    psid.na <- psid[is.na(ind.psid)]
    ind.psid <- na.omit(ind.psid)
    
    if(length(psid.na) == 0)
      psid.na <- NA
    tab.psid <- tab[ind.psid, , drop = FALSE]

    if( ngroups >= 2 )  {
      temp <- levels(factor(phenotype))
      meanM <- sdM <-
        matrix(as.numeric(NA), nrow = nrow(tab.psid), ncol = length(temp))
      for(j in 1:length(temp))  {
        meanM[,j] <- rowMeans(tab.psid[, phenotype == temp[j], drop = FALSE])
        sdM[,j] <- apply(tab.psid[, phenotype == temp[j], drop = FALSE], 1, sd)
      }
      psid.DF <- data.frame(psid.nona, meanM, sdM, statList$tstat[ind.psid],
                            statList$pval[ind.psid])
      na.DF <- as.data.frame(matrix(as.numeric(NA), nrow = length(psid.na),
                                    ncol = 2*length(temp) + 2))
      temp.DF <- data.frame(psid.na, na.DF)
      colnames(psid.DF) <- colnames(temp.DF) <-
        c("Probes", paste("Mean", temp, sep = "_"),
          paste("StDev", temp, sep = "_"),
          ifelse(ngroups == 2, "T-Statistic", "F-Statistic"), "p-value")
    }
    ##else if( ngroups == 2 )  {
    ##  mean0 <- rowMeans(tab.psid[, phenotype == 0, drop = FALSE])
    ##  sd0 <- apply(tab.psid[, phenotype == 0, drop = FALSE], 1, sd)
    ##  mean1 <- rowMeans(tab.psid[, phenotype == 1, drop = FALSE])
    ##  sd1 <- apply(tab.psid[, phenotype == 1, drop = FALSE], 1, sd)
    ##  psid.DF <- data.frame(psid.nona, mean0, mean1, sd0, sd1,
    ##                        statList$tstat[ind.psid], statList$pval[ind.psid])
    ##  temp.DF <- data.frame(psid.na, NA, NA, NA, NA, NA, NA)
    ##  colnames(psid.DF) <- colnames(temp.DF) <-
    ##    c("Probes", "Mean0", "Mean1", "StDev0", "StDev1",
    ##      "T-Statistic", "p-value")
    ##}
    else  {
      mean0 <- rowMeans(tab.psid)
      sd0 <- apply(tab.psid, 1, sd)
      psid.DF <- data.frame(psid.nona, mean0, sd0, statList$rho[ind.psid],
                            statList$tstat[ind.psid], statList$pval[ind.psid])
      temp.DF <- data.frame(psid.na, NA, NA, NA, NA, NA)
      colnames(psid.DF) <- colnames(temp.DF) <-
        c("Probes", "Mean", "StDev", "PearsonCorCoef", "FisherZ", "p-value")
    }
    
    if( keepUnknownProbes == TRUE )  {
      psid.DF <- rbind(psid.DF, temp.DF)
      rownames(psid.DF) <- c(psid.nona, psid.na)
    }

    if(!is.null(annotpkg))  {
      temp.psid <- as.character(psid.DF$Probes)
      temp.an <-
        unlist(mget(temp.psid, get(paste(annotpkg, "ACCNUM", sep = ""))))
      temp.ll <-
        unlist(mget(temp.psid, get(paste(annotpkg, "LOCUSID", sep = ""))))
      temp.gs <-
        unlist(mget(temp.psid, get(paste(annotpkg, "SYMBOL", sep = ""))))

      ## workaround for gene names split up in BioConductor 1.7 annotations
      temp.gn1 <- mget(temp.psid, get(paste(annotpkg, "GENENAME", sep = "")))
      temp.gn <- unlist(lapply(temp.gn1, paste, collapse = ";"))
      
      psid.DF <- cbind(Probes = temp.psid, AccNum = temp.an, GeneID = temp.ll,
                       Symbol = temp.gs, Name = temp.gn, psid.DF[,-1])
    }
    
    ##return(psid.DF)
    resList[[i]] <- psid.DF
  }
  
  return(resList)
}


######################################################################

getPathwayStatistics.NGSk <-
  function(statV, probeID, G, index, keepUnknownProbes = FALSE, annotpkg = NULL)
{
  if(!is.null(annotpkg))
    annotpkg <- .checkInputs.annotpkg(annotpkg)

  if(class(probeID) != "character")
    stop("Invalid input given for parameter 'probeID'\n")
  if(!is.null(names(statV)) && !identical(names(statV), probeID))
    stop("Names within 'statV' do not match 'probeID'\n")
  if(class(G) != "list" || length(G) == 0 || is.null(G[[1]]$probes))
    stop("Invalid input given for parameter 'G'\n")
  if(any(is.na(index)) || !all(is.finite(index)))
    stop("The input vector 'index' cannot contain undefined values\n")
  if(any(index <= 0) || any(index > length(G)))
    stop("'index' needs refer to valid indices in 'G'\n")

  resList <- list()
  for(i in 1:length(index))  {
    psid <- G[[index[i]]]$probes
    ind.psid <- match(psid, probeID)
    psid.nona <- psid[!is.na(ind.psid)]
    psid.na <- psid[is.na(ind.psid)]
    ind.psid <- na.omit(ind.psid)

    if(length(psid.na) == 0)
      psid.na <- NA

    psid.DF <- data.frame(psid.nona, statV[ind.psid])
    temp.DF <- data.frame(psid.na, NA)
    colnames(psid.DF) <- colnames(temp.DF) <- c("Probes", "Statistic")
    
    if( keepUnknownProbes == TRUE )  {
      psid.DF <- rbind(psid.DF, temp.DF)
      rownames(psid.DF) <- c(psid.nona, psid.na)
    }

    if(!is.null(annotpkg))  {
      temp.psid <- as.character(psid.DF$Probes)
      temp.an <-
        unlist(mget(temp.psid, get(paste(annotpkg, "ACCNUM", sep = ""))))
      temp.ll <-
        unlist(mget(temp.psid, get(paste(annotpkg, "LOCUSID", sep = ""))))
      temp.gs <-
        unlist(mget(temp.psid, get(paste(annotpkg, "SYMBOL", sep = ""))))
      
      ## workaround for annotation bug in BioConductor 1.7 annotations
      temp.gn1 <- mget(temp.psid, get(paste(annotpkg, "GENENAME", sep = "")))
      temp.gn <- unlist(lapply(temp.gn1, paste, collapse = ";"))
      
      psid.DF <- cbind(Probes = temp.psid, AccNum = temp.an, GeneID = temp.ll,
                       Symbol = temp.gs, Name = temp.gn,
                       psid.DF[,-1, drop = FALSE])
    }
    
    resList[[i]] <- psid.DF
  }
  
  return(resList)
}

########################################################################
