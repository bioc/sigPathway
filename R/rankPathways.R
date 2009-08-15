#########################################################################
## Name:        rankPathways.R
## Author:      Weil Lai
## Description: ranks the gene sets based on the magnitude of their
##              gene set statistics
##
## Change Log:
##   * April 4, 2006
##       - split sigPathway.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## rankPathways() and rankPathways.NGSk() summarizes the top pathways
## by ranking the magnitude of the gene set statistics.  rankPathways()
## is usually used for combining results from NTk and NEk calculations.
## rankPathways.NGSk() can be used for ranking pathways from one type of
## gene set calculations (e.g., GSEA, NGSk).
#########################################################################

rankPathways <- function(res.A, res.B, G, tab, phenotype, gsList, ngroups,
                         methodNames = NULL, npath = 25, allpathways = FALSE)
{
  if(class(res.A) != "list" || is.null(res.A$ngs))
    stop("Invalid input given for parameter 'res.A'\n")
  if(class(res.B) != "list" || is.null(res.B$ngs))
    stop("Invalid input given for parameter 'res.B'\n")
  if(res.A$ngs != res.B$ngs)
    stop("'res.A' and 'res.B' do not refer to the same number of gene sets.\n")
  
  if(class(G) != "list" || length(G) == 0 || is.null(G[[1]]$probes))
    stop("Invalid input given for parameter 'G'\n")

  .checkInputs.gsList(gsList)
  
  if(is.na(npath) || !is.finite(npath) || npath <= 0 ||
     npath > length(res.A$t.set.new) )
    stop(paste("'npath' must be a positive integer less than or equal to the number of pathways studied (n = ", length(gsList$nprobesV), ").\n", sep = ""))
  if(!is.logical(allpathways))
    stop("'allpathways' must be set to either TRUE or FALSE.\n")

  .checkInputs.tab.phenotype.ngroups(tab, phenotype, ngroups)
  
  ngs <- length(gsList$nprobesV)
  
  mapping <- integer(ngs)
  for( i in 1:ngs )
    mapping[i] <- min(which(res.A$t.set == res.A$t.set[i]))
 
  unique.set.index <- unique(mapping)
  ngs.unique <- length(unique.set.index)

  r.mA <- r.mB <- integer(ngs)

  r.mA[unique.set.index] <-
    ngs.unique + 1 - rank(abs(res.A$t.set.new[unique.set.index]))
  r.mB[unique.set.index] <-
    ngs.unique + 1 - rank(abs(res.B$t.set.new[unique.set.index]))

  r.mA <- r.mA[mapping]
  r.mB <- r.mB[mapping]
  r0 <- r.mA + r.mB

  if(!allpathways)  {
    id.A <- order(abs(res.A$t.set.new), decreasing = TRUE)[1:npath]
    id.B <- order(abs(res.B$t.set.new), decreasing = TRUE)[1:npath]
    id <- unique(c(id.A,id.B))
    idnew <- id[order(r0[id])]
  }else  {
    idnew <- order(r0)[1:npath]
  }

  if(ngroups == 2)  {
    temp <- levels(factor(phenotype))
    meanV.0 <- rowMeans(tab[,phenotype == temp[1]])
    meanV.1 <- rowMeans(tab[,phenotype == temp[2]])
    temp.up <- 1*(meanV.1 > meanV.0)
    
    percentUp <- numeric(length(idnew))
    for(i in 1:length(idnew))  {
      psid <- G[[gsList$indGused[idnew[i]]]]$probes
      ind.psid <- match(psid, rownames(tab))
      psid.nona <- psid[!is.na(ind.psid)]
      psid.na <- psid[is.na(ind.psid)]
      ind.psid <- na.omit(ind.psid)
      
      percentUp[i] <- round(100*sum(temp.up[ind.psid])/length(ind.psid), 2)
    }
  }
  
  temp1 <- cbind(res.A$t.set.new[idnew], res.B$t.set.new[idnew])
  temp1 <- round(temp1, 2)
  
  set.name1 <- set.name2 <- character(0)
  for(i in idnew)  {
    set.name1 <- c(set.name1, G[[gsList$indGused[i]]]$src)
    temp.name <- G[[gsList$indGused[i]]]$title
    if(length(temp.name) == 0)
      set.name2 <- c(set.name2, attr(temp.name, "Term"))
    else
      set.name2 <- c(set.name2, temp.name)
  }
  
  select.set.size <- gsList$nprobesV[idnew]
  
  temp2 <- cbind(res.A$q.value[idnew], res.B$q.value[idnew])
    
  temp3 <- round(cbind(r.mA[idnew], r.mB[idnew]), 2)


  if(ngroups == 2)  {
    tabpath <- data.frame(gsList$indGused[idnew], I(set.name1), I(set.name2),
                          select.set.size, percentUp,
                          temp1[,1], temp2[,1], temp3[,1],
                          temp1[,2], temp2[,2], temp3[,2])
    
    if( is.null(methodNames) || length(methodNames) != 2 )  {
      colnames(tabpath) <- c("IndexG", "Gene Set Category", "Pathway",
                             "Set Size", "Percent Up",
                             "MethodA Stat", "MethodA q-value", "MethodA Rank",
                             "MethodB Stat", "MethodB q-value", "MethodB Rank")
    }else  {
      colnames(tabpath)[1:5] <- c("IndexG", "Gene Set Category", "Pathway",
                                  "Set Size", "Percent Up")
      suffixes <- c("Stat", "q-value", "Rank")
      colnames(tabpath)[6:8] <- paste(methodNames[1], suffixes, sep = " ")
      colnames(tabpath)[9:11] <- paste(methodNames[2], suffixes, sep = " ")
    }
  }else  {
    tabpath <- data.frame(gsList$indGused[idnew], I(set.name1), I(set.name2),
                          select.set.size,
                          temp1[,1], temp2[,1], temp3[,1],
                          temp1[,2], temp2[,2], temp3[,2])
    
    if( is.null(methodNames) || length(methodNames) != 2 )  {
      colnames(tabpath) <- c("IndexG", "Gene Set Category", "Pathway",
                             "Set Size",
                             "MethodA Stat", "MethodA q-value", "MethodA Rank",
                             "MethodB Stat", "MethodB q-value", "MethodB Rank")
    }else  {
      colnames(tabpath)[1:4] <- c("IndexG", "Gene Set Category", "Pathway",
                                  "Set Size")
      suffixes <- c("Stat", "q-value", "Rank")
      colnames(tabpath)[5:7] <- paste(methodNames[1], suffixes, sep = " ")
      colnames(tabpath)[8:10] <- paste(methodNames[2], suffixes, sep = " ")
    }
  }    

  return(tabpath)
}

######################################################################

rankPathways.NGSk <- function(res.NGSk, G, gsList, methodName = "NGSk",
                              npath = 25)
{
  if(class(res.NGSk) != "list" || is.null(res.NGSk$ngs))
    stop("Invalid input given for parameter 'res.NGSk'\n")
  
  if(class(G) != "list" || length(G) == 0 || is.null(G[[1]]$probes))
    stop("Invalid input given for parameter 'G'\n")

  .checkInputs.gsList(gsList)
  
  if(is.na(npath) || !is.finite(npath) || npath <= 0 ||
     npath > length(res.NGSk$t.set.new) )
    stop("'npath' must be a positive integer.\n")

  ngs <- length(gsList$nprobesV)
  
  mapping <- integer(ngs)
  for( i in 1:ngs )
    mapping[i] <- min(which(res.NGSk$t.set == res.NGSk$t.set[i]))
 
  unique.set.index <- unique(mapping)
  ngs.unique <- length(unique.set.index)
  
  r.mNGSk <- integer(ngs)
  r.mNGSk[unique.set.index] <-
    ngs.unique + 1 - rank(abs(res.NGSk$t.set.new[unique.set.index]))
  r.mNGSk <- r.mNGSk[mapping]

  idnew <- order(r.mNGSk)[1:npath]
  
  temp1 <- round(res.NGSk$t.set.new[idnew], 2)
  
  set.name1 <- set.name2 <- character(0)
  for(i in idnew)  {
    set.name1 <- c(set.name1, G[[gsList$indGused[i]]]$src)
    temp.name <- G[[gsList$indGused[i]]]$title
    if(length(temp.name) == 0)
      set.name2 <- c(set.name2, attr(temp.name, "Term"))
    else
      set.name2 <- c(set.name2, temp.name)
  }

  select.set.size <- gsList$nprobesV[idnew]
  temp2 <- res.NGSk$q.value[idnew]
  temp3 <- round(r.mNGSk[idnew], 2)

  tabpath <- data.frame(gsList$indGused[idnew], I(set.name1), I(set.name2),
                        select.set.size, temp1, temp2, temp3)
  colnames(tabpath)[1:4] <-
    c("IndexG", "Gene Set Category", "Pathway", "Set Size")
  suffixes <- c("Stat", "q-value", "Rank")
  colnames(tabpath)[5:7] <- paste(methodName, suffixes, sep = " ")

  return(tabpath)
}
  
######################################################################
