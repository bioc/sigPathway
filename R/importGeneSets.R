#########################################################################
## Name:        importGeneSets.R
## Author:      Weil Lai
## Description: Converts gene sets stored in GMT, GMX, GRP, and XML
##              formats to sigPathway's gene set list format
##
## Change Log:
##   * June 13, 2007
##     - importGeneSets() added
##   * June 7, 2007
##     - gmtToG(), gmxToG(), grpToG(), and xmlToG() added to sigPathway
#########################################################################

importGeneSets <- function(fileNames, verbose = TRUE)  {
  ## Check if the files do really exist
  tmp <- file.exists(fileNames)
  if(any(tmp == FALSE))  {
    stop(paste("File(s) ", fileNames[tmp == FALSE],
               " does/do not exist.\n", sep = ""))
  }
  
  ## Check the file name extensions to see whether they are either
  ## GMT, GMX, GRP, or XML
  extType <- numeric(length(fileNames))
  extType[grep("[.][gG][mM][tT]$", fileNames)] <- 1
  extType[grep("[.][gG][mM][xX]$", fileNames)] <- 2
  extType[grep("[.][gG][rR][pP]$", fileNames)] <- 3
  extType[grep("[.][xX][mM][lL]$", fileNames)] <- 4

  if(any(extType == 0))  {
    stop(paste("File(s) ", fileNames[extType == 0],
               " has/have the wrong extension(s).\n", sep = ""))
  }
  
  ## Read in gene sets based on their file extensions
  G <- list()
  for(i in 1:length(fileNames))  {
    fn <- fileNames[i]
    if(verbose)
      cat("Reading ", fn, "\n", sep = "")

    G.output <- switch(extType[i], gmtToG(fn, FALSE), gmxToG(fn, FALSE),
                       grpToG(fn, FALSE), xmlToG(fn, FALSE))
    G <- c(G, G.output)
  }

  if(verbose)
    cat("Imported ", length(G), " gene set(s)\n", sep = "")
  return(invisible(G))
}

gmtToG <- function(fileNames, verbose = TRUE)  {
  G <- list()
  for(fn in fileNames)  {
    if(verbose)
      cat("Reading ", fn, "\n", sep = "")
    tmpLines <- readLines(fn)
    ngs <- length(tmpLines)
    G.output <- vector("list", length(ngs))
    
    for(i in 1:ngs)  {
      tmp <- unlist(strsplit(tmpLines[i], split = "\t"))
      gs <- unique(tmp[3:length(tmp)])
      gs <- gs[gs != ""]  ## this is for cases where the GMT file was saved with
                          ## extra tabs from Excel
      G.output[[i]] <- list(src = tmp[1],
                            title = ifelse(tmp[2] == "na" |
                              tmp[2] == "NA", "", tmp[2]),
                            probes = gs)
    }
    G <- c(G, G.output)
  }
  if(verbose)
    cat("Imported ", length(G), " gene set(s)\n", sep = "")
  return(invisible(G))
}

gmxToG <- function(fileNames, verbose = TRUE)  {
  G <- list()
  for(fn in fileNames)  {
    if(verbose)
      cat("Reading ", fn, "\n", sep = "")
    tmpMat <- read.delim(fn, header = FALSE, as.is = TRUE)
    ngs <- ncol(tmpMat)
    nr <- nrow(tmpMat)
    
    G.output <- vector("list", length = ngs)
    for(i in 1:ngs)  {
      gs <- unique(tmpMat[3:nr,i])
      gs <- gs[gs != ""]
      G.output[[i]] <- list(src = tmpMat[1,i],
                            title = ifelse(tmpMat[2,i] == "na" |
                              tmpMat[2,i] == "NA", "", tmpMat[2,i]),
                            probes = gs)
    }
    G <- c(G, G.output)
  }
  if(verbose)
    cat("Imported ", length(G), " gene set(s)\n", sep = "")
  return(invisible(G))
}

grpToG <- function(fileNames, verbose = TRUE)  {
  G <- list()
  for(fn in fileNames)  {
    if(verbose)
      cat("Reading ", fn, "\n", sep = "")
    tmpLines <- unique(readLines(fn))
    
    ## GSEA wiki suggests that there could be comments embedded in the GRP file
    idxComment <- grep("^#", tmpLines)
    if(length(idxComment) > 0)
      tmpLines <- tmpLines[-idxComment]
    
    srcName <- sub("[.]grp", "", fileNames)
    G.output <- list()
    G.output[[1]] <- list(src = srcName, title = "", probes = tmpLines)
    G <- c(G, G.output)
  }
  if(verbose)
    cat("Imported ", length(G), " gene set(s)\n", sep = "")
  return(invisible(G))
}

xmlToG <- function(fileNames, verbose = TRUE)  {
  if(require(XML))  {
    G <- list()
    for(fn in fileNames)  {
      if(verbose)
        cat("Reading ", fn, "\n", sep = "")
      doc <- xmlTreeParse(fn)
      r <- xmlRoot(doc)
      ngs <- xmlSize(r)
      
      G.output <- vector("list", length = ngs)
      for(i in 1:ngs)  {
        attrs <-
          xmlAttrs(r[[i]])[c("STANDARD_NAME", "SYSTEMATIC_NAME", "ORGANISM",
                             "CONTRIBUTOR", "DESCRIPTION_BRIEF", "MEMBERS")]
        names(attrs) <- NULL
        G.output[[i]] <-
          list(src = attrs[1],
               title = attrs[5],
               organism = attrs[3],
               name = attrs[2],
               contributor = attrs[4],
               probes = unique(unlist(strsplit(attrs[6], split=","))))    
      }
      G <- c(G, G.output)
    }
    
    if(verbose)
      cat("Imported ", length(G), " gene set(s)\n", sep = "")
    return(invisible(G))
  }else  {
    stop("Unable to load XML package.  Please check that XML has been installed on your system.")
  }
}
