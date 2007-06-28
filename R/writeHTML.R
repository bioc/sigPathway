#########################################################################
## Name:        writeHTML.R
## Author:      Weil Lai
## Description: Functions to write results of pathway analysis to HTML
##              documents
##
## Change Log:
##   * April 4, 2006
##       - new functions: .convertDF2HTML(), writeSP(), writeSigPathway()
#########################################################################

#########################################################################
## (1) .convertDF2HTML() converts data frames into HTML tables
## (2) writeSP() creates multiple HTML documents representing the list of
##     top pathways and their associated probe sets
## (3) writeSigPathway() is a wrapper function for writeSP() using
##     results obtained from runSigPathway()
#########################################################################

.convertDF2HTML <- function(inputDF, digits)  {
  if(!is.data.frame(inputDF))
    stop("'inputDF' must be a data frame\n")

  nr <- nrow(inputDF)
  nc <- ncol(inputDF)
  
  output.body <- character(nr)
  for(i in 1:nc)  {
    temp <- inputDF[[i]]
    if(is.factor(temp))
      temp <- as.character(temp)

    if(!is.character(temp))  {
      temp <- formatC(temp, digits = digits[i], format = "f")
      output.body <- paste(output.body,
                           paste("<TD align=\"right\">", temp, "</TD>"),
                           sep = "")
    }else  {
      output.body <- paste(output.body,
                           paste("<TD>", temp, "</TD>"),
                           sep = "")
    }  
  }

  output.rownames <- paste("<TD align=\"right\">", rownames(inputDF), "</TD>")
  output.colnames <-
    paste(paste("<TH>", c("", colnames(inputDF)), "</TH>"), collapse = "")

  output <-
    c("<TABLE border=1>",
      paste("<TR>",
          c(output.colnames, paste(output.rownames, output.body)),
            "</TR>"),
      "</TABLE>"
      )

  return(output)
}


########################################################################

writeSP <- function(rpDF, gpsList, parameterList = NULL, resDir = getwd(),
                    outputDirName = "sigPathway_results",
                    topIndexFileName = "TopPathwaysTable.html")  {
  if(!is.data.frame(rpDF) && !all(c("IndexG", "Pathway") %in% colnames(rpDF)))
    stop("'rpDF' must be a data frame from the output of rankPathways() or rankPathways.NGSk()\n")
  if(class(gpsList) != "list" && !all(unlist(lapply(gpsList, is.data.frame))))
    stop("'gpsList' must be a list (containing data frames) from getPathwayStatistics or getPathwayStatistics.NGSk()\n")
  if(nrow(rpDF) != length(gpsList))
    stop("Numbers of rows in 'rpDF' does not match length of 'gpsList'\n")

  resDir <- sub("/+$", "", resDir)
  outputDirName <- sub("/+$", "", outputDirName)
  outputDir <- file.path(resDir, outputDirName)
  
  if(!file.exists(resDir))
    stop("The directory specified in 'resDir' does not exist\n")
  if(file.exists(outputDir))
    stop("The directory name specified in 'outputDirName' already exists\n")
  if(!dir.create(outputDir))
    stop(paste("Unable to create the output directory", outputDir, "\n"))
  
  digits.rpDF <- numeric(ncol(rpDF))

  digits.rpDF[grep("Stat", colnames(rpDF))] <- 2
  digits.rpDF[grep("Rank", colnames(rpDF))] <- 1
  digits.rpDF[grep("q-value", colnames(rpDF))] <- 4

  rpDF2 <- rpDF
  rpDF2[,"Pathway"] <- 
    paste("<a href=\"pathways/pathway_", rpDF[,"IndexG"], ".html\">",
          rpDF[,"Pathway"], "</a>", sep = "")
  rpHTML <- c("<p><big><big>List of Top Pathways</big></big>",
              .convertDF2HTML(rpDF2, digits.rpDF))

  ## If provided, make table of parameters and provide link
  if(!is.null(parameterList))  {
    rpHTML <- c(rpHTML, "<a href=\"pathways/parameters.html\">Link to Parameters Used in the Analysis</a>")
    parameterList$phenotype <- paste(parameterList$phenotype, collapse = " ")
    parameterV <- unlist(parameterList)
    parameterDF <- data.frame(I(names(parameterV)), I(parameterV))
    rownames(parameterDF) <- 1:nrow(parameterDF)
    colnames(parameterDF) <- c("Parameter", "Value")
    parameterHTML <- c("<p><big><big>List of Parameters</big></big></p>",
                       .convertDF2HTML(parameterDF, NULL))
  }
    
  writeLines(rpHTML, file.path(outputDir, topIndexFileName))

  
  ## Make HTML tables for each pathway on the list of top pathways
  dir.create(file.path(outputDir, "pathways"))
  if(!is.null(parameterList))  {
    writeLines(parameterHTML,
               file.path(outputDir, "pathways", "parameters.html"))
  }
  
  tempDF <- gpsList[[1]]
  digits.gps <- numeric(ncol(tempDF))
  digits.gps[grep("^Mean", colnames(tempDF))] <- 1
  digits.gps[grep("^StDev", colnames(tempDF))] <- 1
  digits.gps[grep("^[TF][-]Statistic$", colnames(tempDF))] <- 3
  digits.gps[grep("^Statistic$", colnames(tempDF))] <- 3
  digits.gps[grep("^PearsonCorCoef$", colnames(tempDF))] <- 3
  digits.gps[grep("^FisherZ$", colnames(tempDF))] <- 3
  digits.gps[grep("^p-value$", colnames(tempDF))] <- 4

  for(i in 1:nrow(rpDF))  {
    temp.text <-
      c("<a href=\"../TopPathwaysTable.html\">Back to Table of Top Pathways</a>",
        paste("<p><big><big>", rpDF[i,"Pathway"], "</big></big></p>", sep = "")
        )
    tempDF <- gpsList[[i]]
    rownames(tempDF) <- 1:nrow(tempDF)
    temp <- .convertDF2HTML(tempDF, digits.gps)
    writeLines(c(temp.text, temp),
               file.path(outputDir, "pathways",
                         paste("pathway_", rpDF[i,"IndexG"], ".html", sep = "")
                         )
               )
  }

  cat("The results have been saved to the following directory:\n")
  cat(paste(outputDir, "\n"))
}

########################################################################

writeSigPathway <- function(spList, resDir = getwd(),
                            outputDirName = "sigPathway_results",
                            topIndexFileName = "TopPathwaysTable.html")  {
  if(class(spList) != "list")
    stop("'spList' needs to be a list\n")
  
  if(!all(names(spList) %in%
          c("gsList", "list.NTk", "list.NEk", "df.pathways",
            "list.gPS", "parameters"))
     )
    stop("'spList' is missing elements that are normally present in the list generated from runSigPathway()\n")

  writeSP(spList$df.pathways, spList$list.gPS, spList$parameters,
          resDir, outputDirName, topIndexFileName)
}
