\name{writeSigPathway}
\alias{writeSigPathway}
\alias{writeSP}
\title{Write results of pathway analysis to HTML format}
\description{
  Writes the table of top-ranked pathways and their associated probe set
  to HTML files.
}
\usage{
writeSigPathway(spList, resDir = getwd(),
                outputDirName = "sigPathway_results",
                topIndexFileName = "TopPathwaysTable.html")
writeSP(rpDF, gpsList, parameterList = NULL, resDir = getwd(),
        outputDirName = "sigPathway_results",
        topIndexFileName = "TopPathwaysTable.html")
}
\arguments{
  \item{spList}{a list containing the output from the
    \code{runSigPathway} function}
  \item{rpDF}{a data frame of top-ranked pathways from
    \code{rankPathways} or \code{rankPathways.NGSk}}
  \item{gpsList}{a list containing data frames of probes represented in
    gene sets from \code{getPathwayStatistics} or
    \code{getPathwayStatistics.NGSk}}
  \item{parameterList}{a list containing the values of parameters used
    in the analysis}
  \item{resDir}{a character string specifying the file directory to
    write the results}
  \item{outputDirName}{a character string specifying the folder to write
    the results within \code{resDir}}
  \item{topIndexFileName}{a character string specifying the name for the
    HTML file containing the table of top-ranked pathways}
}
\details{
  These functions export the results of the pathway analysis (e.g.,
  \code{runSigPathway}) to several HTML files.  The user can then
  quickly browse through the files for genes of interest within the
  top-ranked genes.
}
\value{
  None returned
}
\note{This function only uses the output of \code{runSigPathway} to
  generate the HTML files.  Please see the help page of
  \code{runSigPathway} for example usage.  The \code{writeSP} function
  should be used for those who have taken calculated the pathway
  statistics separately as shown in the help file of \code{calculate.NTk},
  \code{calculate.NEk}, and \code{calculate.NGSk}
}
\author{Weil Lai}
\keyword{array}
\keyword{htest}
