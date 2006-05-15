\name{runSigPathway}
\alias{runSigPathway}
\title{Perform pathway analysis}
\description{
  Performs pathway analysis
}
\usage{
runSigPathway(G, minNPS = 20, maxNPS = 500,
              tab, phenotype, nsim = 1000,
              weightType = c("constant", "variable"), ngroups = 2,
              npath = 25, verbose = FALSE, allpathways = FALSE,
              annotpkg = NULL, alwaysUseRandomPerm = FALSE)
}
\arguments{
  \item{G}{a list containing the source, title, and probe sets
    associated with each curated pathway}
  \item{minNPS}{an integer specifying the minimum number of probe sets
    in \code{tab} that should be in a gene set}
  \item{maxNPS}{an integer specifying the maximum number of probe sets
    in \code{tab} that should be in a gene set}
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{nsim}{an integer indicating the number of permutations to use}
  \item{weightType}{a character string specifying the type of weight to
    use when calculating NEk statistics}
  \item{ngroups}{an integer indicating the number of groups in the
    matrix}
  \item{npath}{an integer indicating the number of top gene sets to
    consider from each statistic when ranking the top pathways}
  \item{verbose}{a boolean to indicate whether to print debugging messages
    to the R console}
  \item{allpathways}{a boolean to indicate whether to include the top
    npath pathways from each statistic or just consider the top npath
    pathways (sorted by the sum of ranks of both statistics) when
    generating the summary table}
  \item{annotpkg}{a character vector specifying the name of the
    BioConductor annotation package to use to fetch accession numbers,
    Entrez Gene IDs, gene name, and gene symbols}
  \item{alwaysUseRandomPerm}{a boolean to indicate whether the algorithm
    can use complete permutations for cases where \code{nsim} is greater
    than the total number of unique permutations possible with the
    \code{phenotype} vector}
}
\details{
  \code{runSigPathway} is a wrapper function that
  
  (1) Selects the gene sets to analyze using \code{selectGeneSets}
  
  (2) Calculates NTk and NEk statistics using \code{calculate.NTk} and
  \code{calculate.NEK}
  
  (3) Ranks the top \code{npath} pathways from each statistic using
  \code{rankPathways}
  
  (4) Summarizes the means, standard deviation, and individual
  statistics of each probe set in each of the above pathways using
  \code{getPathwayStatistics}
}
\value{
  A list containing
  \item{gsList}{a list containing three vectors from the output of
    the \code{selectGeneSets} function}
  \item{list.NTk}{a list from the output of calculate.NTk}
  \item{list.NEk}{a list from the output of calculate.NEk}
  \item{df.pathways}{a data frame from \code{rankPathways} which
    contains the top pathways' indices in \code{G}, gene set category,
    pathway title, set size, NTk statistics, NEk statistics, the
    corresponding q-values, and the ranks. 
  }
  \item{list.gPS}{a list from \code{getPathwayStatistics} containing
    \code{nrow(df.pathways)} data frames corresponding to the pathways
    listed in \code{df.pathways}.  Each data frame contains the
    name, mean, standard deviation, the test statistic (e.g., t-test),
    and the corresponding unadjusted p-value.  If \code{ngroups} = 1, the
    Pearson correlation coefficient is also returned.  If a valid
    \code{annotpkg} is specified, the probes' accession numbers, Entrez
    Gene IDs, gene name, and gene symbols are also returned.}
  \item{parameters}{a list of parameters (e.g., \code{nsim}) used in the
    analysis}
}
\references{
  Tian L., Greenberg S.A., Kong S.W., Altschuler J., Kohane I.S., Park
  P.J. (2005)  Discovering statistically significant pathways in
  expression profiling studies.  \emph{Proceedings of the National
    Academy of Sciences of the USA}, \bold{102}, 13544-9.

  \url{http://www.pnas.org/cgi/doi/10.1073/pnas.0506577102}
}
\author{Lu Tian, Peter Park, and Weil Lai}
\examples{
## Load in filtered, expression data
data(MuscleExample)

## Prepare the pathways to analyze and run analysis with 1 wrapper function

nsim <- 1000
ngroups <- 2
verbose <- TRUE
weightType <- "constant"
npath <- 25
allpathways <- FALSE
annotpkg <- "hgu133a"

res.muscle <- runSigPathway(G, 20, 500, tab, phenotype, nsim,
                            weightType, ngroups, npath, verbose,
                            allpathways, annotpkg)

## Summarize results
print(res.muscle$df.pathways)

## Get more information about the probe sets' means and other statistics
## for the top pathway in res.pathways
print(res.muscle$list.gPS[[1]])

## Write table of top-ranked pathways and their associated probe sets to
## HTML files
writeSigPathway(res.muscle, tempdir(), "sigPathway_rSP",
                "TopPathwaysTable.html")

}
\keyword{array}
\keyword{htest}
