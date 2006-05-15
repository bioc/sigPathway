\name{calculate.NGSk}
\alias{calculate.NGSk}
\title{Calculate NGSk (NTk-like) statistics with gene label permutation}
\description{
  Calculates the NGSk (NTk-like) statistics with gene label permutation
  and the corresponding p-values and q-values for each selected pathway.
}
\usage{
calculate.NGSk(statV, gsList, nsim = 1000, verbose = FALSE,
               alwaysUseRandomPerm = FALSE)
}
\arguments{
  \item{statV}{a numeric vector of test statistic (not p-values) for
    each individual probe/gene}
  \item{gsList}{a list containing three vectors from the output of
    the \code{selectGeneSets} function}
  \item{nsim}{an integer indicating the number of permutations to use}
  \item{verbose}{a boolean to indicate whether to print debugging messages
    to the R console}
  \item{alwaysUseRandomPerm}{a boolean to indicate whether the algorithm
    can use complete permutations for cases where \code{nsim} is greater
    than the total number of unique permutations possible with the
    \code{phenotype} vector}
}
\details{
  This function is a generalized version of NTk calculations;
  \code{calculate.NTk} calls this function internally.  To use this
  function, the user must specify a vector of test statistics (e.g.,
  t-statistic, Wilcoxon).  Pathways from this function can be ranked
  with \code{rankPathways.NGSk} or with \code{rankPathways} when
  combined with results from another pathway analysis algorithm (e.g.,
  \code{calculate.NEk}).
}
\value{
  A list containing
  \item{ngs}{number of gene sets}
  \item{nsim}{number of permutations performed}
  \item{t.set}{a numeric vector of Tk/Ek statistics}
  \item{t.set.new}{a numeric vector of NTk/NEk statistics}
  \item{p.null}{the proportion of nulls}
  \item{p.value}{a numeric vector of p-values}
  \item{q.value}{a numeric vector of q-values}
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

## Prepare the pathways to analyze
probeID <- rownames(tab)
gsList <- selectGeneSets(G, probeID, 20, 500)

nsim <- 1000
ngroups <- 2
verbose <- TRUE
weightType <- "constant"
methodName <- "NGSk"
npath <- 25
allpathways <- FALSE
annotpkg <- "hgu133a"

statV <- calcTStatFast(tab, phenotype, ngroups)$tstat
res.NGSk <- calculate.NGSk(statV, gsList, nsim, verbose)

## Summarize top pathways from NGSk
res.pathways.NGSk <-
  rankPathways.NGSk(res.NGSk, G, gsList, methodName, npath)
print(res.pathways.NGSk)

## Get more information about the probe sets' means and other statistics
## for the top pathway in res.pathways.NGSk
gpsList <-
  getPathwayStatistics.NGSk(statV, probeID, G, res.pathways.NGSk$IndexG,
                            FALSE, annotpkg)
print(gpsList[[1]])

## Write table of top-ranked pathways and their associated probe sets to
## HTML files
parameterList <-
  list(nprobes = nrow(tab), nsamples = ncol(tab),
       phenotype = phenotype, ngroups = ngroups,
       minNPS = 20, maxNPS = 500, ngs = res.NGSk$ngs,
       nsim.NGSk = res.NGSk$nsim,
       annotpkg = annotpkg, npath = npath, allpathways = allpathways)

writeSP(res.pathways.NGSk, gpsList, parameterList, tempdir(),
        "sigPathway_cNGSk", "TopPathwaysTable.html")

}
\keyword{array}
\keyword{htest}
