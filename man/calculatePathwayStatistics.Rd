\name{calculatePathwayStatistics}
\alias{calculate.NTk}
\alias{calculate.NEk}
\title{Calculate the NTk and NEk statistics}
\description{
  Calculates the NTk and NEk statistics and the corresponding p-values and
  q-values for each selected pathway.
}
\usage{
calculate.NTk(tab, phenotype, gsList, nsim = 1000,
              ngroups = 2, verbose = FALSE, alwaysUseRandomPerm = FALSE)
calculate.NEk(tab, phenotype, gsList, nsim = 1000,
              weightType = c("constant", "variable"),
              ngroups = 2, verbose = FALSE, alwaysUseRandomPerm = FALSE)
}
\arguments{
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{gsList}{a list containing three vectors from the output of
    the \code{selectGeneSets} function}
  \item{nsim}{an integer indicating the number of permutations to use}
  \item{weightType}{a character string specifying the type of weight to
    use when calculating NEk statistics}
  \item{ngroups}{an integer indicating the number of groups in the
    matrix}
  \item{verbose}{a boolean to indicate whether to print debugging messages
    to the R console}
  \item{alwaysUseRandomPerm}{a boolean to indicate whether the algorithm
    can use complete permutations for cases where \code{nsim} is greater
    than the total number of unique permutations possible with the
    \code{phenotype} vector}
}
\details{
  These functions calculate the NTk and NEk statistics and the
  corresponding p-values and q-values for each selected pathway.  The output
  of both functions should be together to rank top pathways with
  the \code{rankPathways} function.
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
\author{Lu Tian,Peter Park, and Weil Lai}
\examples{
## Load in filtered, expression data
data(MuscleExample)

## Prepare the pathways to analyze
probeID <- rownames(tab)
gsList <- selectGeneSets(G, probeID, 20, 500)

## Calculate NTk and weighted NEk for each gene set
## * Use a higher nsim (e.g., 2500) value for more reproducible results
nsim <- 1000
ngroups <- 2
verbose <- TRUE
weightType <- "constant"
methodNames <- c("NTk", "NEk")
npath <- 25
allpathways <- FALSE
annotpkg <- "hgu133a"

res.NTk <- calculate.NTk(tab, phenotype, gsList, nsim, ngroups, verbose)
res.NEk <- calculate.NEk(tab, phenotype, gsList, nsim, weightType,
                         ngroups, verbose)

## Summarize results
res.pathways <- rankPathways(res.NTk, res.NEk, G, tab, phenotype,
                             gsList, ngroups, methodNames, npath, allpathways)
print(res.pathways)

## Get more information about the probe sets' means and other statistics
## for the top pathway in res.pathways
statList <- calcTStatFast(tab, phenotype, ngroups)
gpsList <-
  getPathwayStatistics(tab, phenotype, G, res.pathways$IndexG,
                       ngroups, statList, FALSE, annotpkg)
print(gpsList[[1]])

## Write table of top-ranked pathways and their associated probe sets to
## HTML files
parameterList <-
  list(nprobes = nrow(tab), nsamples = ncol(tab),
       phenotype = phenotype, ngroups = ngroups,
       minNPS = 20, maxNPS = 500, ngs = res.NTk$ngs,
       nsim.NTk = res.NTk$nsim, nsim.NEk = res.NEk$nsim,
       weightType = weightType,
       annotpkg = annotpkg, npath = npath, allpathways = allpathways)

writeSP(res.pathways, gpsList, parameterList, tempdir(), "sigPathway_cPS",
        "TopPathwaysTable.html")

}
\keyword{array}
\keyword{htest}
