\name{calculate.GSEA}
\alias{calculate.GSEA}
\title{Calculate 2-sided statistics based on the GSEA algorithm}
\description{
  Calculates the 2-sided statistics based on the GSEA algorithm.
}
\usage{
calculate.GSEA(tab, phenotype, gsList, nsim = 1000,
               verbose = FALSE, alwaysUseRandomPerm = FALSE)
}
\arguments{
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric or character vector indicating the phenotype}
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
  This function assumes 2 distinct types of phenotypes in the data.
  It calculates a variant of the GSEA statistics (Mootha et
  al.) with the following modifications:
  (a) GSEA was changed from a 1-sided to a 2-sided approach.
  (b) The 2-group t-statistics is used as the difference metric.

  The function also normalizes the GSEA statistic and calculates the
  corresponding q-values for each gene set as described in Tian
  et al. (2005)  The function's output can be used for further analysis
  in other functions such as \code{rankPathways.NGSk} or
  \code{getPathwayStatistics.NGSk}.
}
\value{
  A list containing
  \item{ngs}{number of gene sets}
  \item{nsim}{number of permutations performed}
  \item{t.set}{a numeric vector of Tk statistics}
  \item{t.set.new}{a numeric vector of NTk statistics}
  \item{p.null}{the proportion of nulls}
  \item{p.value}{a numeric vector of p-values}
  \item{q.value}{a numeric vector of q-values}
}
\references{
  Mootha V.K., Lindgren C.M., Eriksson K.F., Subramanian A., Sihag S.,
  Lehar J., Puigserver P., Carlsson E., Ridderstrale M., Laurila E.,
  Houstis N., Daily M.J., Patterson N., Mesirov J.P., Golud T.R., Tamayo
  P., Spiegelman B., Lander E.S., Hirshhorn J.N., Altshuler D., Groop
  L.C.  (2003)  PGC-1alpha-responsive genes involved in oxidative
  phosphorylation are coordinately downregulated in human diabetes.
  \emph{Nature Genetics}, \bold{34}, 267-73.
  
  Tian L., Greenberg S.A., Kong S.W., Altschuler J., Kohane I.S., Park
  P.J. (2005)  Discovering statistically significant pathways in
  expression profiling studies.  \emph{Proceedings of the National
    Academy of Sciences of the USA}, \bold{102}, 13544-9.

  \url{http://www.pnas.org/cgi/doi/10.1073/pnas.0506577102}
}
\author{Lu Tian, Peter Park, and Weil Lai}
\keyword{array}
\keyword{htest}
