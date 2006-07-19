\name{calcTStatFast}
\alias{calcTStatFast}
\title{Compute T-Statistics and Corresponding P-Values}
\description{
  Computes t-statistics and corresponding p-values.
}
\usage{
calcTStatFast(tab, phenotype, ngroups = 2)
}
\arguments{
  \item{tab}{a numeric matrix of expression values, with the rows and
    columns representing probe sets and sample arrays, respectively}
  \item{phenotype}{a numeric (or character if \code{ngroups} >= 2)
    vector indicating the phenotype}
  \item{ngroups}{an integer indicating the number of groups in the expression
    matrix}
}
\details{
  If there are two groups in the matrix, it is recommended to use 0 and
  1 to denote which sample columns belong to which group.  If the
  phenotype is a character vector, then the phenotype ranked first in
  the alphabet is considered as 0.
  
  If \code{ngroups} = 2, the t-test done here is equivalent to a unpaired
  two-sample t-test, assuming unequal variances.  Please note that as of
  version 1.1.6, the sign of the t-statistic is positive when the mean
  of group 1 is greater than the mean of group 0.
  
  If there is only one group in the matrix (e.g., Alzheimer's
  data set as reanalyzed in Tian et al. (2005)), then the phenotype vector
  should consist of continuous values.  In this case, the association
  between phenotype and expression values is first calculated as Pearson
  correlation coefficients, transformed to Fisher's z, and then rescaled
  so that its variance is 1:

  z = 0.5*log((1+rho)/(1-rho))*sqrt(n-3), where n is the number of
  phenotypes.

  If \code{ngroups} > 2, the f-statistics (from 1-way ANOVA) are
  calculated.  The user will need to check that the data have similar
  variances among the groups.
  
}
\value{
  \item{pval}{A vector of unadjusted p-values}
  \item{tstat}{A vector of t-statistics (\code{ngroups} = 2) or rescaled
    Fisher's z (\code{ngroups} = 1)}
  \item{rho}{(Also returned when \code{ngroups} = 1) A vector of Pearson
    correlation coefficients}
}
\author{Weil Lai}
\examples{
## Load inflammatory myopathy data set
data(MuscleExample)
statList <- calcTStatFast(tab, phenotype, ngroups = 2)

## Generate histogram of p-values
hist(statList$pval, xlab = "Unadjusted p-values", ylab = "Frequency")
}
\keyword{array}
\keyword{htest}
