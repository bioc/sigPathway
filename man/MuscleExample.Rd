\name{MuscleExample}
\docType{data}
\alias{MuscleExample}
\alias{tab}
\alias{phenotype}
\alias{G}
\title{Subset of Inflammatory Myopathy Dataset to Demonstrate
  sigPathway}
\description{
  \code{MuscleExample} is an R workspace containing the following objects:
  (1) \code{tab}:        a matrix of 5000 rows and 15 columns
  (2) \code{phenotype}:  a indicator vector which denotes which columns in
  \code{tab} are arrays from normal (NORM) and inclusion body myositis
  (IBM)
  (3) \code{G}:          a list containing the source, title, and the
  probe set IDs associated with 626 pathways

  The full inflammatory myopathway dataset (which includes all probe
  sets and samples, including more NORM, IBM, and dermatomyositis
  arrays) and a more comprehensive pathway annotation list for the
  HG-U133A and other selected array platforms are available at
  \url{http://www.chip.org/~ppark/PNAS05/}

  Although the objects contained in \code{MuscleExample} are subsets of
  the full dataset, the results obtained from running pathway analysis
  with \code{MuscleExample} are comparable to those obtained using the
  full dataset.  This example dataset contains 8 IBM and 7 NORM arrays.
  The 5000 probe sets were selected by considering the variance of the
  expression values of each probe set among the 15 arrays. 
}

\usage{data(MuscleExample)}
\format{1 integer matrix, 1 numeric vector, and 1 list}
\source{\url{http://www.chip.org/~ppark/PNAS05/}}
\references{
  Tian L., Greenberg S.A., Kong S.W., Altschuler J., Kohane I.S., Park
  P.J. (2005)  Discovering statistically significant pathways in
  expression profiling studies.  \emph{Proceedings of the National
    Academy of Sciences of the USA}, \bold{102}, 13544-9.

  \url{http://www.pnas.org/cgi/doi/10.1073/pnas.0506577102}
}
\keyword{datasets}
