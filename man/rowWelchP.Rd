% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differentialMethylation.R
\name{rowWelchP}
\alias{rowWelchP}
\title{rowWelchP}
\usage{
rowWelchP(
  X,
  inds.g1,
  inds.g2 = -inds.g1,
  na.rm = FALSE,
  alternative = "two.sided"
)
}
\arguments{
\item{X}{Matrix on which the test is performed for every row}

\item{inds.g1}{column indices of group 1 members}

\item{inds.g2}{column indices of group 2 members}

\item{na.rm}{Should NAs be removed (logical)}

\item{alternative}{Testing alternative. Must be one of "two.sided" (default),"less","greater" or "all".
in case of "all" a data frome with corresping alternative variables is returned. 
Otherwise the result is a vector.}
}
\value{
vector (or data.frame if alternative=="all") of p-values resulting from the Welch's t-test
}
\description{
performs a two-sided Welch's t-test (unequal variances, unequal sample sizes) on each row of a matrix X with the indices inds.1 vs indices inds.g2 as group assignments.
}
\note{
Requires \code{matrixStats} package
}
\examples{
\donttest{
library(RnBeads.hg19)
data(small.example.object)
logger.start(fname=NA)
meth.mat <- meth(rnb.set.example)
sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
p.vals <- rowWelchP(meth.mat,sample.groups[[1]],sample.groups[[2]])
}
}
\author{
Fabian Mueller
}
