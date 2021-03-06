% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lolaUtils.R
\name{lolaBarPlot}
\alias{lolaBarPlot}
\title{lolaBarPlot}
\usage{
lolaBarPlot(
  lolaDb,
  lolaRes,
  scoreCol = "pValueLog",
  orderCol = scoreCol,
  signifCol = "qValue",
  includedCollections = c(),
  pvalCut = 0.01,
  maxTerms = 50,
  colorpanel = sample(rainbow(maxTerms, v = 0.5)),
  groupByCollection = TRUE,
  orderDecreasing = NULL
)
}
\arguments{
\item{lolaDb}{LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}}

\item{lolaRes}{LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package}

\item{scoreCol}{column name in \code{lolaRes} to be plotted}

\item{orderCol}{column name in \code{lolaRes} which is used for sorting the results}

\item{signifCol}{column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}}

\item{includedCollections}{vector of collection names to be included in the plot. If empty (default), all collections are used}

\item{pvalCut}{p-value cutoff to be employed for filtering the results}

\item{maxTerms}{maximum number of items to be included in the plot}

\item{colorpanel}{colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
vector indicates that black will be used for all bars.}

\item{groupByCollection}{facet the plot by collection}

\item{orderDecreasing}{flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
to increasing). \code{NULL} (default) for automatic determination.}
}
\value{
ggplot object containing the plot
}
\description{
plot a barplot of LOLA enrichment results
}
\examples{
\donttest{
library(RnBeads.hg19)
data(small.example.object)
logger.start(fname=NA)
# compute differential methylation
dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
# download LOLA DB
lolaDest <- tempfile()
dir.create(lolaDest)
lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
# perform enrichment analysis
res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
# select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
# in the example dataset
lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
# plot
lolaBarPlot(res$lolaDb, lolaRes, scoreCol="oddsRatio", orderCol="maxRnk", pvalCut=0.05)
}
}
\author{
Fabian Mueller
}
