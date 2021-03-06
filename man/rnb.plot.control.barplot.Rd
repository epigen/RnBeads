% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/controlPlots.R
\name{rnb.plot.control.barplot}
\alias{rnb.plot.control.barplot}
\title{rnb.plot.control.barplot}
\usage{
rnb.plot.control.barplot(
  rnb.set,
  probe,
  sample.subset = 1:length(samples(rnb.set)),
  writeToFile = FALSE,
  numeric.names = FALSE,
  name.prefix = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{rnb.set}{\code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object with valid
quality control information}

\item{probe}{exact id of the control probe consisting of the control probe type (see \code{\link{rnb.plot.control.boxplot}})}

\item{sample.subset}{an integer vector specifying the subset of samples for which the plotting should be performed}

\item{writeToFile}{flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}}

\item{numeric.names}{if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits}

\item{name.prefix}{in case \code{writeToFile} is \code{TRUE}, a \code{character} singleton specifying a prefix to the variable part of the image file names}

\item{verbose}{if \code{TRUE} additional diagnostic output is generated}

\item{...}{other arguments to \code{\link{createReportPlot}}}
}
\value{
plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class
			\code{\link{ggplot}} otherwise.
}
\description{
Per-sample bar plots of Illumina HumanMethylation control probes
}
\examples{
\donttest{
library(RnBeads.hg19)
data(small.example.object)
control.meta.data <- rnb.get.annotation("controls450")
ctrl.probe<-paste0(unique(control.meta.data[["Target"]])[4], ".3")
print(ctrl.probe) # EXTENSION.3
rnb.plot.control.barplot(rnb.set.example, ctrl.probe)
}

}
\author{
Pavlo Lutsik
}
