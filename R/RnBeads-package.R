#' Analysis of genome-scale DNA methylation data with RnBeads
#'
#' RnBeads facilitates comprehensive analysis of various types of DNA methylation data at the genome scale. It extends
#' previous approaches for such analysis by high throughput capabilities, as well as presenting results in a
#' comprehensive, highly interpretable fashion. 
#'
#' The complete analysis can be performed by calling the function \code{\link{rnb.run.analysis}}. 
#' 
#' @references Yassen Assenov*, Fabian Mueller*, Pavlo Lutsik*, Joern Walter, Thomas Lengauer and Christoph Bock (2014) Compehensive Analysis of DNA Methylation Data with RnBeads, Nature Methods, 11(11):1138-1140.
#' @import methods MASS cluster RColorBrewer fields ggplot2 S4Vectors IRanges GenomicRanges methylumi ff gridExtra limma
#' @importFrom BiocGenerics annotation
#' @importFrom BiocGenerics annotation<-
#' @importFrom illuminaio readIDAT
#' @importFrom gplots colorpanel
#' @importFrom gplots heatmap.2
#' @importFrom plyr rbind.fill
#' @importFrom matrixStats rowVars
#' @importFrom matrixStats colVars
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats colMins
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
#' @importFrom matrixStats rowSds
#' @importFrom matrixStats colSds
#' @importFrom matrixStats rowQuantiles
#' @importFrom matrixStats colQuantiles
#' @importFrom matrixStats rowMedians
#' @importFrom matrixStats colMedians
#' @docType package
#' @name RnBeads
NULL

#' RnBeads option values and restrictions
#'
#' The values of options in RnBeads are stored in dedicated R objects accompanying the package. These objects are named
#' \code{infos}, \code{accepted}, \code{current} and \code{previous}. They should not be loaded or otherwise operated on
#' by users. Please refer to the documentation of \code{\link{rnb.options}} for accessing and modifying option values in
#' \pkg{RnBeads}.
#'
#' @docType data
#' @keywords datasets
#' @name accepted
#' @aliases current infos previous
#' @format \code{infos} is a \code{data.frame} containing information about all options in \pkg{RnBeads}. Row names in
#'         this table are the option names; the column names are \code{"Type"}, \code{"Named"}, \code{"Null"},
#'         \code{"Max"}, \code{"Min"}, \code{"MaxInclusive"} and \code{"MinInclusive"}.
#'         \code{accepted} is a \code{list} containing the sets of accepted values for some of the options.
#'         \code{current} is a \code{list} with current values for all options.
#'         \code{previous} is a \code{list} with previous values for the affected options; this list is only temporarily
#'         used while setting option values through \code{\link{rnb.options}} or \code{\link{rnb.xml2options}}.
#' @author Yassen Assenov
NULL
