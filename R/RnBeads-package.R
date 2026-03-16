#' Analysis of genome-scale DNA methylation data with RnBeads
#'
#' RnBeads facilitates comprehensive analysis of various types of DNA methylation data at the genome scale. It extends
#' previous approaches for such analysis by high throughput capabilities, as well as presenting results in a
#' comprehensive, highly interpretable fashion. 
#'
#' The complete analysis can be performed by calling the function \code{\link{rnb.run.analysis}}. 
#' 
#' @references Yassen Assenov*, Fabian Mueller*, Pavlo Lutsik*, Joern Walter, Thomas Lengauer and Christoph Bock (2014) Compehensive Analysis of DNA Methylation Data with RnBeads, Nature Methods, 11(11):1138-1140.
#' @import methods MASS cluster fields ggplot2 S4Vectors IRanges GenomicRanges methylumi ff limma
#' @importFrom BiocGenerics annotation
#' @importFrom BiocGenerics annotation<-
#' @importFrom grDevices colorRampPalette densCols dev.control dev.off dev2bitmap pdf png rainbow rgb xy.coords
#' @importFrom graphics abline close.screen layout legend lines mtext par plot plot.new polygon screen split.screen
#' @importFrom gridExtra arrangeGrob
#' @importFrom gplots colorpanel heatmap.2
#' @importFrom illuminaio readIDAT
#' @importFrom matrixStats colMaxs colMedians colMins colQuantiles colSds colVars
#' @importFrom matrixStats rowMaxs rowMedians rowMins rowQuantiles rowSds rowVars
#' @importFrom plyr rbind.fill
#' @importFrom stats as.dendrogram as.dist as.formula bartlett.test coef cutree dbeta density dexp dist dnorm ecdf
#' @importFrom stats fisher.test hclust knots kruskal.test lm model.matrix optim p.adjust pbeta pchisq pf pnorm
#' @importFrom stats prcomp predict pt qbeta residuals rnorm rt t.test vcov wilcox.test
#' @importFrom utils browseURL capture.output combn data download.file installed.packages memory.size
#' @importFrom utils read.csv read.delim read.table untar unzip write.table zip
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

#' LUMP Support
#'
#' The sites used by the LUMP algorithm for estimating immune cell content are stored in an object named
#' \code{lump.hg19}. This object should not be loaded or otherwise operated on by users. Please refer to the
#' documentation of \code{\link{rnb.execute.lump}} for information on the algorithm and its implementation in
#' \pkg{RnBeads}.
#'
#' @docType data
#' @keywords datasets
#' @name lump.hg19
#' @format \code{lump.*} is a \code{list} of non-empty \code{integer} matrices, one per supported platform. Every
#'         \code{matrix} contains exactly two columns, denoting chromosome index and chromosome-based index,
#'         respectively. These indices refer to positions within the probe/site annotation table employed by
#'         \pkg{RnBeads} for the corresponding platform.
#' @author Yassen Assenov
NULL

#' LUMP Support (hg38)
#'
#' Those are the same sites as reported in \code{lump.hg19}, but lifted to `hg38` with UCSC's liftOver functionality. This only
#' applies for the CpG-wise sites; i.e. those used for sequencing data sets, since `hg38` is not supported for array-based data
#' sets.
#'
#' @docType data
#' @keywords datasets
#' @name lump.hg38
#' @format \code{lump.*} is a \code{list} of non-empty \code{integer} matrices, one per supported platform. Here, only `CpG` is 
#'          available for BS datasets.
#' @author Michael Scherer
NULL

#' XenoMix EPICv1 (hg38)
#'
#' This is the manifest file for EPICv1, which is required for detecting inter-species probes for estimating mouse contamination
#' with XenoMix
#'
#' @docType data
#' @keywords datasets
#' @name epic_v1_manifest
#' @format \code{data.frame} provided by Illumina capturing information about the probe design.
#' @author Nima Esmaeelpour
NULL

#' XenoMix EPICv2 (hg38)
#'
#' This is the manifest file for EPICv2, which is required for detecting inter-species probes for estimating mouse contamination
#' with XenoMix
#'
#' @docType data
#' @keywords datasets
#' @name epic_v2_manifest
#' @format \code{data.frame} provided by Illumina capturing information about the probe design.
#' @author Nima Esmaeelpour
NULL

#' XenoMix interspeciesprobes (v1)
#'
#' These are the inter-species probes for EPICv1
#'
#' @docType data
#' @keywords datasets
#' @name interspecies_probes_v1
#' @format \code{list} of CGIs for the inter-species probes.
#' @author Nima Esmaeelpour
NULL

#' XenoMix interspeciesprobes (v2)
#'
#' These are the inter-species probes for EPICv2
#'
#' @docType data
#' @keywords datasets
#' @name interspecies_probes_v2
#' @format \code{list} of CGIs for the inter-species probes.
#' @author Nima Esmaeelpour
NULL