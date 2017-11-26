########################################################################################################################
## RnBSet-class.R
## created: 2012-04-06
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBSet class definition.
########################################################################################################################
## GLOBALS

RNBSET.SLOTNAMES<-c("meth.sites", "covg.sites")

##
## ---------------------------------------------------------------------------------------------------------------------
## CLASS DEFINITIONS
## ---------------------------------------------------------------------------------------------------------------------
#' @include bigFf.R
setOldClass(c("ff_matrix"))
setClassUnion("matrixOrff", c("matrix", "ff_matrix"))
setClassUnion("matrixOrffOrBigFfMat", c("matrix", "ff_matrix", "BigFfMat"))
setClassUnion("matrixOrffOrNULL", c("matrix", "ff_matrix", "NULL"))
setClassUnion("matrixOrffOrBigFfMatOrNULL", c("matrix", "ff_matrix", "BigFfMat", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))

#' RnBSet Class
#'
#' Basic class for storing DNA methylation and experimental quality information
#'
#' @details
#' It is a virtual class and objects of type \code{RnBSet} should not be instantiated. Instead, the child classes are
#' used: \code{\linkS4class{RnBeadRawSet}} and \code{\linkS4class{RnBeadSet}} for Infinium HumanMethylation and
#' \code{\linkS4class{RnBiseqSet}} for bisulfite sequencing data
#'
#' @section Slots:
#' \describe{
#'   \item{\code{pheno}}{Sample annotations (phenotypic and processing data) in the form of a \code{data.frame}.}
#'   \item{\code{sites}}{A \code{matrix} object storing the identifiers of the methylation sites for which the
#' 		methylation information is present}
#'   \item{\code{meth.sites}}{\code{matrix} of methylation values. Every row corresponds to a methylation site,
#' 		and every column - to a sample.}
#'   \item{\code{covg.sites}}{\code{matrix} of coverage values. Every row corresponds to a methylation site,
#' 		and every column - to a sample.}
#'   \item{\code{regions}}{\code{list} of all identifiers of methylation sites for which methylation information
#' 		is available.}
#'   \item{\code{meth.regions}}{\code{list} of methylation \code{matrix} objects, one per available region type. Every row in a
#' 		matrix corresponds to a methylation site, and every column - to a sample.}
#'   \item{\code{covg.regions}}{\code{list} of coverage \code{matrix} objects, one per available region type.
#' 		Every row corresponds to a region, and every column - to a sample.}
#' 	 \item{\code{status}}{\code{list} with meta-information about the object.}
#' 	 \item{\code{assembly}}{\code{character} vector of length one, specifying the genome assembly which the object is linked to, e.g. "hg19".}
#'   \item{\code{target}}{\code{character} vector of length one, specifying the feature class:
#' 		\code{"CpG"} for sequencing data, \code{"probes450"} and \code{"probes27"} for
#' 		HumanMethylation450 and HumanMethylation27 microarrays respectively.}
#'   \item{\code{inferred.covariates}}{\code{list} with covariate information.
#' 		Can contain elements \code{"sva"} and \code{"cell.types"}.}
#' 	 \item{\code{version}}{Package version in which the dataset was created.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link[=pheno,RnBSet-method]{pheno}}}{Gets the phenotypic and processing data of the dataset.}
#'   \item{\code{\link[=samples,RnBSet-method]{samples}}}{Gets the identifiers of all samples in the dataset.}
#'   \item{\code{\link[=summarized.regions,RnBSet-method]{summarized.regions}}}{Gets the genomic annotations for
#'   which methylation data is present.}
#' 	 \item{\code{\link[=meth,RnBSet-method]{meth}}}{Gets a \code{matrix} of methylation values in the dataset.}
#' 	 \item{\code{\link[=mval,RnBSet-method]{mval}}}{Gets a \code{matrix} of M values in the dataset.}
#'   \item{\code{\link[=covg,RnBSet-method]{covg}}}{Gets the \code{matrix} of coverage values of the dataset.}
#'   \item{\code{\link[=remove.sites,RnBSet-method]{remove.sites}}}{Removes sites from the dataset.}
#'   \item{\code{\link[=remove.samples,RnBSet-method]{remove.samples}}}{Removes samples from the dataset.}
#'   \item{\code{\link[=addPheno,RnBSet-method]{addPheno,RnBSet-method}}}{Add sample annotation to the dataset.}
#'   \item{\code{\link[BiocGenerics]{combine}}}{Combines two datasets.}
#'   \item{\code{\link{regionMapping,RnBSet-method}}}{Retrieve the sites mapping to a given region type}
#'   \item{\code{\link[=rnb.sample.summary.table,RnBSet-method]{rnb.sample.summary.table}}}{Creates a sample summary table from an RnBSet object.}
#' }
#'
#' @name RnBSet-class
#' @rdname RnBSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBSet
setClass("RnBSet",
		representation(pheno="data.frame",
				sites="matrix",
				meth.sites="matrixOrffOrBigFfMat",
				covg.sites="matrixOrffOrBigFfMatOrNULL",
				regions="list",
				meth.regions="list",
				covg.regions="listOrNULL",
				status="listOrNULL",
				assembly="character",
				target="characterOrNULL",
				inferred.covariates="list",
				version="characterOrNULL"),
		prototype(pheno=data.frame(),
				sites=matrix(nrow=0, ncol=0),
				meth.sites=matrix(nrow=0, ncol=0),
				covg.sites=NULL,
				regions=list(),
				meth.regions=list(),
				covg.regions=NULL,
				status=NULL,
				assembly="hg19",
				target=NULL,
				inferred.covariates=list(),
				version=as.character(packageVersion("RnBeads"))),
		contains = "VIRTUAL",
		package = "RnBeads")

## ---------------------------------------------------------------------------------------------------------------------
## DUMMY CONSTRUCTOR
## ---------------------------------------------------------------------------------------------------------------------
#
#setMethod("initialize", "RnBSet",
#		function(pheno=data.frame(),
#				sites=matrix(),
#				meth.sites=matrix(),
#				covg.sites=NULL,
#				regions=list(),
#				meth.regions=list(),
#				covg.regions=NULL,
#				status=NULL,
#				assembly="hg19",
#				target=NULL,
#				inferred.covariates=list()
#				){
#					.Object@pheno<-pheno
#					.Object@sites<-sites
#					.Object@meth.sites<-betas
#					.Object@covg.sites<-covg.sites
#
#					.Object@status<-status
#
#					.Object@target<-target
#
#
#				})

## ---------------------------------------------------------------------------------------------------------------------
## ACCESSORS
## ---------------------------------------------------------------------------------------------------------------------

if (!isGeneric("pheno")) setGeneric("pheno", function(object) standardGeneric("pheno"))

#' pheno-methods
#'
#' Extracts sample phenotype and/or processing information.
#'
#' @param object Dataset of interest.
#' @return Sample annotation information available for the dataset in the form of a \code{data.frame}.
#'
#' @rdname pheno-methods
#' @docType methods
#' @aliases pheno
#' @aliases pheno,RnBSet-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' pheno(rnb.set.example)
#' }
setMethod("pheno", signature(object="RnBSet"), function(object) object@pheno)

########################################################################################################################

if (!isGeneric("samples")) {
	setGeneric("samples", function(object) standardGeneric("samples"))
}

#' samples-methods
#'
#' Extracts sample identifiers
#'
#' @param object Dataset of interest.
#'
#' @details The column of the sample annotation table which contains identifiers is globally controlled via the
#'  \code{"identifiers.column"} option. In case the latter is \code{NULL} column names of the matrix returned
#' by the \code{meth} method are treated as sample identifiers. In case the latter are also missing, a \code{character}
#' vector with sample numbers is returned.
#'
#' @return \code{character} vector of sample identifiers.
#'
#' @rdname samples-methods
#' @docType methods
#' @aliases samples
#' @aliases samples,RnBSet-method
#' @aliases samples,RnBeadClustering-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' samples(rnb.set.example)
#' }
setMethod("samples", signature(object="RnBSet"),
	function(object) {
		pheno.table <- pheno(object)
		id.column <- rnb.getOption("identifiers.column")
		ids <- NULL
		if (!(is.null(pheno.table) || is.null(id.column))) {
			if (is.character(id.column)) {
				if (id.column %in% colnames(pheno.table)) {
					ids <- pheno.table[, id.column]
				}
			} else if (1L <= id.column && id.column <= ncol(pheno.table)) {
				ids <- pheno.table[, id.column]
			}

			if (is.null(ids) || any(is.na(ids)) || anyDuplicated(ids) != 0) {
				rnb.warning("The supplied identifiers column is not found or is not suitable")
			}
			ids <- as.character(ids)
		} else if (!is.null(colnames(object@meth.sites))) {
			ids <- colnames(object@meth.sites)
		} else {
			ids <- as.character(1:nrow(object@pheno))
		}
		ids
	}
)

########################################################################################################################

if(!isGeneric("sites")) setGeneric("sites",
			function(object) standardGeneric("sites"))

#' sites-methods
#'
#' Methylation sites object information for which is present in the \code{RnBSet} object.
#'
#' @param object Dataset of interest.
#'
#' @return A matrix of type \code{integer} describing the sites, information for which is
#' present in the \code{object}
#'
#' @rdname sites-methods
#' @docType methods
#' @aliases sites
#' @aliases sites,RnBSet-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' sites(rnb.set.example)
#' }
setMethod("sites", signature(object="RnBSet"),
		function(object){
			return(object@sites)
		})

if(!isGeneric("regions")) setGeneric("regions",
			function(object, ...) standardGeneric("regions"))
########################################################################################################################

#' regions-methods
#'
#' Methylation regions, information for which is present in the \code{RnBSet} object.
#'
#' @param object Dataset of interest.
#' @param type   Region type(s) of interest as a \code{character} vector. If this is set to \code{NULL}, all region
#'               types summarized in the object are returned.
#' @return Methylation site and region assignment. If \code{type} is singleton, a \code{matrix} is returned. The first
#' 		   column corresponds to the methylation context index. The second column is the index of the chromosome in
#'         the genome, and the third is the index of the region in the \code{GRanges} object of the region type
#'         annotation. When \code{length(type)>1}, a list of such matrices is returned for each element of \code{type}.
#' 		   If \code{type} is \code{NULL}, matrices for all summarized region types are returned.
#'
#' @note
#' Methylation context index is an integer number denoting the sequence context of the cytosine of interest. Index
#' \code{1} corresponds to \code{CpG}, the only supported index in bisulfite sequencing datasets.
#'
#' @rdname regions-methods
#' @docType methods
#' @aliases regions
#' @aliases regions,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' head(regions(rnb.set.example))
#' }
#' @seealso \code{\link[=summarized.regions,RnBSet-method]{summarized.regions}} for all summarized region types in a dataset;
#'   \code{\link{rnb.get.chromosomes}} listing all supported chromosomes for a given genome assembly
#' @author Pavlo Lutsik
#' @export
setMethod("regions", signature(object="RnBSet"),
		function(object, type=NULL){
			if(!(is.character(type)))
				stop("Invalid argument type")

			if(is.null(object@regions)){
				warning("No region information present, returning NULL")
				return(NULL)
			}
			if(!is.null(type)){
				if(!all(type %in% names(object@regions)))
					stop(sprintf("No information for type %s",type))
				if(length(type==1))
				return(object@regions[[type]]) else
				return(object@regions[type])
			}else{
				return((object@regions))
			}
		})

########################################################################################################################

if (!isGeneric("summarized.regions")) {
	setGeneric("summarized.regions", function(object) standardGeneric("summarized.regions"))
}

#' summarized.regions-methods
#'
#' Gets the genomic annotations for which methylation data is present in the \code{RnBSet} object.
#'
#' @param object Methylation dataset of interest.
#'
#' @return \code{character} vector listing all genomic annotations summarized in the given dataset. If the dataset
#'         contains methylation in sites only, an empty vector is returned.
#'
#' @seealso \code{\link[=summarize.regions,RnBSet-method]{summarize.regions}} for calculating region-wise methylation in a dataset;
#'          \code{\link{rnb.set.annotation}} for adding or replacing a region annotation table
#'
#' @rdname summarized.regions-methods
#' @docType methods
#' @aliases summarized.regions
#' @aliases summarized.regions,RnBSet-method
#' @author Yassen Assenov
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' summarized.regions(rnb.set.example)
#' }
setMethod("summarized.regions", signature(object = "RnBSet"),
	function(object) {
		result <- names(object@regions)
		if (is.null(result)) {
			result <- character()
		}
		result
	}
)

########################################################################################################################

## get.dataset.matrix
##
## Extracts a specific data matrix from the given methylation dataset and sets row names if necessary.
##
## @param object     Methylation dataset as an object of class inheriting \code{RnBSet}.
## @param type       Site type (e.g. \code{"sites"} or \code{"probes450"}) for site/probe matrix, or region name for
##                   the corresponding region-based matrix.
## @param row.names  Flag indicating if row names must be generated.
## @param mm.sites   Data matrix for the site level.
## @param mm.regions List of data matrices, one per supported region type.
## @param i          indices of sites/regions to be retrieved (index or logical). retrieves all if \code{NULL} (default).
## @param j          indices of samples to be retrieved (index or logical). retrieves all if \code{NULL} (default).
## @return Requested data matrix. Note that this might be \code{NULL}.
## @author Pavlo Lutsik
get.dataset.matrix <- function(object, type, row.names, mm.sites, mm.regions, i=NULL, j=NULL) {
	if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
		stop("invalid value for type")
	}
	if (!parameter.is.flag(row.names)) {
		stop("invalid value for row.names; expected TRUE or FALSE")
	}
	if (!is.element(class(i), c("NULL", "integer", "numeric", "logical"))) {
		stop("invalid value for i; expected NULL, index or logical")
	}
	if (!is.element(class(j), c("NULL", "integer", "numeric", "logical", "character"))) {
		stop("invalid value for j; expected NULL, index, character or logical")
	}
	if (is.character(j)){
		j <- match(j, samples(object))
		if (any(is.na(j))){
			stop("invalid sample names")
		}
	}
	if (type %in% c("sites", object@target)) {
		if (is.null(mm.sites)) {
			return(NULL)
		}
		if("ff" %in% class(mm.sites)){
			open(mm.sites)
		}
		if (is.null(i) && is.null(j)){
			result <- mm.sites[, , drop = FALSE]
		} else if(is.null(i)){
			result <- mm.sites[, j, drop = FALSE]
		} else if(is.null(j)){
			result <- mm.sites[i, , drop = FALSE]
		} else {
			result <- mm.sites[i, j, drop = FALSE]
		}
	} else if (!(type %in% names(object@regions))) {
		stop("unsupported region type")
	} else if (is.null(mm.regions[[type]])) {
		return(NULL)
	} else {
		if (is.null(i) && is.null(j)){
			result <- mm.regions[[type]][, , drop = FALSE]
		} else if(is.null(i)){
			result <- mm.regions[[type]][, j, drop = FALSE]
		} else if(is.null(j)){
			result <- mm.regions[[type]][i, , drop = FALSE]
		} else {
			result <- mm.regions[[type]][i, j, drop = FALSE]
		}
	}
	if (is.null(j)){
		colnames(result) <- samples(object)
	} else {
		colnames(result) <- samples(object)[j]
	}
	if (row.names) {
		if (is.null(i)){
			rownames(result) <- get.row.names(object, type)
		} else {
			rownames(result) <- get.row.names(object, type)[i]
		}
	} else {
		rownames(result) <- NULL
	}
	return(result)
}

########################################################################################################################

if(!isGeneric("mval")) setGeneric("mval", function(object, ...) standardGeneric("mval"))
#' mval-methods
#'
#' Extracts DNA methylation information (M values) for a specified set of genomic features.
#'
#' @param object 	dataset of interest.
#' @param type 		\code{character} singleton. If this is set to \code{"sites"} (default), DNA methylation information
#'                  for each available site is returned. Otherwise, this should be one of region types for for which
#'                  summarized DNA methylation information is computed in the given dataset.
#' @param row.names	Flag indicating of row names are to be generated in the result.
#' @param epsilon   Threshold of beta values to use when adjusting for potential M values close to +infinity or
#'                  -infinity. See \code{\link{rnb.beta2mval}} for more details.
#'
#' @return \code{matrix} with methylation M values.
#'
#' @seealso \code{\link[=meth,RnBSet-method]{meth}} for extracting methylation beta values
#' @rdname mval-methods
#' @docType methods
#' @aliases mval
#' @aliases mval,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' ## per-site M-value matrix
#' mm<-mval(rnb.set.example, row.names=TRUE)
#' head(mm)
#' ## M-values for each covered gene
#' gmm<-mval(rnb.set.example, type="gene", row.names=TRUE)
#' head(gmm)
#' }
#' @export
setMethod("mval", signature(object = "RnBSet"),
		  function(object, type = "sites", row.names = FALSE, epsilon = 0) {
		  	beta.values <- get.dataset.matrix(object, type, row.names, object@meth.sites, object@meth.regions)
		  	rnb.beta2mval(beta.values, epsilon)
		  }
)

if(!isGeneric("meth")) setGeneric("meth", function(object, ...) standardGeneric("meth"))
#' meth-methods
#'
#' Extracts DNA methylation information (beta values) for a specified set of genomic features.
#'
#' @param object 	dataset of interest.
#' @param type 		\code{character} singleton. If this is set to \code{"sites"} (default), DNA methylation information
#'                  for each available site is returned. Otherwise, this should be one of region types for for which
#'                  summarized DNA methylation information is computed in the given dataset.
#' @param row.names	flag indicating if row names are to be generated in the result.
#' @param i     	indices of sites/regions to be retrieved. By default (\code{NULL}), all will be retrieved.
#' @param j     	indices of samples to be retrieved. By default (\code{NULL}), all will be retrieved.
#'
#' @return \code{matrix} with methylation beta values.
#'
#' @seealso \code{\link[=mval,RnBSet-method]{mval}} for calculating M values
#' @rdname meth-methods
#' @docType methods
#' @aliases meth
#' @aliases meth,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' ## per-site beta-value matrix
#' mm<-meth(rnb.set.example, row.names=TRUE)
#' head(mm)
#' ## beta-values for each covered gene
#' gmm<-meth(rnb.set.example, type="gene", row.names=TRUE)
#' head(gmm)
#' }
#' @export
setMethod("meth", signature(object = "RnBSet"),
	function(object, type="sites", row.names=FALSE, i=NULL, j=NULL) {
		get.dataset.matrix(object, type, row.names, object@meth.sites, object@meth.regions, i=i, j=j)
	}
)

if(!isGeneric("hasCovg")) setGeneric("hasCovg", function(object,...) standardGeneric("hasCovg"))
#' hasCovg-methods
#'
#' Returns \code{TRUE} if the \code{RnBSet} object contains coverage information for sites or the specified region type.
#'
#' @param object 		\code{RnBSet} of interest.
#' @param type 			\code{character} singleton. If \code{sites} or a region type summarized in the object
#'
#' @return \code{TRUE} if the \code{RnBSet} object contains coverage information for sites or the specified region type. \code{FALSE} otherwise
#'
#' @rdname hasCovg-methods
#' @docType methods
#' @export
#' @aliases hasCovg
#' @aliases hasCovg,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' ## per-site beta-value matrix
#' hasCovg(rnb.set.example)
#' }
setMethod("hasCovg", signature(object="RnBSet"),
	function (object, type="sites") {
		if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
			stop("invalid value for type")
		}
		if (type %in% c("sites", object@target)) {
			result <- !is.null(object@covg.sites)
		} else if (!(type %in% names(object@regions))) {
			stop("unsupported region type")
		} else {
			result <- !is.null(object@covg.regions[[type]])
		}
		return(result)
	}
)


if(!isGeneric("covg")) setGeneric("covg", function(object,...) standardGeneric("covg"))
#' covg-methods
#'
#' Extract coverage information from an object of \code{RnBSet} class.
#'
#' @param object 		Dataset of interest.
#' @param type 			\code{character} singleton. If \code{sites} DNA methylation information per each available
#' 						site is returned. Otherwise should be one of region types for for which the summarized
#' 						coverage information is available
#' @param row.names	    Flag indicating of row names are to be generated in the result.
#' @param i     	indices of sites/regions to be retrieved. By default (\code{NULL}), all will be retrieved.
#' @param j     	indices of samples to be retrieved. By default (\code{NULL}), all will be retrieved.
#'
#' @return coverage information available for the dataset in the form of a \code{matrix}.
#'
#' @rdname covg-methods
#' @docType methods
#' @export
#' @aliases covg
#' @aliases covg,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' ## per-site beta-value matrix
#' cvg<-covg(rnb.set.example, row.names=TRUE)
#' head(cvg)
#' }
setMethod("covg", signature(object="RnBSet"),
	function (object, type="sites", row.names=FALSE, i=NULL, j=NULL) {
		m<-get.dataset.matrix(object, type, row.names, object@covg.sites, object@covg.regions, i=i, j=j)
		m
	}
)

if(!isGeneric("nsites")) setGeneric("nsites", function(object, ...) standardGeneric("nsites"))
#' nsites-methods
#'
#' Returns the number of sites/regions for a given \code{RnBSet} object
#'
#' @param object 	\code{RnBSet} of interest.
#' @param type 		\code{character} singleton. If this is set to \code{"sites"} (default), the number of sites is returned.
#'                  Otherwise, this should be one of region types for for which the number of regions is returned.
#'
#' @return \code{integer} stating the number of sites/regions. \code{NA} if the regions have not been summarized yet.
#'
#' @seealso \code{\link[=meth,RnBSet-method]{meth}} Retrieving the matrix of methylation values
#' @rdname nsites-methods
#' @docType methods
#' @aliases nsites
#' @aliases nsites,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' nsites(rnb.set.example)
#' }
#' @export
setMethod("nsites", signature(object = "RnBSet"),
	function(object, type="sites") {
		if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
			stop("invalid value for type")
		}
		if (type %in% c("sites", object@target)) {
			result <- nrow(object@meth.sites)
		} else if (!(type %in% names(object@regions))) {
			stop("unsupported region type")
		} else if (is.null(object@meth.regions[[type]])) {
			result <- NA
		} else {
			result <- nrow(object@meth.regions[[type]])
		}
		return(result)
	}
)

########################################################################################################################

if (!isGeneric("assembly")) {
	setGeneric("assembly", function(object) standardGeneric("assembly"))
}

#' assembly-methods
#'
#' Extracts information about assembly
#'
#' @param object Dataset of interest.
#' @return Sample annotation information available for the dataset in the form of a \code{data.frame}.
#'
#' @rdname assembly-methods
#' @docType methods
#' @aliases assembly
#' @aliases assembly,RnBSet-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' assembly(rnb.set.example) # "hg19"
#' }
setMethod("assembly", signature(object="RnBSet"),
		  function(object){
		  	return(object@assembly)
		  })


## ---------------------------------------------------------------------------------------------------------------------
## MODIFIERS
## ---------------------------------------------------------------------------------------------------------------------

if (!isGeneric("updateRegionSummaries")) {
	setGeneric("updateRegionSummaries", function(object) standardGeneric("updateRegionSummaries"))
}

#' updateRegionSummaries
#'
#' Updates the region information present in an RnBSet by invoking summarize.regions on all region types
#' present in the object
#'
#' @param object Dataset of interest.
#' @return Sample annotation information available for the dataset in the form of a \code{data.frame}.
#'
#' @rdname updateRegionSummaries
#' @docType methods
#' @aliases updateRegionSummaries
#' @aliases updateRegionSummaries,RnBSet-method
#' @export
setMethod("updateRegionSummaries", signature(object="RnBSet"),
		function(object){
			if (length(object@meth.regions) != 0) {
				region.types <- names(object@meth.regions)
				aggregations <- sapply(object@meth.regions, attr, "aggregation")
				for (i in 1:length(region.types)) {
					object <- summarize.regions(object, region.types[i], aggregations[i])
				}
			}
			object
		}
)

########################################################################################################################

if (!isGeneric("remove.sites")) {
	setGeneric("remove.sites", function(object, probelist, verbose = TRUE) standardGeneric("remove.sites"))
}

#' remove.sites-methods
#'
#' Removes the specified probes from the dataset.
#'
#' @param object    Dataset of interest.
#' @param probelist List of probes to be removed in the form of a \code{logical}, \code{integer} or \code{character}
#'                  vector. If this parameter is \code{logical}, it is not recycled; its length must be equal to the
#'                  number of probes in \code{object}. If it is \code{integer} or \code{character}, it must list only
#'                  probes that exist in the dataset. Specifying probe indices larger than the number of probes, or
#'                  non-existent probe identifiers results in an error.
#' @param verbose	if \code{TRUE} additional diagnostic output is generated
#'
#' @return The modified dataset.
#'
#' @seealso \code{\link[=remove.samples,RnBSet-method]{remove.samples}} for removing samples from a methylation dataset
#'
#' @rdname remove.sites-methods
#' @aliases remove.sites
#' @aliases remove.sites,RnBSet-method
#' @docType methods
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' print(rnb.set.example)
#' ## remove 100 random sites
#' s2r<-sample.int(nrow(sites(rnb.set.example)), 100)
#' rnb.set.f<-remove.sites(rnb.set.example, s2r)
#' print(rnb.set.f)
#' }
setMethod("remove.sites", signature(object = "RnBSet"),
		function(object, probelist, verbose=FALSE) {
			inds <- get.i.vector(probelist, rownames(object@sites))
			if(verbose) {
				rnb.logger.start("Removing sites")
			}
			## Delete methylation sites
			if(length(inds) != 0) {
					object@sites <- object@sites[-inds, ]
					if(!is.null(object@status) && object@status$disk.dump){
						doBigFf <- !is.null(object@status$disk.dump.bigff)
						bff.finalizer <- NULL
						if (doBigFf) doBigFf <- object@status$disk.dump.bigff
						if (doBigFf) bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")
						nSites.new <- nrow(object@meth.sites) - length(inds)
						nSamples <- length(samples(object))
						# methylation
						newMat <- NULL
						if (doBigFf){
							newMat <- BigFfMat(row.n=nSites.new, col.n=nSamples, row.names=NULL, col.names=samples(object), finalizer=bff.finalizer)
						} else {
							newMat <- ff(NA, dim=c(nSites.new, nSamples), dimnames=list(NULL, samples(object)), vmode="double")
						}
						for (j in 1:nSamples){
							newMat[,j] <- object@meth.sites[-inds,j]
						}
						if(isTRUE(object@status$discard.ff.matrices)){
							delete(object@meth.sites)
						}
						object@meth.sites <- newMat

						# coverage
						if(!is.null(object@covg.sites)) {
							newMat <- NULL
							if (doBigFf){
								newMat <- BigFfMat(row.n=nSites.new, col.n=nSamples, row.names=NULL, col.names=samples(object), na.prototype=as.integer(NA), finalizer=bff.finalizer)
							} else {
								newMat <- ff(NA_integer_, dim=c(nSites.new, nSamples), dimnames=list(NULL, samples(object)))
							}
							for (j in 1:nSamples){
								newMat[,j] <- object@covg.sites[-inds,j]
							}
							if(isTRUE(object@status$discard.ff.matrices)){
								delete(object@covg.sites)
							}
							object@covg.sites <- newMat
						}
					} else {
						object@meth.sites <- object@meth.sites[-inds, ,drop=FALSE]
						if(!is.null(object@covg.sites)) {
							object@covg.sites <- object@covg.sites[-inds, ,drop=FALSE]
						}
					}
			}

			## Update region methylation
			if(length(object@meth.regions) != 0){
				region.types <- names(object@meth.regions)
				aggregations <- sapply(object@meth.regions, attr, "aggregation")
				for(i in 1:length(region.types)){
					if(verbose){
						rnb.status(c("summarizing regions:",region.types[i]))
					}
					object <- summarize.regions(object, region.types[i], aggregations[i])
				}
			}

			## Remove information on inferred covariates (they are likely to change when sites are removed)
			if (.hasSlot(object, "inferred.covariates")) {
				i.covariates <- setdiff(names(object@inferred.covariates), "gender")
				if (length(i.covariates) != 0) {
					object@inferred.covariates[i.covariates] <- NULL
					if(verbose){
						rnb.info("removed information on inferred covariates")
					}
				}
			}
			if(verbose){
				rnb.logger.completed()
			}
			object
		}
)

########################################################################################################################

if (!isGeneric("mask.sites.meth")) {
	setGeneric("mask.sites.meth", function(object, mask, verbose=FALSE) standardGeneric("mask.sites.meth"))
}

#' mask.sites.meth-methods
#'
#' Given a logical matrix, sets corresponding entries in the methylation table to NA (masking).
#' Low memory footprint
#'
#' @param object    Dataset of interest.
#' @param mask      logical matrix indicating which sites should be masked
#' @param verbose	if \code{TRUE} additional diagnostic output is generated
#'
#' @return The modified dataset.
#'
#' @rdname mask.sites.meth-methods
#' @aliases mask.sites.meth
#' @aliases mask.sites.meth,RnBSet-method
#' @docType methods
setMethod("mask.sites.meth", signature(object = "RnBSet"),
	function(object, mask, verbose=FALSE) {
		if(!is.null(object@status) && object@status$disk.dump){
			nSamples <- length(samples(object))
			for (j in 1:nSamples){
				object@meth.sites[mask[,j],j] <- NA
			}
		} else {
			object@meth.sites[,][mask] <- NA
			if(inherits(object, "RnBeadRawSet")){
				object@M[,][mask] <- NA
				object@U[,][mask] <- NA
				if(!is.null(object@M0)){
					object@M0[,][mask] <- NA
				}
				if(!is.null(object@U0)){
					object@U0[,][mask] <- NA
				}
			}
		}
		object
	}
)

########################################################################################################################

if (!isGeneric("remove.samples")) {
	setGeneric("remove.samples", function(object, samplelist) standardGeneric("remove.samples"))
}

#' remove.samples-methods
#'
#' Removes the specified samples from the dataset.
#'
#' @param object     Dataset of interest.
#' @param samplelist List of samples to be removed in the form of a \code{logical}, \code{integer} or \code{character}
#'                   vector. If this parameter is \code{logical}, it is not recycled; its length must be equal to the
#'                   number of samples in \code{object}. If it is \code{integer} or \code{character}, it must list only
#'                   samples that exist in the dataset. Specifying sample indices larger than the number of samples, or
#'                   non-existent sample identifiers results in an error.
#' @return The modified dataset.
#'
#' @seealso \code{\link[=remove.sites,RnBSet-method]{remove.sites}} for removing sites or probes from a methylation dataset
#'
#' @rdname remove.samples-methods
#' @aliases remove.samples
#' @aliases remove.samples,RnBSet-method
#' @docType methods
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' samples(rnb.set.example)
#' ## remove 3 random samples
#' s2r<-sample.int(length(samples(rnb.set.example)), 3)
#' rnb.set.f<-remove.samples(rnb.set.example, s2r)
#' samples(rnb.set.f)
#' }
setMethod("remove.samples", signature(object = "RnBSet"),
		function(object, samplelist) {
			object.old <- object
			inds <- get.i.vector(samplelist, samples(object))
			bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")
			if (length(inds) != 0) {
				if(object@status$disk.dump){
					doBigFf <- !is.null(object@status$disk.dump.bigff)
					if (doBigFf) doBigFf <- object@status$disk.dump.bigff

					mat <- object@meth.sites[,]
					new.matrix <- mat[,-inds, drop=FALSE]
					# delete(object@meth.sites)
					if(isTRUE(object@status$discard.ff.matrices)){
						delete(object@meth.sites)
					}
					if (doBigFf){
						object@meth.sites <- BigFfMat(new.matrix, finalizer=bff.finalizer)
					} else {
						object@meth.sites <- convert.to.ff.matrix.tmp(new.matrix)
					}
				}else{
					object@meth.sites <- object@meth.sites[,-inds, drop=FALSE]
				}

				if (!is.null(object@pheno)) {
					object@pheno <- object@pheno[-inds, ,drop=FALSE]
				}
				if (!is.null(object@covg.sites)) {
					if(object@status$disk.dump){
						doBigFf <- !is.null(object@status$disk.dump.bigff)
						if (doBigFf) doBigFf <- object@status$disk.dump.bigff

						mat <- object@covg.sites[,]
						new.matrix <- mat[,-inds, drop=FALSE]
						# delete(object@covg.sites)
						if(isTRUE(object@status$discard.ff.matrices)){
							delete(object@covg.sites)
						}
						if (doBigFf){
							object@covg.sites <- BigFfMat(new.matrix, finalizer=bff.finalizer)
						} else {
							object@covg.sites <- convert.to.ff.matrix.tmp(new.matrix)
						}
					}else{
						object@covg.sites <- object@covg.sites[,-inds, drop=FALSE]
					}
				}
				for (region in names(object@regions)) {
					if(object@status$disk.dump){
						doBigFf <- !is.null(object@status$disk.dump.bigff)
						if (doBigFf) doBigFf <- object@status$disk.dump.bigff

						mat <- object@meth.regions[[region]][,]
						meth.matrix <- mat[, -inds, drop=FALSE]
						if(isTRUE(object@status$discard.ff.matrices)){
							delete(object@meth.regions[[region]])
						}
						if (doBigFf){
							object@meth.regions[[region]] <- BigFfMat(meth.matrix, finalizer=bff.finalizer)
						} else {
							object@meth.regions[[region]] <- convert.to.ff.matrix.tmp(meth.matrix)
						}
						if(!is.null(object@covg.regions)){
							mat <- object@covg.regions[[region]][,]
							covg.matrix <- mat[, -inds, drop=FALSE]
							if(isTRUE(object@status$discard.ff.matrices)){
								delete(object@covg.regions[[region]])
							}
							if (doBigFf){
								object@covg.regions[[region]] <- BigFfMat(covg.matrix, finalizer=bff.finalizer)
							} else {
								object@covg.regions[[region]] <- convert.to.ff.matrix.tmp(covg.matrix)
							}
						}
						# delete(object@meth.regions[[region]])
						# delete(object@covg.regions[[region]])
					}else{
						object@meth.regions[[region]] <- object@meth.regions[[region]][, -inds, drop=FALSE]
						if(!is.null(object@covg.regions)){
							object@covg.regions[[region]] <- object@covg.regions[[region]][, -inds, drop=FALSE]
						}
					}

					attr(object@meth.regions[[region]], "aggregation")<-attr(object.old@meth.regions[[region]], "aggregation")
				}

				## Remove information on inferred covariates (they are likely to change when samples are removed)
				if (.hasSlot(object, "inferred.covariates")) {
					i.covariates <- setdiff(names(object@inferred.covariates), "gender")
					if (length(i.covariates) != 0) {
						## FIXME: Wouldn't it make more sense to simply take the samples out?
						object@inferred.covariates[i.covariates] <- NULL
					}
				}
			}
			object
		}
)

########################################################################################################################

if (!isGeneric("mergeSamples")) {
	setGeneric("mergeSamples", function(object, ...) standardGeneric("mergeSamples"))
}

#' mergeSamples
#'
#' Take an RnBSet object and merge methylation and phenotype information given a grouping column in the pheno table
#' coverage is combined by taking the sum of coverages
#' pheno is combined by concatenating entries from all samples
#' @param object input RnBSet object
#' @param grp.col a column name (string) of \code{pheno(rnb.set)} that contains unique identifiers for sample groups/replicates
#' 		  to be combined
#' @param cov.weighted weight methylation values using coverage
#' @return the modified RnBSet object
#' @details combines phenotype information, coverage information and methylation information
#' methylation is combined by taking the average. Detection p-values are combined using Fisher's method.
#' For methylation arrays, bead counts are currently not taken into account.
#' objects of class \code{RnBeadRawSet} are automatically converted to \code{RnBeadSet}.
#' @note Requires the packages \pkg{foreach} and \pkg{doParallel}.
#'
#' @rdname mergeSamples-methods
#' @aliases mergeSamples
#' @aliases mergeSamples,RnBSet-method
#' @docType methods
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.example
#' rnb.set.merged <- mergeSamples(rnb.set.example,"Cell_Line")
#' rnb.set.merged
#' pheno(rnb.set.merged)
#' }
# TODOs:
# - incorporate weighted methylation average (coverage)
setMethod("mergeSamples", signature(object = "RnBSet"),
	function(object, grp.col, cov.weighted = FALSE){
		ph <- pheno(object)
		if (!is.element(grp.col,colnames(ph))){
			stop("Could not merge samples: phenotype column does not exist")
		}
		res <- object
		replicate.list <- getMergeList(object, grp.col)
		num.replicates <- sapply(replicate.list,length)
		phm <- sapply(ph, format, trim=TRUE, justify="none") #fomat to matrix, avoiding padded whitespaces
		ph.t <- t(phm)
		mf.pheno <- function(X.sub){
			sapply(1:nrow(X.sub),FUN=function(i){
				if (length(unique(X.sub[i,]))==1 && sum(is.na(X.sub[i,]))==0) {
					return(X.sub[i,1])
				} else if (all(is.na(X.sub[i,]))) {
					return(NA)
				} else {
					return(paste(X.sub[i,],collapse=";"))
				}
			})
		}
		pheno.new <- t(mergeColumns(ph.t,replicate.list,mergeFun=mf.pheno))
		pheno.new <- cbind(pheno.new,num.replicates)
		colnames(pheno.new) <- c(colnames(ph),"rnb_number_merged_samples")

		if (class(object) == "RnBiseqSet"){
            if(cov.weighted){
                merge_fun=function(X.sub,iis){rowMeans((X.sub * cov.site.new[,iis]), na.rm=TRUE)/rowSums(cov.site.new[,iis])}
            }else{
                merge_fun=function(X.sub,iis){rowMeans(X.sub, na.rm=TRUE)}
            }
			meth.site.new <- mergeColumns(meth(object,type="sites",row.names=FALSE),replicate.list)
			covg.site.new <- NULL
			if (!is.null(object@covg.sites)){
				covg.site.new <- mergeColumns(covg(object,type="sites"),replicate.list,mergeFun=function(X.sub){rowSums(X.sub,na.rm=TRUE)})
			}
			# res <- new("RnBiseqSet",
			# 		pheno=data.frame(pheno.new),
			# 		sites=object@sites,
			# 		meth.sites=meth.site.new,
			# 		covg.sites=covg.site.new,
			# 		region.types=summarized.regions(object),
			# 		assembly=object@assembly)
			aa <- annotation(object,"sites")
			sites.obj <- data.frame(chrom=as.character(aa$Chromosome),start=aa$Start,strand=as.character(aa$Strand),stringsAsFactors=FALSE)
			doBigFf <- !is.null(object@status$disk.dump.bigff)
			if (doBigFf) doBigFf <- object@status$disk.dump.bigff
			res <- RnBiseqSet(
				pheno=data.frame(pheno.new),
				sites=sites.obj,
				meth=meth.site.new,
				covg=covg.site.new,
				region.types=summarized.regions(object),
				assembly=object@assembly,
				useff=object@status$disk.dump,
				usebigff=doBigFf
			)
		} else if (is.element(class(object),c("RnBeadSet","RnBeadRawSet"))) {
			meth.site.new <- mergeColumns(meth(object,type="sites",row.names=TRUE),replicate.list)
			p.vals <- NULL
			if (!is.null(object@pval.sites)){
				p.vals <- mergeColumns(dpval(object,row.names=TRUE),replicate.list,
					mergeFun=function(X.sub){
						apply(X.sub,1,function(x){combineTestPvalsMeth(na.omit(x),correlated=FALSE)})
					}
				)
			}
			b.counts <- NULL
			if(object@target=="probesEPIC"){
				platform<-"EPIC"
			}else if (object@target=="probes450"){
				platform<-"450k"
			}else if(object@target=="probes27"){
				platform<-"27k"
			}
			# res <- new("RnBeadSet",
			# 	data.frame(pheno.new),
			# 	meth.site.new,
			# 	p.values=p.vals,
			# 	bead.counts=b.counts,
			# 	platform=platform,
			# 	region.types=summarized.regions(object)
			# )
			res <- RnBeadSet(
				pheno=data.frame(pheno.new),
				betas=meth.site.new,
				p.values=p.vals,
				bead.counts=b.counts,
				platform=platform,
				region.types=summarized.regions(object),
				useff=object@status$disk.dump
			)
		} else {
			stop("Could not merge samples: Invalid class of object")
		}
		return(res)
	}
)
########################################################################################################################

#' combine-methods
#'
#' Combine two objects inheriting from \code{\linkS4class{RnBSet}} class
#'
#' @param x,y 		\code{\linkS4class{RnBeadSet}}, \code{\linkS4class{RnBeadRawSet}}
#' 					or \code{\linkS4class{RnBiseqSet}} object
#'
#' @details The sample sets of \code{x} and \code{y} should be unique.
#' Sample annotation information is merged only for columns which have identical names in both objects.
#' CpG sites of the new object are a union of those present in both objects.
#'
#' @return combined \code{\linkS4class{RnBeadSet}}, \code{\linkS4class{RnBeadRawSet}} or
#' \code{\linkS4class{RnBiseqSet}} object
#'
#' @rdname combine-methods
#' @docType methods
#' @export
#' @aliases combine,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' r1 <- rnb.set.example
#' r1 <- remove.samples(r1,samples(rnb.set.example)[1:5])
#' i <- which(r1@@sites[,2] == 15 | r1@@sites[,2] == 21)
#' sites.rem.r1 <- union(sample(1:nrow(meth(rnb.set.example)),500),i)
#' r1 <- remove.sites(r1,sites.rem.r1)
#' r2 <- rnb.set.example
#' r2 <- remove.samples(r2,samples(rnb.set.example)[6:12])
#' sites.rem.r2 <- sample(1:nrow(meth(rnb.set.example)),800)
#' r2 <- remove.sites(r2,sites.rem.r2)
#' rc <- combine(r1,r2)
#' #assertion: check the number of sites
#' sites.rem.c <- intersect(sites.rem.r1,sites.rem.r2)
#' (nrow(meth(rnb.set.example))-length(sites.rem.c)) == nrow(meth(rc))
#' }
setMethod("combine", signature(x="RnBSet",y="RnBSet"),
		function(x,y){
			if (class(x) != class(y)){
				stop("Could not combine RnBSet objects: incompatible classes")
			}
			if (assembly(x) != assembly(y)){
				stop("Could not combine RnBSet objects: incompatible assemblies")
			}
			if (x@target != y@target){
				stop("Could not combine RnBSet objects: incompatible platforms")
			}
			common.samples <- intersect(samples(x),samples(y))
			if (length(common.samples)>0){
				stop(paste0("Could not combine RnBSet objects: the following samples overlap in both objects: ", paste(common.samples,collapse=",")))
			}
			useff <- x@status$disk.dump
			if (x@status$disk.dump != y@status$disk.dump){
				warning(paste0("disk dump status of the two objects to combine disagree. Using disk dump: ", useff))
			}
			usebigff <- useff
			if (usebigff) usebigff <- !is.null(x@status$disk.dump.bigff)
			if (usebigff) usebigff <- x@status$disk.dump.bigff
			if(usebigff){
				bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")
			}
			# prepare a new object
			new.set<-x
			# if(nrow(pheno(x))>=nrow(pheno(y))){
			# 	new.set<-y
			# }else{
			# 	new.set<-x
			# }

			new.set@pheno <- plyr::rbind.fill(pheno(x),pheno(y))

			# combine sites
			sites1<-x@sites
			sites2<-y@sites

			common.chr<-union(unique(sites1[,2]), unique(sites2[,2]))

			subs1<-list()
			subs2<-list()
			common.sites<-list()

			for(chr in common.chr){
				sts<-sort(union(sites1[sites1[,2]==chr,3],sites2[sites2[,2]==chr,3]))
				subs1[[chr]]<-match(sites1[sites1[,2]==chr,3], sts)
				subs2[[chr]]<-match(sites2[sites2[,2]==chr,3], sts)
				common.sites[[chr]]<-cbind(rep(1,length(sts)), rep(chr,length(sts)), sts)
			}

			total.sites<-sum(sapply(common.sites, nrow))

			if("ff_matrix" %in% c(class(sites1), class(sites2))){
				new.sites <- ff(vmode="integer", dim=c(total.sites,3))
				ixx<-1
				for(sts in common.sites){
					new.sites[ixx:(ixx+nrow(sts)),]<-sts
					ixx<-ixx+nrow(sts)+1
				}
			}else{
				new.sites<-do.call("rbind", common.sites)
			}

			colnames(new.sites)<-NULL
			new.set@sites<-new.sites

			slot.names<-RNBSET.SLOTNAMES

			if(inherits(x, "RnBeadSet")){
				slot.names<-c(slot.names, RNBEADSET.SLOTNAMES)
			}
			if(inherits(x, "RnBeadRawSet")){
				slot.names<-c(slot.names, RNBRAWSET.SLOTNAMES)
			}

			for(sl in slot.names){
				if(all(!is.null(slot(x,sl)),!is.null(slot(y,sl)))){
					if(useff){
						#new.matrix<-ff(vmode=vmode(slot(x,sl)), dim=c(total.sites,nrow(pheno(new.set))))
						if (usebigff){
							new.matrix <- BigFfMat(row.n=total.sites, col.n=nrow(pheno(new.set)), vmode=vmode(slot(x,sl)), finalizer=bff.finalizer)
						} else {
							new.matrix <- create.empty.ff.matrix.tmp(vm=vmode(slot(x,sl)), dim=c(total.sites,nrow(pheno(new.set))))
						}
					}else{
						new.matrix<-matrix(NA, nrow=total.sites, ncol=nrow(pheno(new.set)))
					}
					for(chr in common.chr){
						#new.matrix[new.sites[,2]==chr,1:nrow(pheno(x))]<-slot(x,sl)[sites1[,2]==chr,][subs1[[chr]],]
						#new.matrix[new.sites[,2]==chr,(nrow(pheno(x))+1):nrow(pheno(new.set))]<-slot(y,sl)[sites2[,2]==chr,][subs2[[chr]],]
						ix<-which(new.sites[,2]==chr)
						new.matrix[ix[subs1[[chr]]],1:nrow(pheno(x))]<-slot(x,sl)[sites1[,2]==chr,,drop=FALSE]
						new.matrix[ix[subs2[[chr]]],(nrow(pheno(x))+1):nrow(pheno(new.set))]<-slot(y,sl)[sites2[,2]==chr,,drop=FALSE]
					}
					#colnames(new.matrix)<-c(colnames(slot(x,sl)), colnames(slot(y,sl)))
					slot(new.set, sl)<-new.matrix
					rm(new.matrix)
					rnb.cleanMem()
				}else{
					slot(new.set, sl)<-NULL
				}

				if(x@status$disk.dump && isTRUE(x@status$discard.ff.matrices)){
					delete(slot(x, sl))
				}
				if(y@status$disk.dump && isTRUE(y@status$discard.ff.matrices)){
					delete(slot(y, sl))
				}
			}

			if(inherits(x,"RnBeadSet")){
				if(all(!is.null(qc(x)),!is.null(qc(y)))){
					cpn<-intersect(rownames(qc(x)$Cy3), rownames(qc(x)$Cy3))
					cy3.new<-cbind(qc(x)$Cy3[cpn,], qc(y)$Cy3[cpn,])
					cy5.new<-cbind(qc(x)$Cy5[cpn,], qc(y)$Cy5[cpn,])
					colnames(cy3.new)<-NULL
					colnames(cy5.new)<-NULL
					new.set@qc<-list(Cy3=cy3.new, Cy5=cy5.new)
				}else{
					new.set@qc<-NULL
				}
			}

			new.set@status<-list()
			if(inherits(new.set, "RnBeadSet")){
				new.set@status$normalized<-"none"
				new.set@status$background<-"none"
			}
			new.set@status$disk.dump<-useff
			new.set@status$disk.dump.bigff<-usebigff

			for (region.type in union(summarized.regions(x), summarized.regions(y))) {
				if (region.type %in% rnb.region.types(assembly(new.set))) {
					new.set <- summarize.regions(new.set, region.type)
				}
			}
			new.set@inferred.covariates<-list()

			rm(common.sites, sites1, sites2, subs1, subs2, x, y)
			rnb.cleanMem()
			new.set
		}
)

########################################################################################################################
if (!isGeneric("addPheno")) {
	setGeneric("addPheno", function(object, ...) standardGeneric("addPheno"))
}
#' addPheno
#'
#' Adds phenotypic or processing information to the sample annotation table of the given \code{RnBSet} object.
#'
#' @param object \code{\linkS4class{RnBSet}} of interest.
#' @param trait   Trait as a non-empty \code{vector} or \code{factor}. The length of this vector must be equal to the
#'                number of samples in \code{object}, the i-th element storing the value for the i-th sample. Note that
#'                names, if present, are ignored.
#' @param header  Trait name given as a one-element \code{character}. This is the heading to be used for the sample
#'                annotation table. This method fails if such a trait already exists; in other words, if
#'                \code{header \%in\% names(pheno(object))}.
#' @return The modified dataset as an object of type \code{\linkS4class{RnBSet}}.
#'
#' @author Fabian Mueller
#' @export
#' @docType methods
#' @rdname addPheno-RnBSet-methods
#' @aliases addPheno
#' @aliases addPheno,RnBSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' is.hiPSC <- pheno(rnb.set.example)[, "Sample_Group"]=="hiPSC"
#' rnb.set.mod <- addPheno(rnb.set.example, is.hiPSC, "is_hiPSC")
#' pheno(rnb.set.mod)
#' }
setMethod("addPheno", signature(object="RnBSet"),
	function(object, trait, header) {
		if (!((is.vector(trait) || is.factor(trait)) && length(trait) == nrow(pheno(object)))) {
			stop(paste("invalid value for trait; expected vector of length", nrow(pheno(object))))
		}
		if (!(is.character(header) && length(header) == 1 && (!is.na(header)))) {
			stop("invalid value for header; expected one-element character")
		}
		if (is.element(header, names(pheno(object)))) {
			stop(paste("trait", header, "already exists in the sample annotation table"))
		}

		object@pheno[[header]] <- trait
		return(object)
	}
)
########################################################################################################################


if (!isGeneric("summarize.regions")) {
	setGeneric("summarize.regions", function(object, ...) standardGeneric("summarize.regions"))
}
#' summarize.regions-methods
#'
#' Summarize DNA methylation information for which is present in the \code{RnBSet} object.
#'
#' @param object Dataset of interest.
#' @param region.type Type of the region annotation for which the summarization will be performed or \code{"strands"} for summarizing the methylation values from both strands
#' @param aggregation Operation to summarize the methylation values. Currently supported values are \code{"mean"}, \code{"median"}, \code{"min"}, \code{"max"} and \code{"coverage.weighted"}
#' @param overwrite If \code{TRUE} the existing region-level information for \code{region.type} is discarded
#'
#' @return object of the same class as the supplied one containing the summarized methylation information for the specified region types
#'
#' @rdname summarize.regions-methods
#' @docType methods
#' @aliases summarize.regions
#' @aliases summarize.regions,RnBSet-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.set.summarized<-summarize.regions(rnb.set.example, "genes", overwrite=TRUE)
#' head(meth(rnb.set.summarized, type="genes", row.names=TRUE))
#' }
setMethod("summarize.regions", signature(object="RnBSet"),
		function(object, region.type, aggregation = rnb.getOption("region.aggregation"), overwrite = TRUE) {
			if (!(is.character(region.type) && length(region.type) == 1 && (!is.na(region.type)))) {
				stop("invalid value for region.type")
			}
			if (!(is.character(aggregation) && length(aggregation) == 1 && (!is.na(aggregation)))) {
				stop("invalid value for aggregation; expected single character")
			}
			## FIXME: Some of these aren't implemented; and I need them (min and max in particular)
			##        Is there a measurable improvement over the simple get(...) implementation that was dropped?
			aggregation <- aggregation[1]
			if (!(aggregation %in% c("min", "max", "mean", "median", "sum", "coverage.weighted"))) {
				stop("invalid value for aggregation; expected one of \"min\", \"max\", \"mean\", \"median\", \"sum\" or \"coverage.weighted\"")
			}
			if (overwrite == FALSE && region.type %in% names(object@meth.regions)) {
				stop("invalid region type; methylation data already present")
			}
#			aggregate.f <- get(aggregation)
#			aggregate.function <- function(x) {
#				tryCatch(aggregate.f(x, na.rm = TRUE), warning = function(w) { as.double(NA) })
#			}
			## Extract the full annotation tables for the regions and the sites
			if (!(region.type %in% c(rnb.region.types(object@assembly),"strands"))){
				stop("unsupported region type")
			}

			if (region.type =="strands" && !inherits(object, "RnBiseqSet")){
				stop("cannot summarize the strand-specific information for objects other than RnBiseqSet")
			}

			if (aggregation == "coverage.weighted" && !inherits(object, "RnBiseqSet")){
				stop("coverage.weighted aggregation is allowed only for objects of type RnBiseqSet")
			}

			if (aggregation == "coverage.weighted" && is.null(object@covg.sites)){
				stop("cannot apply coverage.weighted aggregation method to an RnBiseqSet object with
						missing coverage information")
			}
			bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")

			if(region.type=="strands"){
				annot.sizes <- rnb.annotation.size(assembly=object@assembly)
				mapping <- sapply(names(rnb.get.chromosomes(assembly=object@assembly)), function(chr){
					num.sites <- annot.sizes[[chr]]
					#TODO:this is not really robust
					IRanges(start=(1:(num.sites/2))*2-1, width=2, names=(1:(num.sites/2))*2-1)
				})
			}else{
				mapping <- rnb.get.mapping(region.type, object@target, object@assembly)
			}

			chromInds <- unique(object@sites[,2])
			#construct the overlap data structure for retrieving other information
			regMap.ov.str <- lapply(chromInds, function(chr.id){
				chr.map <- object@sites[,2]==chr.id
				names(chr.map) <- NULL
				site.ranges <- IRanges(start=object@sites[chr.map,3], width=1)
				chr.name <- names(rnb.get.chromosomes(assembly=object@assembly))[chr.id]
				mapping.contains.chrom <- chr.name %in% names(mapping)
				if(!mapping.contains.chrom){
					return(NULL)
				}
				chr.mapping.ind <- match(chr.name,names(mapping))
				olap <- IRanges::as.matrix(findOverlaps(mapping[[chr.mapping.ind]], site.ranges))

				if(nrow(olap)<1) return(NULL)
				return(list(
					chr.id=chr.id,
					chr.name=chr.name,
					chr.mapping.ind=chr.mapping.ind,
					chr.match.inds=which(chr.map),
					olap=olap
				))
			})
			# logger.info(c("DEBUG:","Generated mapping structure for all chromosomes"))
			region.indices <- do.call("rbind", lapply(regMap.ov.str, function(x){
				if (is.null(x)) return(NULL)
				indOnChrom <- unique(x$olap[,1])
				regInd <- as.integer(names(mapping[[x$chr.mapping.ind]][indOnChrom]))
				cbind(rep(1, length(regInd)), rep(x$chr.id, length(regInd)), regInd)
			}))
			# logger.info(c("DEBUG:","Generated region index data frame"))
			regions2sites <- unlist(lapply(regMap.ov.str, function(x){
				if (is.null(x)) return(list())
				tapply(x$chr.match.inds[x$olap[,2]], factor(x$olap[,1], levels=unique(x$olap[,1])), list)
			}), recursive=FALSE)
			names(regions2sites) <- NULL
			# regions2sites.tab <- do.call("rbind",lapply(1:length(regions2sites), FUN=function(i){
			# 	cbind(rep(i, length(regions2sites[[i]])), regions2sites[[i]])
			# }))
			# regions2sites.tab.fac <- factor(regions2sites.tab[,1], levels=unique(regions2sites.tab[,1]))
			# logger.info(c("DEBUG:","Generated mapping of regions to sites"))

			nSamples <- length(samples(object))

			aggr.f <- NULL
			if (aggregation=="mean"){
				aggr.f <- function(siteInds, siteVec, covgVec=NULL){
					mean(siteVec[siteInds], na.rm=TRUE)
					# 0.666
				}
			} else if (is.element(aggregation, c("min", "max", "mean", "median", "sum"))){
				aggr.f <- function(siteInds, siteVec, covgVec=NULL){
					do.call(aggregation, list(siteVec[siteInds], na.rm=TRUE))
				}
			} else if (aggregation=="coverage.weighted"){
				aggr.f <- function(siteInds, siteVec, covgVec){
					cTotal <- sum(covgVec[siteInds], na.rm=TRUE)
					sum(siteVec[siteInds]*covgVec[siteInds], na.rm=TRUE)/cTotal
				}
			}
			site.meth <- object@meth.sites
			site.covg <- object@covg.sites
			aggr.meth.sample <- function(j){
				siteVec <- site.meth[,j]
				covgVec <- NULL
				if (aggregation=="coverage.weighted") covgVec <- site.covg[,j]
				vapply(regions2sites, aggr.f, numeric(1), siteVec=siteVec, covgVec=covgVec)
			}
			aggr.covg.sample <- function(j){
				siteVec <- site.covg[,j]
				vapply(regions2sites, function(siteInds){
					sum(siteVec[siteInds], na.rm=TRUE)
				}, numeric(1))
			}

			## Assign the resulting matrices to the object
			if (region.type=="strands"){
				if(!is.null(object@status) && object@status$disk.dump){
					doBigFf <- !is.null(object@status$disk.dump.bigff)
					if (doBigFf) doBigFf <- object@status$disk.dump.bigff

					# delete(object@meth.sites)
					if (doBigFf) {
						object@meth.sites <- BigFfMat(row.n=nrow(region.indices), col.n=nSamples, col.names=samples(object), finalizer=bff.finalizer)
						# logger.info(c("DEBUG:","Created BigFfMat for meth"))
					} else {
						object@meth.sites <- convert.to.ff.matrix.tmp(matrix(numeric(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object))))
					}
				} else{
					object@meth.sites <- matrix(numeric(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object)))
				}
				for (j in 1:nSamples){
					# logger.info(c("DEBUG:","Aggregating methylation for sample",j))
					object@meth.sites[,j] <- aggr.meth.sample(j)
				}
				if (!is.null(object@covg.sites)) {
					if(!is.null(object@status) && object@status$disk.dump){
						doBigFf <- !is.null(object@status$disk.dump.bigff)
						if (doBigFf) doBigFf <- object@status$disk.dump.bigff

						# delete(object@covg.sites)
						if (doBigFf) {
							object@covg.sites <- BigFfMat(row.n=nrow(region.indices), col.n=nSamples, col.names=samples(object), finalizer=bff.finalizer)
						} else {
							object@covg.sites <- convert.to.ff.matrix.tmp(matrix(integer(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object))))
						}
					} else {
						object@covg.sites <- matrix(integer(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object)))
					}
					for (j in 1:nSamples){
						# logger.info(c("DEBUG:","Aggregating coverage for sample",j))
						object@covg.sites[,j] <- aggr.covg.sample(j)
					}
				} else {
					object@covg.sites <- NULL
				}
				object@sites <- region.indices
			} else if(!is.null(region.indices)){
				if(!is.null(object@status) && object@status$disk.dump){
					doBigFf <- !is.null(object@status$disk.dump.bigff)
					if (doBigFf) doBigFf <- object@status$disk.dump.bigff
					# if(!is.null(object@meth.regions[[region.type]])){
					# 	delete(object@meth.regions[[region.type]])
					# }
					if(rnb.getOption("enforce.destroy.disk.dumps")){
						delete(object@meth.regions[[region.type]])
					}
					if (doBigFf){
						object@meth.regions[[region.type]] <- BigFfMat(row.n=nrow(region.indices), col.n=nSamples, col.names=samples(object), finalizer=bff.finalizer)
						# logger.info(c("DEBUG:","Created BigFfMat for meth"))
					} else {
						object@meth.regions[[region.type]] <- convert.to.ff.matrix.tmp(matrix(numeric(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object))))
					}
				} else {
					object@meth.regions[[region.type]] <- matrix(numeric(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object)))
				}
				for (j in 1:nSamples){
					# logger.info(c("DEBUG:","Aggregating methylation for sample",j))
					object@meth.regions[[region.type]][,j] <- aggr.meth.sample(j)
				}
				if(!is.null(object@covg.sites)) {
					if(!is.null(object@status) && object@status$disk.dump){
						doBigFf <- !is.null(object@status$disk.dump.bigff)
						if (doBigFf) doBigFf <- object@status$disk.dump.bigff
						# if(!is.null(object@covg.regions[[region.type]])) {
						# 	delete(object@covg.regions[[region.type]])
						# }
						if(rnb.getOption("enforce.destroy.disk.dumps")){
							delete(object@covg.regions[[region.type]])
						}
						if (doBigFf){
							if (is.null(object@covg.regions)) object@covg.regions <- list()
							object@covg.regions[[region.type]] <- BigFfMat(row.n=nrow(region.indices), col.n=nSamples, col.names=samples(object), finalizer=bff.finalizer)
						} else {
							object@covg.regions[[region.type]] <- convert.to.ff.matrix.tmp(matrix(integer(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object))))
						}
					}else{
						object@covg.regions[[region.type]] <- matrix(integer(0), nrow=nrow(region.indices), ncol=nSamples, dimnames=list(NULL,samples(object)))
					}
					for (j in 1:nSamples){
						# logger.info(c("DEBUG:","Aggregating coverage for sample",j))
						object@covg.regions[[region.type]][,j] <- aggr.covg.sample(j)
					}
				}else{
					object@covg.regions <- NULL
				}

				attr(object@meth.regions[[region.type]], "aggregation") <- aggregation
				object@regions[[region.type]] <- region.indices
			}else{ #no valid regions found
				object@meth.regions[[region.type]] <- matrix(0L, nrow=0, ncol=ncol(object@meth.sites))
				if(!is.null(object@covg.sites)) object@covg.regions[[region.type]] <- matrix(0L, nrow=0, ncol=ncol(object@meth.sites))
				attr(object@meth.regions[[region.type]], "aggregation") <- aggregation
				object@regions[[region.type]] <- matrix(0L, nrow=0, ncol=3)
			}
			rm(site.meth) #for ff and BigFfMat, the finalizer should be "delete" and thus the objects should be deleted from disk when this function terminates
			rm(site.covg)
			object
		}
)
########################################################################################################################

if (!isGeneric("remove.regions")) {
	setGeneric("remove.regions", function(object, ...) standardGeneric("remove.regions"))
}

#' remove.regions-methods
#'
#' Remove the summarized methylation information for a given region type from an \code{RnBSet} object.
#'
#' @param object Dataset of interest.
#' @param region.type Type of the region annotation for which the summarization should be removed
#'
#' @return object of the same class as the supplied one without the summarized methylation information for the specified region type
#'
#' @rdname remove.regions-methods
#' @docType methods
#' @aliases remove.regions
#' @aliases remove.regions,RnBSet-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' summarized.regions(rnb.set.example)
#' rnb.set.reduced<-remove.regions(rnb.set.example, "genes")
#' summarized.regions(rnb.set.reduced)
#' }
setMethod("remove.regions", signature(object="RnBSet"),
	function(object, region.type) {
		object@regions[[region.type]] <- NULL
		object@meth.regions[[region.type]] <- NULL
		if(!is.null(object@covg.sites)) object@covg.regions[[region.type]] <- NULL
		return(object)
	}
)

########################################################################################################################

if (!isGeneric("regionMapping")) {
	setGeneric("regionMapping", function(object, ...) standardGeneric("regionMapping"))
}

#' regionMapping-methods
#'
#' get the mapping of regions in the RnBSet object to methylation site indices in the RnBSet object
#'
#' @param object Dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param region.type  region type. see \code{\link{rnb.region.types}} for possible values
#' @return A list containing for each region the indices (as integers) of sites that belong to that region
#'
#' @rdname regionMapping-methods
#' @docType methods
#' @aliases regionMapping
#' @aliases regionMapping,RnBSet-method
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' promoter.probe.list <- regionMapping(rnb.set.example,"promoters")
#' #get the number of CpGs per promoter in the dataset:
#' sapply(promoter.probe.list,length)
#' }
setMethod("regionMapping", signature(object = "RnBSet"),
function(object, region.type) {

			if (!inherits(object, "RnBSet")) {
				stop("invalid value for object; expected RnBSet")
			}
			if (!(is.character(region.type) && length(region.type) == 1 && (!is.na(region.type)))) {
				stop("invalid value for type")
			}

			if (!(region.type %in% rnb.region.types(object@assembly))) {
				stop(paste0("unsupported annotation type (annotation): ",region.type))
			}
			if (!(region.type %in% names(object@regions))) {
				stop(paste0("unsupported annotation type (RnBSet): ",region.type))
			}
			chrom.maps <- rnb.get.mapping(region.type, object@target, object@assembly)

			chrom.integer2name <- names(rnb.get.chromosomes(assembly=object@assembly))
			obj.sites <- data.frame(object@sites)
			region.map <- object@regions[[region.type]]
			chr.inds.reg <- unique(region.map[,2])

			obj.sites[,2] <- factor(chrom.integer2name[obj.sites[,2]],levels=chrom.integer2name[unique(obj.sites[,2])])
#			obj.sites[,2] <- factor(chrom.integer2name[obj.sites[,2]],levels=chrom.integer2name)
#			obj.sites[,2] <- as.factor(chrom.integer2name[obj.sites[,2]])
			chrom.site.inds <- tapply(obj.sites[,3],obj.sites[,2],FUN=function(x){
				IRanges(start=x,width=1)
			})

			chrom.offsets <- sapply(chrom.site.inds,length)
			chrom.offsets <-cumsum(c(0,chrom.offsets[-length(chrom.offsets)]))
			names(chrom.offsets) <- names(chrom.site.inds)

			result <- lapply(chr.inds.reg,FUN=function(chr){

				curChromName <- chrom.integer2name[chr]
				rnbs.regs <- region.map[region.map[,2]==chr,3]
				rnbs.regs.char <- format(rnbs.regs,trim=TRUE,scientific=FALSE)
				rrRanges <- chrom.maps[[curChromName]]
				#only take the regions that are also in the RnBSet object
				if (!all(rnbs.regs.char %in% names(rrRanges))) {stop(paste("Not all regions in RnBSet are present in the annotation (",curChromName,")"))}
				rrRanges <- rrRanges[rnbs.regs.char,]
				olap<-as.matrix(findOverlaps(chrom.site.inds[[curChromName]], rrRanges))

				olap[,1]<-olap[,1]+chrom.offsets[curChromName]
				res<-tapply(olap[,1], olap[,2], list)

				return(res)
			})

			result<-unlist(result, recursive=FALSE)
			names(result)<-NULL
			if (dim(region.map)[1] != length(result)){
				stop("regionMapping failed")
			}
			return(result)

		}
)

########################################################################################################################
#' annotation-methods
#'
#' Genomic annotation of the methylation sites or regions covered in the supplied dataset.
#'
#' @param object dataset as an object of type inheriting \code{RnBSet}.
#' @param type   loci or regions for which the annotation should be obtained. If the value of this parameter is
#'               \code{"sites"} (default), individual methylation sites are annotated. Otherwise, this must be one of
#'               the available region types, as returned by \code{\link{rnb.region.types}}.
#' @param add.names flag specifying whether the unique site identifiers should be used as row names of the
#' 					resulting data frame
#' @param include.regions if \code{TRUE} one additional column is added to the returned annotation dat frame
#' 						  for each of the available region types, giving the indices of the
#'
#' @return Annotation table in the form of a \code{data.frame}.
#'
#' @rdname annotation-methods
#' @docType methods
#' @aliases annotation
#' @aliases annotation,RnBSet-method
#' @author Pavlo Lutsik
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' ## show present sites
#' head(annotation(rnb.set.example, add.names=TRUE))
#' ## show promoters
#' ann.prom<-annotation(rnb.set.example, type="promoters", add.names=TRUE)
#' head(ann.prom)
#' }
setMethod("annotation", signature(object = "RnBSet"),
		function(object, type="sites", add.names=FALSE, include.regions=FALSE) {

			if (!inherits(object, "RnBSet")) {
				stop("invalid value for object; expected RnBSet")
			}
			if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
				stop("invalid value for type")
			}

			if (type == "sites") {
				type <- object@target
				subsets <- object@sites
			} else {
				if (!(type %in% rnb.region.types(object@assembly))) {
					stop(paste0("unsupported annotation type (annotation): ",type))
				}
				if (!(type %in% names(object@regions))) {
					## This region type is not initialized with summarize.regions
					## FIXME: Report an error or initialize.
					stop(paste0("unsupported annotation type (RnBSet): ",type))
				}
				subsets <- object@regions[[type]]
			}

			annot <- rnb.get.annotation(type, object@assembly)
			ind.shift<-rnb.annotation.size(type, object@assembly)
			ind.shift<-cumsum(c(0,ind.shift[-length(ind.shift)]))
			subsets.full<-subsets[,3]+ind.shift[subsets[,2]]

			result<-rnb.annotation2data.frame(annot, add.names=add.names)[subsets.full,]

			if(include.regions){
				dump<-sapply(names(object@regions), function(rt){
					result[,rt]<<-rep(0L,nrow(result))
					map<-regionMapping(object, rt)
					index_map<-lapply(1:length(map), function(ix) rep(ix, length(map[[ix]])))
					result[unlist(map),rt]<<-unlist(index_map)
				})
			}

			return(result)
		}
)
########################################################################################################################

#if (!isGeneric("save.matrices")) {
	setGeneric("save.matrices", function(object, path, ...) standardGeneric("save.matrices"))
#}

setMethod("save.matrices", signature(object="RnBSet", path="character"),
		function(object, path){

			if(!is.null(object@status) && object@status$disk.dump){
				if("ff" %in% class(object@meth.sites)){
					ffmatrix <- object@meth.sites
					ffsave(ffmatrix,file=file.path(path, "rnb.meth"),rootpath=getOption('fftempdir'))
					rm(ffmatrix)
				} else if("BigFfMat" %in% class(object@meth.sites)){
					save.bigFfMat(object@meth.sites, file=file.path(path, "rnb.meth"), rootpath=getOption('fftempdir'))
				}

				if("ff" %in% class(object@covg.sites)){
					ffmatrix <- object@covg.sites
					ffsave(ffmatrix, file=file.path(path, "rnb.covg"),rootpath=getOption('fftempdir'))
					rm(ffmatrix)
				} else if("BigFfMat" %in% class(object@covg.sites)){
					save.bigFfMat(object@covg.sites, file=file.path(path, "rnb.covg"), rootpath=getOption('fftempdir'))
				}

				if(length(object@regions) != 0){

					for(rgn in 1:length(object@regions)){

						rgnpath<-file.path(path,rgn)
						if(!file.exists(rgnpath)){
							dir.create(rgnpath)
						}

						if("ff" %in% class(object@meth.regions[[rgn]])){
							ffmatrix<-object@meth.regions[[rgn]]
							ffsave(ffmatrix, file=file.path(path, rgn, "rnb.meth"),rootpath=getOption('fftempdir'))
							rm(ffmatrix)
						} else if("BigFfMat" %in% class(object@meth.regions[[rgn]])){
							save.bigFfMat(object@meth.regions[[rgn]], file=file.path(path, rgn, "rnb.meth"), rootpath=getOption('fftempdir'))
						}

						if("ff" %in% class(object@covg.regions[[rgn]])){
							ffmatrix<-object@covg.regions[[rgn]]
							ffsave(ffmatrix, file=file.path(path, rgn, "rnb.covg"),rootpath=getOption('fftempdir'))
							rm(ffmatrix)
						} else if("BigFfMat" %in% class(object@covg.regions[[rgn]])){
							save.bigFfMat(object@covg.regions[[rgn]], file=file.path(path, rgn, "rnb.covg"), rootpath=getOption('fftempdir'))
						}

					}

				}
			}


		})

########################################################################################################################

setGeneric("load.matrices",
		function(object, path, ...) standardGeneric("load.matrices"))

setMethod("load.matrices", signature(object="RnBSet", path="character"),
		function(object, path, temp.dir=tempdir()){
			doBigFf <- !is.null(object@status)
			if (doBigFf) doBigFf <- !is.null(object@status$disk.dump.bigff)
			if (doBigFf) doBigFf <- object@status$disk.dump.bigff

			if (doBigFf){
				object@meth.sites <- load.bigFfMat(file.path(path, "rnb.meth"), rootpath=getOption("fftempdir"))
				if(!is.null(object@covg.sites)){
					object@covg.sites <- load.bigFfMat(file.path(path, "rnb.covg"), rootpath=getOption("fftempdir"))
				}
			} else {
				if(sum(grepl("rnb.meth", list.files(path)))==2){
					load_env<-new.env()
					suppressMessages(ffload(file=file.path(path, "rnb.meth"), envir=load_env,rootpath=getOption("fftempdir")))
					object@meth.sites<-get("ffmatrix", envir=load_env)
					rm(load_env)
				}

				if(sum(grepl("rnb.covg", list.files(path)))==2){
					load_env<-new.env()
					suppressMessages(ffload(file=file.path(path, "rnb.covg"), envir=load_env,rootpath=getOption("fftempdir")))
					object@covg.sites<-get("ffmatrix", envir=load_env)
					rm(load_env)
				}
			}

			rgns <- names(object@regions)
			if(!is.null(rgns)){
				if (.hasSlot(object, 'version')) {
					rgns <- 1:length(rgns)
				}

				for(rgn in rgns){
					if (doBigFf){
						object@meth.regions[[rgn]] <- load.bigFfMat(file.path(path, rgn, "rnb.meth"), rootpath=getOption("fftempdir"))
						if(!is.null(object@covg.regions[[rgn]])){
							object@covg.regions[[rgn]] <- load.bigFfMat(file.path(path, rgn, "rnb.covg"), rootpath=getOption("fftempdir"))
						}
					} else {
						if(sum(grepl("rnb.meth",list.files(file.path(path, rgn))))==2){
							load_env<-new.env()
							suppressMessages(ffload(file=file.path(path, rgn, "rnb.meth"), envir=load_env, rootpath=getOption("fftempdir")))
							object@meth.regions[[rgn]]<-get("ffmatrix", envir=load_env)
							rm(load_env)
						}
						if(sum(grepl("rnb.covg",list.files(file.path(path, rgn))))==2){
							load_env<-new.env()
							suppressMessages(ffload(file=file.path(path, rgn, "rnb.covg"), envir=load_env, rootpath=getOption("fftempdir")))
							object@covg.regions[[rgn]]<-get("ffmatrix", envir=load_env)
							rm(load_env)
						}
					}
				}
			}

			return(object)

		})

########################################################################################################################

#' save.rnb.set
#'
#' Consistent saving of an \code{RnBSet} objects with large matrices of type \link{ff}.
#'
#' @param object     \code{RnBSet}-inheriting object.
#' @param path	      the name of the output file (or directory if \code{archive} is \code{FALSE})
#' 					  without an extension. If only the file name is given the object will be saved
#' 					  in the current working directory.
#' @param archive     if \code{TRUE} (default value) the output is a ZIP-file.
#'
#' @details 		  The saved object can be reloaded with the \link{load.rnb.set} function.
#'
#' @return 			  invisibly, the full path to the ZIP file (if \code{archive} is \code{TRUE}),
#' 					  or to the output directory (otherwise)
#'
#' @author Pavlo Lutsik
#' @export
save.rnb.set<-function(object, path, archive=TRUE){

	## Validate parameters
	if (!inherits(object, "RnBSet")) {
		stop("invalid value for object")
	}
	if (!(is.character(path) && length(path) == 1 && isTRUE(!grepl("^[/\\.]*$", path)))) {
		stop("invalid value for path")
	}
	if (!parameter.is.flag(archive)) {
		stop("invalid value for archive")
	}

	if(object@status$disk.dump && .Platform$OS == "windows" && Sys.getenv("R_ZIPCMD")==""){
		rnb.warning(c("Zip not found on this Windows system, this RnBSet object will not be saved.",
				"See the instructions for installing ZIP on Windows in the FAQ section of the RnBeads website."))
		return(invisible(path))
	}

	## Get the full path of the file or directory to be created
	fullpath <- normalizePath(gsub("/$", "", gsub("\\", "/", path, fixed = TRUE)), winslash = "/", mustWork = FALSE)

	## Create or overwrite a directory to store the files
	if (unlink(fullpath, recursive = TRUE) == 1) {
		stop("Specified path already exists and cannot be overwritten")
	}
	if (archive) {
		if (unlink(paste0(fullpath, ".zip"), recursive = TRUE) == 1) {
			stop("Specified path already exists and cannot be overwritten")
		}
	}
	if (!dir.create(fullpath, showWarnings = FALSE, recursive = TRUE)) {
		stop("Could not create output directory")
	}

	## Save all data structures
	save.matrices(object, fullpath)
	save(object, file=file.path(fullpath, "rnb.set.RData"))

	## Create a ZIP archive of the whole directory
	if(archive){
		currdir <- setwd(fullpath)
		zip(paste0(fullpath, ".zip"), dir(), flags = "-rm9X")
		while(length(list.files(path))>0){
			TRUE;
		}
		setwd(currdir)
		if (unlink(fullpath, recursive = TRUE) == 1) {
			rnb.warning("Could not clean output directory after zipping")
		}
		fullpath <- paste0(fullpath, ".zip")
	}

	return(invisible(fullpath))
}

########################################################################################################################

#' load.rnb.set
#'
#' Loading of the \code{RnBSet} objects with large matrices of type \pkg{ff}.
#'
#' @param path			full path of the file or directory. If \code{archive} is \code{FALSE})
#' 					  	without an extension.
#' @param temp.dir		\code{character} singleton which specifies temporary directory, used while loading
#'
#' @return				Loaded object
#'
#' @author Pavlo Lutsik
#' @export
load.rnb.set<-function(path, temp.dir=tempdir()){

	## Validate parameters
	if (!(is.character(path) && length(path) == 1 && isTRUE(!grepl("^[/\\.]*$", path)))) {
		stop("invalid value for path")
	}
	if (!(is.character(temp.dir) && length(temp.dir) == 1 && isTRUE(!grepl("^[/\\.]*$", temp.dir)))) {
		stop("invalid value for temp.dir")
	}
	if (!file.exists(path)) {
		stop("invalid value for path; the path does not exist")
	}
	if (!isTRUE(file.info(temp.dir)[1, "isdir"])) {
		stop("invalid value for temp.dir; the path does not exist or is not a directory")
	}

	if(.Platform$OS == "windows" && Sys.getenv("R_ZIPCMD")==""){
	 	method="internal"
	}else{
		method="unzip"
	}

	if(!file.info(path)[["isdir"]]){
		td<-tempfile("extraction", temp.dir)
		unzip(path, exdir=td, unzip=method)
	}else{
		td<-path
	}

	load_env<-new.env(parent=emptyenv())
	load(file.path(td, "rnb.set.RData"),envir=load_env)
	load.matrices(get("object", load_env), td, temp.dir=temp.dir)
}

########################################################################################################################

if (!isGeneric("destroy")) setGeneric("destroy", function(object) standardGeneric("destroy"))
#' destroy-methods
#'
#' Remove tables stored to disk from the file system. Useful for cleaning up disk dumped objects.
#'
#' @param object object inheriting from \code{\linkS4class{RnBSet}}
#' @return Nothing of particular interest
#'
#' @rdname destroy-methods
#' @docType methods
#' @aliases destroy
#' @aliases destroy,RnBSet-method
#' @export
setMethod("destroy", signature(object="RnBSet"),
		function(object){

			if(object@status$disk.dump){
				delete(object@meth.sites)

				if(!is.null(object@covg.sites)){
					delete(object@covg.sites)
				}

				if(!is.null(object@regions)){
					for(rgn in names(object@regions)){
						delete(object@meth.regions[[rgn]])
						if(!is.null(object@covg.regions))
						{
							delete(object@covg.regions[[rgn]])
						}
					}
				}
			}
			return(invisible(TRUE))
		}
)

########################################################################################################################

## meth.matrices
##
## Creates a list of methylation value (beta) matrices for the given dataset.
##
## @param object        Methylation dataset object of type that inherits \code{RnBSet}.
## @param include.sites Flag indicating if the methylation matrix of sites or probes is to be included in the result.
## @return Non-empty \code{list} of matrices of beta values. If \code{include.sites} is \code{TRUE}, the first matrix in
##         the list is the one based on sites or probes. Other matrices store region-based methylation for (some of) the
##         regions addressed in the option \code{"region.types"}.
## @author Yassen Assenov
meth.matrices <- function(object, include.sites = rnb.getOption("analyze.sites")) {
	result <- list()
	if (include.sites) result[["sites"]] <- meth(object)
	for (rtype in rnb.region.types.for.analysis(object)) {
		X <- tryCatch(meth(object, rtype), error = function(e) { NULL })
		if (!is.null(X)) {
			result[[rtype]] <- X
		}
	}
	return(result)
}

########################################################################################################################

## get.row.names
##
## Generates row names based on the genomic location.
##
## @param object \code{RnBSet} object.
## @return \code{character} vector of row names.
## @author Pavlo Lutsik
get.row.names<-function(object, type="sites"){
	if(type=="sites"){
		target<-object@target
		subsets<-object@sites
	}else if(type %in% names(object@regions)){
		target<-type
		subsets<-object@regions[[type]]
	}else stop("unsupported region type")

	loc.info<-annotation(object, type=type, add.names=TRUE)
	if ("ID" %in% colnames(loc.info) && anyDuplicated(loc.info[, "ID"]) == 0) {
		result <- loc.info[,"ID"]
	} else if (!is.null(rownames(loc.info))) {
		result <- rownames(loc.info)
	} else {
		result <- paste(loc.info[,"Chromosome"], loc.info[,"Start"], as.character(loc.info[,"Strand"]), sep=".")
	}
	result
}

########################################################################################################################

## rnb.get.row.token
##
## Gets the methylation target, that is, the basic methylation feature of a dataset based on its platform.
##
## @param object Methylation dataset of interest, an object of type inheriting \code{MethyLumiSet} or \code{RnBSet}.
## @param plural Flag, indicating if the plural form of the word.
## @return Word or phrase denoting the term for a single target of the platform.
## @author Pavlo Lutsik
rnb.get.row.token<-function(object, plural = FALSE){

	if (is.character(object)) {
		result <- ifelse(object %in% c("RnBiseqSet", "RnBSet"), "site", "probe")
	} else if (inherits(object, "MethyLumiSet")){
		result <- "probe"
	} else if (object@target == "CpG") {
		result <- "site"
	} else { # object@target == "probes450"
		result <- "probe"
	}

	ifelse(plural, paste0(result, "s"), result)
}

########################################################################################################################

## rnb.get.covg.token
##
## Gets the measure of coverage of a dataset based on its platform.
##
## @param object  Methylation dataset of interest, an object of type inheriting \code{MethyLumiSet} or \code{RnBSet}.
## @param capital Flag, indicating if the first letter of the returned phrase should be capitalized.
## @return Word or phrase denoting the term for depth of coverage.
## @author Pavlo Lutsik
rnb.get.covg.token<-function(object, capital=FALSE){

	if (is.character(object)) {
		result <- ifelse(object %in% c("RnBiseqSet", "RnBSet"), "coverage", "bead counts")
	} else if (inherits(object, "MethyLumiSet")) {
		result <- "bead counts"
	} else if (object@target == "CpG") {
		result <- "coverage"
	} else { # object@target == "probes450"
		result <- "bead counts"
	}

	ifelse(capital, capitalize(result), result)
}

########################################################################################################################
if(!isGeneric("sampleMethApply")) setGeneric("sampleMethApply", function(object, ...) standardGeneric("sampleMethApply"))
#' sampleMethApply-methods
#'
#' Applies a function over the methylation values for all samples in an \code{RnBSet} using a low memory footprint.
#' 
#' @param object object inheriting from \code{\linkS4class{RnBSet}}
#' @param fn function to be applied
#' @param type \code{character} singleton. Specify "sites" (default) or a region type over which the function is applied
#' @param ... arguments passed on to the function
#' @return Result analogous to \code{apply(meth(rnbSet, type), 2, FUN=FUN)}
#'
#' @seealso \code{\link[=meth,RnBSet-method]{meth}} Retrieving the matrix of methylation values
#' @rdname sampleMethApply-methods
#' @docType methods
#' @aliases sampleMethApply
#' @aliases sampleMethApply,RnBSet-method
setMethod("sampleMethApply", signature(object = "RnBSet"),
	function(object, fn, type="sites", ...) {
		if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
			stop("invalid value for type")
		}
		if (type %in% c("sites", object@target)) {
			result <- nrow(object@meth.sites)
		} else if (!(type %in% names(object@regions))) {
			stop("unsupported region type")
		}
		res <- sapply(1:length(samples(object)), FUN=function(j){
			fn(meth(object, type=type, j=j), ...)
		})
		return(res)
	}
)
if(!isGeneric("sampleCovgApply")) setGeneric("sampleCovgApply", function(object, ...) standardGeneric("sampleCovgApply"))
#' sampleCovgApply-methods
#'
#' Applies a function over the coverage values for all samples in an \code{RnBSet} using a low memory footprint.
#' @param object object inheriting from \code{\linkS4class{RnBSet}}
#' @param fn function to be applied
#' @param type \code{character} singleton. Specify "sites" (default) or a region type over which the function is applied
#' @param ... arguments passed on to the function
#' @return Result analogous to \code{apply(covg(rnbSet, type), 2, FUN=FUN)}
#'
#' @seealso \code{\link[=meth,RnBSet-method]{covg}} Retrieving the matrix of coverage values
#' @rdname sampleCovgApply-methods
#' @docType methods
#' @aliases sampleCovgApply
#' @aliases sampleCovgApply,RnBSet-method
setMethod("sampleCovgApply", signature(object = "RnBSet"),
	function(object, fn, type="sites", ...) {
		if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
			stop("invalid value for type")
		}
		if (type %in% c("sites", object@target)) {
			result <- nrow(object@covg.sites)
		} else if (!(type %in% names(object@regions))) {
			stop("unsupported region type")
		}
		res <- sapply(1:length(samples(object)), FUN=function(j){
			fn(covg(object, type=type, j=j), ...)
		})
		return(res)
	}
)
########################################################################################################################
if(!isGeneric("getNumNaMeth")) setGeneric("getNumNaMeth", function(object, ...) standardGeneric("getNumNaMeth"))
#' getNumNaMeth-methods
#'
#' for each site/region, the getNumNaMeth retrieves the number of NA values accross all samples.
#' Does this efficiently by breaking down the methylation matrix into submatrices
#' @param object object inheriting from \code{\linkS4class{RnBSet}}
#' @param type "sites" or region type
#' @param chunkSize size of each submatrix (performance tuning parameter)
#' @param mask logical matrix. its entries will also be considered NAs in counting
#' @return vector containing the number of NAs per site/region
#'
#' @rdname getNumNaMeth-methods
#' @docType methods
#' @aliases getNumNaMeth
#' @aliases getNumNaMeth,RnBSet-method
setMethod("getNumNaMeth", signature(object = "RnBSet"),
	function(object, type="sites", chunkSize=1e5, mask=NULL) {
		if (!(is.character(type) && length(type) == 1 && (!is.na(type)))) {
			stop("invalid value for type")
		}
		if (!(type %in% c("sites", object@target, names(object@regions)))) {
			stop("unsupported region type")
		}
		#get start and end indices for the chunks
		n <- nsites(object, type)
		indStarts <- seq(1,n,by=chunkSize)
		indEnds <- c(indStarts[-1]-1, n)
		#apply to each chunk
		res <- unlist(lapply(1:length(indStarts), FUN=function(i){
			indsCur <- indStarts[i]:indEnds[i]
			mm <- meth(object, type=type, i=indsCur)
			isNaMat <- is.na(mm)
			if (!is.null(mask)) isNaMat <- isNaMat | mask[indsCur,]
			return(as.integer(rowSums(isNaMat)))
		}))
		
		return(res)
	}
)
########################################################################################################################
#' rnb.sample.summary.table
#'
#' Creates a sample summary table from an RnBSet object
#'
#' @param rnbSet \code{\linkS4class{RnBSet}} of interest.
#' @return a summary table (as data.frame) with the following variables for each sample (rows):
#' \item{sampleName}{Name of the sample}
#' \item{*_num (* can be 'sites' or a region type)}{Number of sites or regions with coverage in the sample}
#' \item{*_covgMean (\code{RnBiseqSet} only)}{Mean coverage of sites or regions in the sample}
#' \item{*_covgMedian (\code{RnBiseqSet} only)}{Median coverage of sites or regions in the sample}
#' \item{*_covgPerc25 (\code{RnBiseqSet} only)}{25 percentile of coverage of sites or regions in the sample}
#' \item{*_covgPerc75 (\code{RnBiseqSet} only)}{75 percentile of coverage of sites or regions in the sample}
#' \item{*_numCovg5,10,30,60 (\code{RnBiseqSet} only)}{Number of sites or regions with coverage greater or equal to 5,10,30,60}
#' \item{sites_numDPval5em2,1em2,1em3 (\code{RnBeadSet} only)}{Number of sites with a detection p-value smaller than 0.05,0.01,0.001}
#' \item{**_numSitesMean (** is any region type)}{Mean number of sites in a region}
#' \item{**_numSitesMedian}{Median number of sites in a region}
#' \item{**_numSites2,5,10,20}{Number of regions with at least 2,5,10,20 sites with valid methylation measurements}
#' @author Fabian Mueller
#' @aliases rnb.sample.summary.table,RnBSet-method
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.sample.summary.table(rnb.set.example)
#' }
rnb.sample.summary.table <- function(rnbSet) {
	is.biseq <- "RnBiseqSet" %in% class(rnbSet)
	is.beads <- "RnBeadSet" %in% class(rnbSet)
	df.empty <- data.frame(matrix(nrow=length(samples(rnbSet)),ncol=0))
	rownames(df.empty) <- samples(rnbSet)
	tt <- data.frame(df.empty,sampleName=samples(rnbSet),stringsAsFactors=FALSE)
	reg.types.regions <- summarized.regions(rnbSet)
	reg.types <- c("sites",reg.types.regions)
	for (rr in reg.types){
		# logger.status(c("Region type:",rr))
		# rnb.cleanMem()
		tt.cur <- df.empty
		tt.cur$num <- sampleMethApply(rnbSet, function(x){sum(!is.na(x))}, type=rr)
		if (is.biseq){
			covgStats <- do.call("rbind", lapply(1:length(samples(rnbSet)), FUN=function(j){
				# logger.status(c("  Sample:",j))
				mm <- as.vector(meth(rnbSet, rr, j=j))
				cc <- as.vector(covg(rnbSet, rr, j=j))
				cc[cc==0] <- NA
				cc[is.na(mm)] <- NA
				qq <- quantile(cc, probs = c(0.25,0.75), na.rm=TRUE)
				res <- c(
					mean(cc, na.rm=TRUE),
					median(cc, na.rm=TRUE),
					qq[1],
					qq[2],
					sum(cc>=5, na.rm=TRUE),
					sum(cc>=10, na.rm=TRUE),
					sum(cc>=30, na.rm=TRUE),
					sum(cc>=60, na.rm=TRUE)
				)
				return(res)
			}))
			colnames(covgStats) <- c("covgMean", "covgMedian", "covgPerc25", "covgPerc75", "numCovg5", "numCovg10", "numCovg30", "numCovg60")
			tt.cur <- cbind(tt.cur, covgStats)

		}
		if (is.beads){
			if (rr == "sites"){
				pp <- dpval(rnbSet,type=rr)
				if (!is.null(pp)) {
					tt.cur$numDPval5em2 <- colSums(pp < 5e-2, na.rm=TRUE)
					tt.cur$numDPval1em2 <- colSums(pp < 1e-2, na.rm=TRUE)
					tt.cur$numDPval1em3 <- colSums(pp < 1e-3, na.rm=TRUE)
				}
			}
		}

		if (rr %in% reg.types.regions){
			regions2sites <- regionMapping(rnbSet,region.type=rr)
			#compute the number of sites per region and sample
			nsamples <- length(samples(rnbSet))
			num.sites <- sapply(1:nsamples,function(i){
				# logger.status(c("  Sample:",i))
				mm.s.nna <- !is.na(as.vector(meth(rnbSet, j=i)))
				sapply(1:nsites(rnbSet, rr),function(j){
					sum(mm.s.nna[regions2sites[[j]]])
				})
			})
			# num.sites2 <- t(sapply(1:nsites(rnbSet, rr),function(i){
			# 	# logger.status(c("  Site/Region:",i))
			# 	colSums(!is.na(meth(rnbSet, i=regions2sites[[i]])))
			# })) # a bit slower, but more memory effective, if to include later, check the code again
			tt.cur$numSitesMean   <- colMeans(num.sites, na.rm=TRUE)
			tt.cur$numSitesMedian <- colMedians(num.sites, na.rm=TRUE)
			tt.cur$numSites2  <- colSums(num.sites>=2, na.rm=TRUE)
			tt.cur$numSites5  <- colSums(num.sites>=5, na.rm=TRUE)
			tt.cur$numSites10 <- colSums(num.sites>=10,na.rm=TRUE)
			tt.cur$numSites20 <- colSums(num.sites>=20,na.rm=TRUE)
		}
		colnames(tt.cur) <- paste(rr,colnames(tt.cur),sep="_")

		tt <- data.frame(tt,tt.cur)
	}
	return(tt)
}
########################################################################################################################
