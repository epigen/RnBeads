########################################################################################################################
## RnBeadSet-class.R
## created: 2012-04-06
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBeadSet class definition.
########################################################################################################################
## GLOBALS

RNBEADSET.SLOTNAMES<-c("pval.sites")

## ---------------------------------------------------------------------------------------------------------------------
## CLASS DEFINITIONS
## ---------------------------------------------------------------------------------------------------------------------

#' RnBeadSet Class
#'
#' Stores the preprocessed information from HumanMethylation experiments
#'
#' @details
#' There are multiple ways to create an object of type \code{RnBeadSet}:
#' \describe{
#'  \item{Loading from files}{Dataset can be loaded from text or binary files. See the function
#'       \code{\link{rnb.execute.import}} for more details.}
#'  \item{Downloading from GEO}{See the function \code{\link{rnb.read.geo}} for details.}
#'  \item{Converting from \code{MethyLumiSet}}{...}
#' } 
#'
#' @section Slots:
#' \describe{
#'   \item{\code{pval.sites}}{\code{matrix} of detection p-values with the same dimensions as \code{betas}, or
#'        \code{NULL} if the detection p-values are not available.}
#'   \item{\code{pval.regions}}{\code{list} of methylation \code{matrix} objects, one per available region type. Every row in a 
#' 		matrix corresponds to a methylation site, and every column - to a sample.}
#'   \item{\code{covg.sites}}{\code{matrix} of bead counts per probe with the same dimensions as \code{betas}, or
#'        \code{NULL} if this data are not available.}
#'   \item{\code{qc}}{Quality control probe information in the form of a \code{list} of two elements - \code{"Cy3"} and
#'        \code{"Cy5"}, storing intensities of probes on the green and red channels, respectively. This slot's value is
#'        \code{NULL} if no control probe information is available.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link[=samples,RnBSet-method]{samples}}}{Gets the identifiers of all samples in the dataset.}
#'   \item{\code{\link[=pheno,RnBSet-method]{pheno}}}{Gets the phenotypic and processing data of the dataset.}
#'   \item{\code{\link[=meth,RnBSet-method]{meth}}}{Gets the \code{matrix} of methylation beta-values of the dataset.}
#'   \item{\code{\link[=dpval,RnBeadSet-method]{dpval}}}{Gets the \code{matrix} of detection p-values of the dataset.}
#'   \item{\code{\link[=covg,RnBSet-method]{covg}}}{Gets the \code{matrix} of bead counts of the dataset.}
#'   \item{\code{\link[=qc,RnBeadSet-method]{qc}}}{Gets the intensities of the quality control probes.}
#'   \item{\code{\link[=remove.sites,RnBSet-method]{remove.sites}}}{Removes probes from the dataset.}
#'   \item{\code{\link[=remove.samples,RnBeadSet-method]{remove.samples}}}{Removes samples from the dataset.}
#'   \item{\code{\link[BiocGenerics]{combine}}}{Combines two datasets.}
#' }
#'
#' @name RnBeadSet-class
#' @rdname RnBeadSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBeadSet
setClass("RnBeadSet",
         representation(pval.sites="matrixOrffOrNULL",
				 		pval.regions="listOrNULL",
				 		qc="listOrNULL"#,
						#status="listOrNULL"
		#				bead.counts="listOrNULL"
		),
		 contains="RnBSet",
         prototype(#pheno=data.frame(),
				 #meth.sites=matrix(),
				 pval.sites=NULL,
				 pval.regions=NULL,
				 qc=NULL#, 
				 #status=NULL
#				 bead.counts=NULL
 		 ), package = "RnBeads"
)

## ---------------------------------------------------------------------------------------------------------------------
## VALIDITY
## ---------------------------------------------------------------------------------------------------------------------

##
## validRnBeadSetObject
##
## Valididty check for RnBeadSet-class
##
validRnBeadSetObject<-function(object){

  if(is.null(object@pheno) || is.null(object@meth.sites)) return("Pheno and betas slots can not be NULL")
  if(dim(object@pheno)[1L]!=dim(object@meth.sites)[2L]) return("The number of samples is ambiguous")

#  if(!is.null(object@pval.sites)){
#    if(!setequal(rownames(object@pval.sites), rownames(object@meth.sites))) return("Rows of the beta and p-value tables do not match")
#    if(!setequal(colnames(object@pval.sites), colnames(object@meth.sites))) return("Columns of the beta and p-value tables do not match")
#
#  }
#
#  if(!is.null(object@covg.sites)){
#    if(!setequal(rownames(object@covg.sites), rownames(object@meth.sites))) return("Rows of the beta and bead count tables do not match")
#    if(!setequal(colnames(object@covg.sites), colnames(object@meth.sites))) return("Columns of the beta and bead count tables do not match")
#  }

  return(TRUE)
}

setValidity("RnBeadSet", method=validRnBeadSetObject)

## ---------------------------------------------------------------------------------------------------------------------
## CONSTRUCTORS
## ---------------------------------------------------------------------------------------------------------------------

setMethod("initialize", "RnBeadSet",
		function(.Object,
				pheno=data.frame(),
				sites=matrix(nrow=0, ncol=0),
				meth.sites=matrix(nrow=0, ncol=0),
				pval.sites=NULL,
				covg.sites=NULL,
				qc=NULL,
				target="probes450",
				status=list(normalized=FALSE, background=FALSE, disk.dump=FALSE)
		) {
			
#			.Object@target<-target
#			
#			.Object@pheno<-pheno
#			.Object@sites<-sites
#			.Object@meth.sites<-betas
			.Object@pval.sites<-pval.sites
			#.Object@covg.sites<-bead.counts
			#.Object@regions<-list()
			
			status<-status
			
			#.Object@inferred.covariates <- list()
			
			.Object@qc<-qc
			
			callNextMethod(.Object,
					pheno=pheno,
					sites=sites,
					meth.sites=meth.sites,
					covg.sites=covg.sites,
					status=status,
					assembly=ifelse(target=="probesMMBC", "mm10", ifelse(target=="probesEPICv2", "hg38", "hg19")),
					target=target
			)
			
		})

########################################################################################################################

#' Wrapper function RnBeadSet
#'
#'
#' @param pheno       		Phenotypic data.
#' @param probes			\code{character} vector of Infinium(R) probe identifiers
#' @param betas      		\code{matrix} or \code{ff_matrix} of beta values. If \code{probes} are missing should contain Infinium probe identifiers as row names.
#' @param p.values    		\code{matrix} or \code{ff_matrix} of detection p-values.
#' @param bead.counts       ...
#' @param qc                ...
#' @param platform	   		\code{character} singleton specifying the microarray platform: \code{"450k"} corresponds to HumanMethylation450 microarray, and \code{"27k"} stands for HumanMethylation27.
#' @param summarize.regions ...
#' @param region.types		A \code{character} vector specifying the region types, for which the methylation infromation will be summarized.
#' @param useff		  		If \code{TRUE} the data matrices will be stored as \code{ff} objects
#' 
#' @return an object of class RnBeadSet
#' 
#' @name RnBeadSet
#' @rdname RnBeadSet-class
#' @aliases initialize,RnBeadSet-method
#' @export
RnBeadSet<-function(
		pheno,
		probes,
		betas,
		p.values = NULL, 
		bead.counts = NULL,
		qc = NULL,
		platform = "450k",
		summarize.regions = TRUE,
		region.types = rnb.region.types.for.analysis(ifelse(platform=="MMBC", "mm10", ifelse(target=="probesEPICv2", "hg38", "hg19"))),
		useff=rnb.getOption("disk.dump.big.matrices")
		){
		
		if(missing(pheno)){
			stop("argument pheno should be supplied")
		}
		
		if(missing(betas)){
			#warning("Valid RnBeadSet object was not created: pheno and/or betas missing")
			stop("argument betas should be supplied")
		}
		
		if(missing(probes)){
			if(!is.null(rownames(betas))){
				probes<-rownames(betas)
			}else{
				stop("If probes are not supplied, betas should have probe identifers as row names")
			}
		}
		
		if(!is.data.frame(pheno)){
			stop("invalid value for pheno: should be a data frame")
		}
		
		if(!any(c('matrix', 'ff_matrix') %in% class(betas))){
			stop("invalid value for betas: should be a matrix or an ff_matrix")
		}
	
		if(!is.null(p.values) && !any(c('matrix', 'ff_matrix') %in% class(p.values))){
			stop("invalid value for p.values: should be a matrix or an ff_matrix")
		}
		
		if(!is.null(bead.counts) && !any(c('matrix', 'ff_matrix') %in% class(bead.counts))){
			stop("invalid value for bead.counts: should be a matrix or an ff_matrix")
		}
		
		if(!is.null(qc)){
			if(!is.list(qc)){
				stop("invalid value for qc: should be a list")
			}
		}
		
		if(!is.character(platform) || length(platform)!=1L){
			stop("invalid value for platform: should be a character of length one")
		}
		
		if(!is.character(region.types)){
			stop("invalid value for region types: should be a character vector")
		}
		
		if(!is.logical(summarize.regions) || length(summarize.regions)!=1L){
			stop("invalid value for summarize.regions: should be a logical of length one")
		}
		
		if(!is.logical(useff) || length(useff)!=1L){
			stop("invalid value for useff: should be a logical of length one")
		}
		
		if (platform == "EPIC") {
			target <- "probesEPIC"
			assembly <- "hg19"
		}else if (platform == "EPICv2") {
			target <- "probesEPICv2"
			assembly <- "hg38"
		}else if (platform == "450k") {
			target <- "probes450"
			assembly <- "hg19"
		} else if(platform == "27k"){ 
			target <- "probes27"
			assembly <- "hg19"
		}else{
			stop("Invalid value for platform")
		}
		mr<-match.probes2annotation(probes, target, assembly)
		sites<-mr[[1]]
		site.ids<-mr[[2]]
		
		if(useff){
			if(!"ff" %in% class(betas)){
				betas<-convert.to.ff.matrix.tmp(betas[site.ids,,drop=FALSE])
			}
			if(!is.null(p.values) && !"ff" %in% class(p.values)){
				p.values<-convert.to.ff.matrix.tmp(p.values)
			}
			if(!is.null(bead.counts) && !"ff" %in% class(bead.counts)){
				bead.counts<-convert.to.ff.matrix.tmp(bead.counts)
			}
		}else{
			betas<-betas[site.ids,,drop=FALSE]
			if(!is.null(p.values)){
				p.values<-p.values[site.ids,,drop=FALSE]
			}
			if(!is.null(bead.counts)){
				bead.counts<-bead.counts[site.ids,,drop=FALSE]
			}
		}
		status<-list()
		status[["normalized"]]<-"none"
		status[["background"]]<-"none"
		status[["disk.dump"]]<-useff
		
		object<-new("RnBeadSet",
				target = target,
				pheno,
				sites=sites,
				meth.sites=betas,
				pval.sites=p.values,
				covg.sites=bead.counts,
				qc=qc,
				status=status
				)
				
		if(summarize.regions){
			for (region.type in region.types) {
				if (region.type %in% rnb.region.types("hg19")) {
						object <- summarize.regions(object, region.type)
					}
				}
			}
			
		return(object)
}

########################################################################################################################

## Prints a summary of a RnBeadSet dataset. This is used in the \code{show} methods.
## 
## @param object The RnBeadSet object to be presented.
##
rnb.show.rnbeadset <- function(object) {
	cat("Object of class ", class(object), "\n", sep = "")
	cat(sprintf("%8d samples\n", nrow(object@pheno)))
	cat(sprintf("%8d probes\n", nrow(object@sites)))
	probe.types <- rownames(object@sites)
	if (!is.null(probe.types)) {
		probe.types <- sapply(c("^cg", "^ch", "^rs"), function(type) { sum(grepl(type, probe.types)) })
		cat(sprintf("\tof which: %g CpG, %g CpH, and %g rs\n", probe.types[1], probe.types[2], probe.types[3]))
	}
	if (!is.null(object@regions)) {
		cat("Region types:\n")
		for (rn in names(object@regions)) {
			cat(sprintf("\t%8d regions of type %s\n", nrow(object@regions[[rn]]), rn))
		}
	}
	cat(sprintf("Intensity information is %s\n", ifelse(inherits(object, "RnBeadRawSet"), "present", "absent")))
	cat(sprintf("Detection p-values are %s\n", ifelse(is.null(object@pval.sites), "absent", "present")))
	cat(sprintf("Bead counts are %s\n", ifelse(is.null(object@covg.sites), "absent", "present")))
	cat(sprintf("Quality control information is %s\n", ifelse(is.null(qc(object)), "absent", "present")))
	cat(sprintf("Summary of normalization procedures:\n"))
	
	if(!is.null(object@status) && !is.null(object@status$normalized)){
		if(object@status$normalized=="none"){
			cat(sprintf("\tThe methylation data was not normalized.\n"))
		}else{
			cat(sprintf("\tThe methylation data was normalized with method %s.\n", object@status$normalized))			
		}
	}
	
	if(!is.null(object@status) && !is.null(object@status$background)){
		if(object@status$background=="none"){
			cat(sprintf("\tNo background correction was performed.\n"))
		}else{
			cat(sprintf("\tBackground correction was performed with method %s.\n", object@status$background))			
		}
	}
	if(isImputed(object)){
	  cat("\tData set was imputed.\n")
	}
}

setMethod("show", "RnBeadSet", rnb.show.rnbeadset)

########################################################################################################################
## MethyLumi2RnBeadSet
##
## Extracts all methylation data from a MethyLumiSet object, necessary to instantiate RnBeadSet
##
## @param methySet \code{MethylumiSet} instance to be converted.
## @return \code{list} with up to 5 elements: \code{"pheno"}, \code{"betas"}, \code{"bead.counts"}, \code{"p.values"},
##         \code{"qc"}
##
## @author Pavlo Lutsik
MethyLumiSet2RnBeadSet <- function(methySet) {
	result <- list()
	result[["pheno"]] <- as(phenoData(methySet), "data.frame")
	beta.vals <- betas(methySet)

	if(all(c("methylated.N", "unmethylated.N") %in% assayDataElementNames(methySet))){

		## Set beta values with low coverage (bead.counts) to NA
		result[["bead.counts"]]<-methylated.N(methySet)+unmethylated.N(methySet)
		beta.vals[methylated.N(methySet)<2]<-NA
		beta.vals[unmethylated.N(methySet)<2]<-NA
		
	}
	result[["probes"]] <- featureNames(methySet)
	result[["betas"]] <- beta.vals
	result[["p.values"]] <- pvals(methySet)

	if(!is.null(controlData(methySet))){

		## Extract quality control data
		cd<-controlData(methySet)
		green<-intensitiesByChannel(cd,"Cy3")
		red<-intensitiesByChannel(cd,"Cy5")

		if(length(grep('[A-z]', rownames(green)))==dim(green)[1L]){

			if("ProbeID" %in% varLabels(featureData(cd)))
				probe.id.col<-"ProbeID" else if ("Address" %in% varLabels(featureData(cd))) probe.id.col<-"Address"

			probe.mapping<-as(featureData(cd)[,probe.id.col], "data.frame")
			rownames(green)<-probe.mapping[rownames(green), probe.id.col]
			rownames(red)<-probe.mapping[rownames(red), probe.id.col]
			
		}

		result[["qc"]]<-list(Cy3=green, Cy5=red)
	}
	result
}

#' as("RnBeadSet", "MethyLumiSet")
#'
#' Convert a \code{\linkS4class{RnBeadSet}} object to \code{\linkS4class{MethyLumiSet}}
#' 
#' @name coercion-methods
#' 
setAs("MethyLumiSet", "RnBeadSet",

	function(from, to){

		if(!inherits(from,"MethyLumiSet")){
			stop("not a MethyLumiSet object:", deparse(substitute(methylumi.set)))
		}

		m.data <- MethyLumiSet2RnBeadSet(from)
		object<-RnBeadSet(
				pheno=m.data$pheno,
				probes=m.data$probes,
				betas=m.data$betas,
				p.values=m.data$p.values,
				bead.counts=m.data$bead.counts)
		if ("qc" %in% names(m.data)) {
			qc(object)<-m.data[["qc"]]
		}

		object
})

## ---------------------------------------------------------------------------------------------------------------------
## GETTERS
## ---------------------------------------------------------------------------------------------------------------------

if(!isGeneric("dpval")) setGeneric('dpval',
           function(object, ...) standardGeneric('dpval'))

#' dpval-methods
#'
#' Extract detection p-values from an object of \code{\linkS4class{RnBeadSet}} class.
#'
#' @param object 		\code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBeadRawSet}} object
#' @param type 			\code{character} singleton. If \code{sites} detection p-values per each available 
#' 						site is returned. Otherwise should be one of region types for for which the summarized 
#' 						p-values are available
#' @param row.names	    Flag indicating of row names are to be generated in the result.
#' @param i     	    Indices of sites/regions to be retrieved. By default (\code{NULL}), all will be retrieved.
#' @param j     	    Indices of samples to be retrieved. By default (\code{NULL}), all will be retrieved.
#' 
#' @return detection p-values available for the dataset in the form of a \code{matrix}.
#'
#' @rdname dpval-methods
#' @docType methods
#' @export
#' @aliases dpval
#' @aliases dpval,RnBeadSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' dp<-dpval(rnb.set.example, row.names=TRUE)
#' head(dp)
#' }
setMethod("dpval", signature(object = "RnBeadSet"),
	function(object, type = "sites", row.names = FALSE, i = NULL, j = NULL) {
		get.dataset.matrix(object, type, row.names, object@pval.sites, object@meth.regions, i, j)
	}
)

########################################################################################################################
  
setGeneric("qc", function(object) standardGeneric("qc"))

#' qc-methods
#'
#' Extracts HumanMethylation quality control information
#'
#' @param object Dataset of interest.
#' @return Quality control information available for the dataset in the form of a \code{list} with two elements:
#' \code{Cy3} and \code{Cy5}.
#'
#' @rdname qc-methods
#' @aliases qc
#' @aliases qc,RnBeadSet-method
#' @docType methods
#' @export 
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' qcinf<-dpval(rnb.set.example, row.names=TRUE)
#' head(qcinf$Cy3)
#' head(qcinf$Cy5)
#' }
setMethod("qc", signature(object="RnBeadSet"), function(object){ object@qc })

setGeneric("qc<-", function(object, value) standardGeneric("qc<-"))
setMethod("qc<-", signature(object = "RnBeadSet", value = "listOrNULL"),
	function(object, value) {
		## TODO: Add validation and export, or simply use @qc in our code
		object@qc <- value
		object
	}
)

## ---------------------------------------------------------------------------------------------------------------------
## MODIFIERS
## ---------------------------------------------------------------------------------------------------------------------

if (!isGeneric("remove.sites")) {
	setGeneric("remove.sites", function(object, probelist, verbose = TRUE) standardGeneric("remove.sites"))
}

#' @rdname remove.sites-methods
#' @aliases remove.sites,RnBeadSet-method
#' @docType methods
#' @export
setMethod("remove.sites", signature(object = "RnBeadSet"),
	function(object, probelist, verbose = TRUE) {
		inds <- get.i.vector(probelist, rownames(object@sites))
		if (length(inds) != 0) {
			if (!is.null(object@pval.sites)) {
				if(object@status$disk.dump){
					new.matrix<-object@pval.sites[-inds,,drop=FALSE]
					if(isTRUE(object@status$discard.ff.matrices)){
						delete(object@pval.sites)
					}
					object@pval.sites <- convert.to.ff.matrix.tmp(new.matrix)
					rm(new.matrix); rnb.cleanMem()
				}else{
					object@pval.sites <- object@pval.sites[-inds, ,drop=FALSE]
				}
			}
		}

		callNextMethod()
	}
)

########################################################################################################################

if (!isGeneric("remove.samples")) {
	setGeneric("remove.samples", function(object, samplelist) standardGeneric("remove.samples"))
}

#' @rdname remove.samples-methods
#' @aliases remove.samples,RnBeadSet-method
#' @docType methods
#' @export
setMethod("remove.samples", signature(object = "RnBeadSet"),
	function(object, samplelist) {
		inds <- get.i.vector(samplelist, samples(object))
		if (length(inds) != 0) {
			if (!is.null(object@pval.sites)) {
				if(object@status$disk.dump){
					new.matrix<-object@pval.sites[,-inds, drop=FALSE]
					if(isTRUE(object@status$discard.ff.matrices)){
						delete(object@pval.sites)
					}
					object@pval.sites <- convert.to.ff.matrix.tmp(new.matrix)
					rm(new.matrix); rnb.cleanMem()
				}else{
					object@pval.sites <- object@pval.sites[,-inds, drop=FALSE]
				}
			}
			if (!is.null(object@qc)) {
				object@qc$Cy3 <- object@qc$Cy3[,-inds, drop=FALSE]
				object@qc$Cy5 <- object@qc$Cy5[,-inds, drop=FALSE]
			}
		}
		callNextMethod()
	}
)

########################################################################################################################

setMethod("save.matrices", signature(object="RnBeadSet", path="character"),
		function(object, path){
			
			if(!is.null(object@pval.sites)){
				
				if(!is.null(object@status) && object@status$disk.dump){
					
					if("ff" %in% class(object@pval.sites)){
						ffmatrix<-object@pval.sites
						ffsave(ffmatrix, file=file.path(path, "rnb.pvals"),rootpath=getOption('fftempdir'))
						rm(ffmatrix)
						
					}
				}

				if(!is.null(object@regions)){
					
					for(rgn in 1:length(object@regions)){
						
						rgnpath<-file.path(path,rgn)
						if(!file.exists(rgnpath)){
							dir.create(rgnpath)	
						}

						if("ff" %in% class(object@pval.regions[[rgn]])){
							
							ffmatrix<-object@meth.regions[[rgn]]
							ffsave(ffmatrix, file=file.path(path, rgn, "rnb.pval"),rootpath=getOption('fftempdir'))
							rm(ffmatrix)
							
						}			
					}
				}
			}
			callNextMethod(object, path)
			
		})

########################################################################################################################
		
setMethod("load.matrices", signature(object="RnBeadSet", path="character"),
		
		function(object, path, temp.dir=tempdir()){
			
		if(!is.null(object@pval.sites)){
			
			if(sum(grepl("rnb.pvals",list.files(path)))==2){
				load_env<-new.env()
				suppressMessages(ffload(file=file.path(path, "rnb.pvals"), envir=load_env,rootpath=getOption("fftempdir")))
				object@pval.sites<-get("ffmatrix", envir=load_env)
				rm(load_env)
			}
			
			rgns <- object@regions
			if(!is.null(rgns)){
				if (.hasSlot(object, 'version')) {
					rngs <- 1:length(rgns)
				}
				for(rgn in rgns){
					if(sum(grepl("rnb.pvals",list.files(file.path(path, rgn))))==2){
						load_env<-new.env()
						suppressMessages(ffload(file=file.path(path, rgn, "rnb.pvals"), envir=load_env,rootpath=getOption("fftempdir")))
						object@pval.regions[[rgn]]<-get("ffmatrix", envir=load_env)
						rm(load_env)
					}		
				}
			}
		}
			
			callNextMethod(object=object, path=path, temp.dir=temp.dir)
			
		})

########################################################################################################################

#' @rdname destroy-methods
#' @aliases destroy,RnBeadSet-method
#' @docType methods
#' @export
setMethod("destroy", signature(object="RnBeadSet"),
		function(object){
			
			if(object@status$disk.dump){
				
				if(!is.null(object@pval.sites)){
					delete(object@pval.sites)
				}
			}
			callNextMethod()
			
		}
)

## ---------------------------------------------------------------------------------------------------------------------
## HELPER ROUTINES
## ---------------------------------------------------------------------------------------------------------------------
get.relative.covg<-function(m, design.vector){
	
	des1.quant<-sort(m[design.vector=="I",])[round(length(which(design.vector=="I"))*ncol(m)*0.001)]
	des2.quant<-sort(m[design.vector=="II",])[round(length(which(design.vector=="II"))*ncol(m)*0.001)]
	m[design.vector=="I",]<-m[design.vector=="I",]/as.double(des1.quant)
	m[design.vector=="II",]<-m[design.vector=="II",]/as.double(des2.quant)
	m
	
}
########################################################################################################################
summarize.bead.counts<-function(bead.counts.M, bead.counts.U, method="min"){
	
	nrow.M<-nrow(bead.counts.M); ncol.M<-ncol(bead.counts.M)
	nrow.U<-nrow(bead.counts.U); ncol.U<-ncol(bead.counts.U)
	
	if(nrow.M!=nrow.U || ncol.M!=ncol.U){
		stop("Dimensions of bead count matrices differ")
	}else{
		bead.counts<-matrix(NA,nrow.M,ncol.M)
		rownames(bead.counts)<-rownames(bead.counts.M)
	}
	
	if(method=="min"){
		index1<-bead.counts.M<=bead.counts.U
		index1[is.na(index1)]<-FALSE
		bead.counts[index1]<-bead.counts.M[index1]
		index2<-bead.counts.U<bead.counts.M
		index2[is.na(index2)]<-FALSE
		bead.counts[index2]<-bead.counts.U[index2]
		
	}
	
	bead.counts
	
}
########################################################################################################################
combine.sites<-function(sites1, sites2){
	
	common.chr<-intersect(unique(sites1[,2]), unique(sites2[,2]))
	
	common.sites<-lapply(common.chr, function(chr){
				
				sts<-intersect(sites1[sites1[,2]==chr,3],sites2[sites2[,2]==chr,3])				
				
				rbind(rep(1,length(sts)), rep(chr,length(sts)), sts)
				
			})
	
		
	if("ff_matrix" %in% c(class(sites1), class(sites2))){
	
		new.sites<-ff(vmode="integer", dim=c(sum(sapply(common.sites, nrow)),3))

		ixx<-1
		for(sts in common.sites){
			new.sites[,][ixx:(ixx+nrow(sts)),]<-sts
			ixx<-ixx+nrow(sts)+1
		}
		
	}else{
		new.sites<-do.call("rbind", common.sites)
	}
	
	
	new.sites
		
}
########################################################################################################################
match.probes2annotation<-function(probes, target="probes450", assembly="hg19"){
	
	if(!is.character(probes)){
		stop("wrong value for probes")
	}

	probe.identifiers.pattern <- "^cg|^ch|^rs"
	## EPIC v2 has a new probe type "nv"
	if (target == "probesEPICv2") {
		probe.identifiers.pattern <- paste(probe.identifiers.pattern, "|^nv", sep="")
	}
	
	if(!all(grepl(probe.identifiers.pattern, probes))){
		stop("probes contains invalid IDs")
	}
	
	dupl<-duplicated(probes)
	if(any(dupl)){
		stop(sprintf("some supplied probe IDs are duplicated: e.g. %s", paste(probes[which(dupl)[1:3]],collapse=", ")))
	}
		
	## Load probe annotation table
	probe.annotation <- rnb.get.annotation(target, assembly)
    ### TODO: remove after fixing the annotation
    if(target=="probesMMBC"){
        probe.annotation<-endoapply(probe.annotation, function(ag) {
                    names(ag)<-mcols(ag)[["ID"]]
                    mcols(ag)[["Context"]]<-factor(c("cg"="CG", "ch"="CH", "rs"="Other")[substr(names(ag), 1,2)], levels=c("CG", "CH", "Other"))
                    return(ag)})
    }
	## Construct the matrix of site indices
	#x.data<-GenomicRanges::as.data.frame(probe.annotation)
	
	chrs<-rnb.get.chromosomes(assembly)
	
	annotated.probes<-do.call("c", lapply(probe.annotation, names))
	if(length(which(probes %in% annotated.probes))<2){
		err<-"Annotations could be found for less than two rows from the supplied beta value table"
		rnb.error(err)		
	}
	
	if(!all(probes %in% annotated.probes)){
		warn<-"Some of the supplied probes are missing annotation and will be discarded"
		rnb.warning(warn)			
	}
	
	x.data<-rnb.annotation2data.frame(probe.annotation)
	rownames(x.data)<-annotated.probes
	rm(annotated.probes)
	
	#site.ids<-character()
	site.ids<-list()
	p.infos <- lapply(unique(x.data[["Chromosome"]]), 
			function(chr) {
				chr.map<-which(x.data[["Chromosome"]]==chr)
				chr.ids<-x.data[chr.map, "ID"]
				table<-match(probes, chr.ids)
				present<-!is.na(table)
				po<-order(table[present])
				site.ids[[chr]]<<-which(present)[po]
				table<-table[present][po]
				matrix(c(
						as.integer(x.data[chr.map[table], "Context"]), 
						match(x.data[chr.map[table],"Chromosome"], names(chrs)), 
						table),
						ncol = 3, 
						dimnames = list(x.data[chr.map[table],"ID"], c("context", "chr", "index")))
			})
	
#	p.infos<-matrix(c(as.integer(x.data[site.ids, "Context"]), match(x.data[site.ids,"Chromosome"], names(chrs)), unlist(p.indices)),
#			ncol = 3, dimnames = list(site.ids, c("context", "chr", "index")))
	
	p.infos<-do.call("rbind", p.infos)
	site.ids<-do.call("c", site.ids)
			
	rownames(p.infos)<-probes[site.ids]

	fm <- FALSE
	if (length(site.ids)==nrow(x.data)){
		fm <- all(site.ids==1:nrow(x.data))
	}
	
	return(list(probe.infos=p.infos, ids=site.ids, full.match=fm))
	
}
