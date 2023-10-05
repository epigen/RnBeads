########################################################################################################################
## RnBeadRawSet-class.R
## created: 2013-xx-xx
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBeadRawSet class definition.
########################################################################################################################

## ---------------------------------------------------------------------------------------------------------------------
## CLASS DEFINITIONS
## ---------------------------------------------------------------------------------------------------------------------

#' RnBeadRawSet-class 
#' 
#' Main class for storing HumanMethylation micorarray data which includes intensity information 
#' 
#' @section Slots:
#' \describe{
#'	\item{\code{pheno}}{Phenotypic data.}
#'	\item{\code{M}}{\code{matrix} of intensities for the probes measuring the abundance of methylated molecules.} 
#'	\item{\code{U}}{\code{matrix} of intensities for the probes measuring the abundance of unmethylated molecules.} 
#' 	\item{\code{M0}}{\code{matrix} of "out-of-band" intensities for the probes measuring the abundance of methylated molecules.} 
#' 	\item{\code{U0}}{\code{matrix} of "out-of-band" intensities for the probes measuring the abundance of unmethylated molecules.}
#' 	\item{\code{bead.counts.M}}{\code{matrix} of bead counts per probe.}
#' 	\item{\code{bead.counts.U}}{\code{matrix} of bead counts per probe.}
#' }
#' 
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{samples}}{Gets the identifiers of all samples in the dataset.}
#'   \item{\code{\link[=M,RnBeadRawSet-method]{M}}}{Get the matrix of intensities for the probes measuring the abundance of methylated molecules.}
#'   \item{\code{\link[=U,RnBeadRawSet-method]{U}}}{Get the matrix of intensities for the probes measuring the abundance of unmethylated molecules.}
#'   \item{\code{\link{intensities.by.color}}}{Get probe intensities in each color channel.}
#' } 
#' 
#' @name RnBeadRawSet-class
#' @rdname RnBeadRawSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBeadRawSet
#' @include RnBeadSet-class.R
setClass("RnBeadRawSet",
		representation(M="matrixOrffOrNULL",
				U="matrixOrffOrNULL",
				M0="matrixOrffOrNULL",
				U0="matrixOrffOrNULL",
				bead.counts.M="matrixOrffOrNULL",
				bead.counts.U="matrixOrffOrNULL"
		),
		contains="RnBeadSet",
		prototype(#pheno=data.frame(),
				#betas=matrix(),
				#meth.sites=matrix(),
				#pval.sites=NULL,
				#pval.regions=NULL,
				#qc=NULL, 
				#status=NULL,
				M=matrix(),
				U=matrix(),
				M0=NULL,
				U0=NULL,
				bead.counts.M=NULL,
				bead.counts.U=NULL
		),
		package = "RnBeads"
)

RNBRAWSET.SLOTNAMES<-c("M","U","M0","U0","bead.counts.M", "bead.counts.U")

########################################################################################################################

## initialize.RnBeadRawSet
##
## Direct slot filling.
## 
## #docType methods
## #rdname RnBeadRawSet-class
setMethod("initialize", "RnBeadRawSet",
		function(.Object,
				pheno = data.frame(),
				sites = matrix(ncol=0, nrow=0),
				meth.sites = matrix(ncol=0, nrow=0),
				M = matrix(ncol=0, nrow=0),
				U = matrix(ncol=0, nrow=0),
				M0 = NULL,
				U0 = NULL,
				bead.counts.M = NULL, 
				bead.counts.U = NULL,
				covg.sites = NULL,
				pval.sites = NULL,
				qc = NULL,
				target="probes450",
				status=list(normalized=FALSE, background=FALSE, disk.dump=FALSE)
				) {
			
			.Object@target<-target
					
			#.Object@pheno<-pheno
			#.Object@sites<-sites
			.Object@M<-M
			.Object@U<-U
			.Object@M0<-M0
			.Object@U0<-U0
			.Object@bead.counts.M<-bead.counts.M
			.Object@bead.counts.U<-bead.counts.U
			
			#.Object@meth.sites<-betas
			#.Object@pval.sites<-p.values
			#.Object@covg.sites<-bead.counts
			#.Object@regions<-list()
			
			#.Object@status<-list()
			
			#.Object@status[["normalized"]]<-"none"
			#.Object@status[["background"]]<-"none"
			#.Object@status[["disk.dump"]]<-useff
			
			#.Object@inferred.covariates <- list()
			
			#.Object@qc<-qc
			
			#.Object
			callNextMethod(.Object,
					pheno=pheno,
					sites=sites,
					meth.sites=meth.sites, 
					pval.sites=pval.sites,
					covg.sites=covg.sites,
					target=target,
					status=status,
					qc=qc)
		})

########################################################################################################################		

#' Wrapper function RnBeadRawSet
#'
#' @param pheno       		Phenotypic data.
#' @param probes			\code{character} vector of Infinium(R) probe identifiers
#' @param M       	  		Matrix of intensities for the probes measuring the abundance of methylated molecules 
#' @param U       	  		Matrix of intensities for the probes measuring the abundance of unmethylated molecules 
#' @param M0       	  		Matrix of "out-of-band" intensities for the probes measuring the abundance of methylated molecules 
#' @param U0       	  		Matrix of "out-of-band" intensities for the probes measuring the abundance of unmethylated molecules
#' @param bead.counts.M 	Matrix of bead counts per probe.
#' @param bead.counts.U 	Matrix of bead counts per probe.
#' @param p.values    		Matrix of detection p-values.
#' @param qc                ...
#' @param platform	   		\code{character} singleton specifying the microarray platform: \code{"450k"} corresponds to HumanMethylation450 microarray, and \code{"27k"} stands for HumanMethylation27.
#' @param region.types		A \code{character} vector specifying the region types, for which the methylation infromation will be summarized.
#' @param beta.offset		A regularization constant which is added to the denominator at beta-value calculation
#' @param summarize.bead.counts	If \code{TRUE} the coverage slot is filled by summarizing the \code{bead.counts.M} and \code{bead.counts.U} matrices. For type I probes the summarization is done using \code{min} operation, while for type II probes the bead counts should be identical in both supplied matrices
#' @param summarize.regions ...
#' @param useff		  		If \code{TRUE} the data matrices will be stored as \code{ff} objects
#' @param ffcleanup		  	If \code{TRUE} and disk dumping has been enabled the data of the input \code{ff} objects will be deleted 
#'
#' @return an object of class RnBeadRawSet
#' 
#' @name RnBeadRawSet
#' @rdname RnBeadRawSet-class
#' @aliases initialize,RnBeadRawSet-method
#' @export
RnBeadRawSet<-function(
		pheno,
		probes,
		M,
		U,
		M0=NULL,
		U0=NULL,
		bead.counts.M = NULL, 
		bead.counts.U = NULL,
		p.values=NULL,
		qc = NULL,
		platform = "450k",
		beta.offset=100,
		summarize.bead.counts=TRUE,
		summarize.regions=TRUE,
		region.types = rnb.region.types.for.analysis(ifelse(platform=="MMBC", "mm10", ifelse(target=="probesEPICv2", "hg38", "hg19"))),
		useff=rnb.getOption("disk.dump.big.matrices"),
		ffcleanup=FALSE){
		
		if(missing(pheno)){
			stop("argument pheno should be supplied")
		}
		
		if(missing(probes)){
			if(!is.null(rownames(M))){
				probes<-rownames(M)
			}else if(!is.null(rownames(U))){
				probes<-rownames(U)
			}else{
				stop("If probes are not supplied, betas should have probe identifers as row names")
			}
		}
		
		if(!is.data.frame(pheno)){
			stop("invalid value for pheno: should be a data frame")
		}
		
		if(!any(c('matrix', 'ff_matrix') %in% class(M))){
			stop("invalid value for M: should be a matrix or an ff_matrix")
		}
		
		if(!any(c('matrix', 'ff_matrix') %in% class(U))){
			stop("invalid value for U: should be a matrix or an ff_matrix")
		}
		
		if(!is.null(p.values) && !any(c('matrix', 'ff_matrix') %in% class(p.values))){
			stop("invalid value for p.values: should be a matrix or an ff_matrix")
		}
		
		if(!is.null(bead.counts.M) && !any(c('matrix', 'ff_matrix') %in% class(bead.counts.M))){
			stop("invalid value for bead.counts.M: should be a matrix or an ff_matrix")
		}
		
		if(!is.null(bead.counts.U) && !any(c('matrix', 'ff_matrix') %in% class(bead.counts.U))){
			stop("invalid value for bead.counts.U: should be a matrix or an ff_matrix")
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
		
		if(!is.numeric(beta.offset) || length(beta.offset)!=1L || beta.offset<0){
			stop("invalid value for beta.offset: should be a positive numeric of length one")
		}
		
		if(!is.logical(summarize.regions) || length(summarize.regions)!=1L){
			stop("invalid value for summarize.regions: should be a logical of length one")
		}
		
		if(!is.logical(summarize.bead.counts) || length(summarize.bead.counts)!=1L){
			stop("invalid value for summarize.bead.counts: should be a logical of length one")
		}
		
		if(!is.logical(useff) || length(useff)!=1L){
			stop("invalid value for useff: should be a logical of length one")
		}
		
		if (platform =="EPIC") {
			target <- "probesEPIC"
			assembly <- "hg19" ## TODO: EPICv1 will be hg38 compatible
		}else if (platform =="EPICv2") {
			target <- "probesEPICv2"
			assembly <- "hg38"
		}else if (platform =="450k") {
			target <- "probes450"
			assembly <- "hg19"
		} else if(platform == "27k"){
			target <- "probes27"
			assembly <- "hg19"
        }else if(platform == "MMBC"){
            target <- "probesMMBC"
            assembly <- "mm10"
		}else{
			rnb.error("Invalid value for platform")
		}
		

		res<-match.probes2annotation(probes, target, assembly)
		sites<-res[[1]]
		site.ids<-res[[2]]
		fullmatch<-res[[3]]
		
		M<-prepare.slot.matrix(M, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		
		U<-prepare.slot.matrix(U, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		
		if(!is.null(M0)){
			M0<-prepare.slot.matrix(M0, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		}
		if(!is.null(U0)){
			U0<-prepare.slot.matrix(U0, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		}
		if(!is.null(bead.counts.M)){
			bead.counts.M<-prepare.slot.matrix(bead.counts.M, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		}
		if(!is.null(bead.counts.U)){
			bead.counts.U<-prepare.slot.matrix(bead.counts.U, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		}
		if(!is.null(p.values)){
			p.values<-prepare.slot.matrix(p.values, useff=useff, full.match=fullmatch,  subset=site.ids, cleanup=ffcleanup)
		}
		
		rownames(M)<-NULL			
		rownames(U)<-NULL
		
		if(!is.null(M0)){
			rownames(M0)<-NULL
		}
		if(!is.null(U0)){
			rownames(U0)<-NULL
		}
		if(!is.null(bead.counts.M)){
			rownames(bead.counts.M)<-NULL
			M[,][bead.counts.M[,]<1]<-NA
			if(!is.null(M0)){
				M0[,][bead.counts.M[,]<1]<-NA
			}
		}
		if(!is.null(bead.counts.U)){
			rownames(bead.counts.U)<-NULL
			U[,][bead.counts.U[,]<1]<-NA
			if(!is.null(U0)){
				U0[,][bead.counts.U[,]<1]<-NA
			}
		}
		
		if(!is.null(bead.counts.M) && !is.null(bead.counts.U) && summarize.bead.counts){
			bead.counts<-summarize.bead.counts(bead.counts.M[,,drop=FALSE],bead.counts.U[,,drop=FALSE])
		}else{
			bead.counts<-NULL
		}
		
		betas<-beta.value(M[,,drop=FALSE],U[,,drop=FALSE], beta.offset)
		
		if(useff){
			if(!is.null(bead.counts)){
				bead.counts<-convert.to.ff.matrix.tmp(bead.counts)
			}
		}
		
		if(useff){
			betas<-convert.to.ff.matrix.tmp(betas)
		}
		
		status<-list()
		status[["normalized"]]<-"none"
		status[["background"]]<-"none"
		status[["disk.dump"]]<-useff
		
		object<-new("RnBeadRawSet",
			pheno=pheno,
			sites=sites,
			meth.sites=betas,
			M=M,
			U=U,
			M0=M0,
			U0=U0,
			bead.counts.M=bead.counts.M,
			bead.counts.U=bead.counts.U,
			covg.sites=bead.counts,
			pval.sites=p.values,
			qc=qc,
			target=target,
			status=status
		)
		
		if(summarize.regions){
			for (region.type in region.types) {
				if (region.type %in% rnb.region.types(assembly)) {
					object <- summarize.regions(object, region.type)
				}
			}
		}
		return(object)
}
########################################################################################################################
		
## FIXME: dummy validity method
validRnBeadRawSetObject<-function(object){
	return(TRUE)
}

setValidity("RnBeadRawSet", method=validRnBeadRawSetObject)

########################################################################################################################

setMethod("show", "RnBeadRawSet", rnb.show.rnbeadset)

########################################################################################################################

#' Conversion to/from RnBeadRawSet
#' 
#' The \code{"as"} method can be used for the following conversions:
#' \itemize{
#'   \item{}{\code{MethyLumiSet} (in package \pkg{methylumi}) to \code{\linkS4class{RnBeadRawSet}}}
#'   \item{}{\code{RnBeadRawSet} to \code{MethyLumiSet}}
#'   \item{}{\code{RGChannelSet} (in package \pkg{minfi}) to \code{\linkS4class{RnBeadRawSet}}}
#' }
#' 
#' @name as.RnBeadRawSet
setAs("MethyLumiSet", "RnBeadRawSet",
		
		function(from, to){
			
			if(!inherits(from,"MethyLumiSet")){
				stop("not a MethyLumiSet object:", deparse(substitute(methylumi.set)))
			}
			
			m.data <- MethyLumiSet2RnBeadSet(from)
			
			if("methylated.N" %in% ls(from@assayData) && "unmethylated.N" %in% ls(from@assayData) ){
				meth.N.element<-"methylated.N"
				umeth.N.element<-"unmethylated.N"
			}else if("Avg_NBEADS_A" %in% ls(from@assayData) && "Avg_NBEADS_B" %in% ls(from@assayData)){
				meth.N.element<-"Avg_NBEADS_B"	
				umeth.N.element<-"Avg_NBEADS_A"
			}else{
				meth.N.element<-NULL
				umeth.N.element<-NULL
			}
			
			if("methylated.OOB" %in% ls(from@assayData) && "unmethylated.OOB" %in% ls(from@assayData)){
				meth.oob.element<-"methylated.OOB"
				umeth.oob.element<-"unmethylated.OOB"
			}else{
				meth.oob.element<-NULL
				umeth.oob.element<-NULL
			}
			
			if(annotation(from)=="IlluminaMethylationEPIC"){
				platform="EPIC"
			}else if(annotation(from)=="IlluminaMethylationEPIC"){
				platform="EPICv2"
			}else if(annotation(from)=="IlluminaHumanMethylation450k"){
				platform="450k"
			}else if(annotation(from)=="IlluminaHumanMethylation27k"){
				platform="27k"
			}
			
			object<-RnBeadRawSet(
					pheno=m.data$pheno,
					probes=m.data$probes,
					M=methylated(from),
					U=unmethylated(from),
					p.values=m.data$p.values,
					M0=if(!is.null(meth.oob.element)) get(meth.oob.element,from@assayData) else NULL,
					U0=if(!is.null(umeth.oob.element)) get(umeth.oob.element,from@assayData) else NULL,
					bead.counts.M=if(!is.null(meth.N.element)) get(meth.N.element,from@assayData) else NULL,
					bead.counts.U=if(!is.null(umeth.N.element)) get(umeth.N.element,from@assayData) else NULL,
					platform=platform
					)
			
			if ("qc" %in% names(m.data)) {
				qc(object)<-m.data[["qc"]]
			}
						
			object
			
		})

########################################################################################################################
		
setAs("RnBeadRawSet","MethyLumiSet",
		
		function(from, to){
			
			if(!inherits(from,"RnBeadRawSet")){
				stop("not a RnBeadRawSet object:", deparse(substitute(methylumi.set)))
			}
				
			assd<-new.env()
			
			assign("betas", meth(from),  envir=assd)
			assign("pvals", dpval(from),  envir=assd)
			assign("methylated", M(from),  envir=assd)
			assign("unmethylated", U(from),  envir=assd)
			
			if(!is.null(M0(from))){
				assign("methylated.OOB", M0(from),  envir=assd)
			}
			if(!is.null(U0(from))){
				assign("unmethylated.OOB", U0(from),  envir=assd)
			}
			if(!is.null(bead.counts.M(from))){
				assign("methylated.N", bead.counts.M(from),  envir=assd)
			}
			if(!is.null(bead.counts.U(from))){
				assign("unmethylated.N", bead.counts.U(from),  envir=assd)
			}
			
			pd<-pheno(from)
			rownames(pd)<-colnames(meth(from))
			
			mset<-new(to, assd, as(pd, "AnnotatedDataFrame"))
			
			rm(assd)
			
			ann<-annotation(from, add.names=TRUE)
			featureData(mset)<-as(data.frame(COLOR_CHANNEL=ann$Color, rownames=rownames(ann)), "AnnotatedDataFrame")
			featureNames(mset)<-ann[["ID"]]
			
			if (!is.null(qc(from))){
				
				assd<-new.env()
				
				assign("methylated", qc(from)$Cy3,  envir=assd)
				assign("unmethylated", qc(from)$Cy5,  envir=assd)
				
				mset@QC<-new("MethyLumiQC", assd)
					
				if(from@target == "probesEPIC"){
					probeIDs<-rnb.get.annotation("controlsEPIC")[,"Target"]
					## TODO remove this after annotation has been fixed
					index<-rnb.update.controlsEPIC.enrich(rnb.get.annotation("controlsEPIC"))[,"Index"]
					probeIDs<-paste(probeIDs, index, sep=".")
				}else if(from@target == "probesEPICv2"){
					probeIDs<-rnb.get.annotation("controlsEPICv2")[,"Target"]
					probeIDs<-paste(probeIDs, unlist(sapply(table(probeIDs)[unique(probeIDs)], seq, from=1 )), sep=".")
				}else if(from@target == "probes450"){
					probeIDs<-rnb.get.annotation("controls450")[,"Target"]
					probeIDs<-paste(probeIDs, unlist(sapply(table(probeIDs)[unique(probeIDs)], seq, from=1 )), sep=".")
				}else if(from@target == "probes27"){
					probeIDs<-rnb.get.annotation("controls27")[,"Name"]
				}
				
				featureData(mset@QC)<-as(data.frame(Address=rownames(qc(from)$Cy3), rownames=probeIDs), "AnnotatedDataFrame")
				featureNames(mset@QC)<-probeIDs
								
				if(from@target == "probesEPIC"){
					annotation(mset@QC) <- "IlluminaMethylationEPIC"
				}else if(from@target == "probesEPICv2"){
					annotation(mset@QC) <- "IlluminaMethylationEPICv2"
				}else if(from@target == "probes450"){
					annotation(mset@QC) <- "IlluminaHumanMethylation450k"
				}else if(from@target == "probes27"){
					annotation(mset) <- "IlluminaHumanMethylation27k"
				}
			}
			
			if(from@target == "probesEPIC"){
				annotation(mset) <- "IlluminaMethylationEPIC"
			}else if(from@target == "probesEPICv2"){
				annotation(mset) <- "IlluminaMethylationEPICv2"
			}else if(from@target == "probes450"){
				annotation(mset) <- "IlluminaHumanMethylation450k"
			}else if(from@target == "probes27"){
				annotation(mset) <- "IlluminaHumanMethylation27k"
			}
			
			mset
			
		})

########################################################################################################################

setAs("RGChannelSet", "RnBeadRawSet", function(from, to) {

		## Get Illumina assay (platform)
		target.info <- getManifest(from)
		assay.name <- target.info@annotation
		if (!(is.character(assay.name) && length(assay.name) == 1 && isTRUE(assay.name != ""))) {
			stop("Unsupported platform; expected one-element character")
		}
		if (assay.name == "IlluminaHumanMethylationEPIC") {
			assay.name <- "probesEPIC"
			platform.name <- "EPIC"
		} else if (assay.name == "IlluminaHumanMethylationEPICv2") {
			assay.name <- "probesEPICv2"
			platform.name <- "EPICv2"
		} else if (assay.name == "IlluminaHumanMethylation450k") {
			assay.name <- "probes450"
			platform.name <- "450k"
		} else if (assay.name == "IlluminaHumanMethylation27k") {
			assay.name <- "probes27"
			platform.name <- "27k"
		} else {
			stop(paste("Unsupported platform", assay.name))
		}
	
		## Use RnBeads' mapping from probe IDs to addresses
		probes.all <- rnb.get.annotation(assay.name, ifelse(assay.name=="probesEPICv2", "hg38", "hg19"))
		probes.all <- lapply(probes.all, function(x) {
				result <- as.data.frame(mcols(x)[, c("Design", "Color", "AddressA", "AddressB")])
				rownames(result) <- names(x)
				result
			}
		)
		probes.all <- do.call(rbind, unname(probes.all))
		controls.all <- rnb.get.annotation(sub("^probes", "controls", assay.name), ifelse(assay.name=="probesEPICv2", "hg38", "hg19"))
		controls.all <- controls.all[, "ID"]

		## Extract data on signals
		mm.green <- getGreen(from)
		mm.red <- getRed(from)
		addresses.raw <- as.integer(rownames(mm.green))
		probes.supported <- probes.all[, 3] %in% addresses.raw | probes.all[, 4] %in% addresses.raw
		probes.all <- probes.all[probes.supported, ]
		rm(target.info, assay.name, probes.supported)

		## Initialize matrices of signal intensities and bead counts
		M <- matrix(NA_real_, nrow = nrow(probes.all), ncol = ncol(mm.red))
		U <- matrix(NA_real_, nrow = nrow(probes.all), ncol = ncol(mm.red))
		M0 <- matrix(NA_real_, nrow = nrow(probes.all), ncol = ncol(mm.red))
		U0 <- matrix(NA_real_, nrow = nrow(probes.all), ncol = ncol(mm.red))
		if (inherits(from, "RGChannelSetExtended")) {
			mm.beads <- get("NBeads", pos = from@assayData)
			beads.M <- matrix(0L, nrow = nrow(probes.all), ncol = ncol(mm.red))
			beads.U <- matrix(0L, nrow = nrow(probes.all), ncol = ncol(mm.red))
		} else {
			beads.M <- NULL
			beads.U <- NULL
		}

		## Map signals on addresses to probe signal intensities
		ii <- which(probes.all[, "Design"] == "II")
		jj <- findInterval(probes.all[ii, "AddressA"], addresses.raw)
		M[ii, ] <- mm.green[jj, ]
		U[ii, ] <- mm.red[jj, ]
		if (!is.null(beads.M)) {
			beads.M[ii, ] <- mm.beads[jj, ]
			beads.U[ii, ] <- mm.beads[jj, ]
		}
		ii <- which(probes.all[, "Design"] == "I" & probes.all[, "Color"] == "Grn")
		jj <- findInterval(probes.all[ii, "AddressB"], addresses.raw)
		M[ii, ] <- mm.green[jj, ]
		M0[ii, ] <- mm.red[jj, ]
		if (!is.null(beads.M)) {
			beads.M[ii, ] <- mm.beads[jj, ]
		}
		jj <- findInterval(probes.all[ii, "AddressA"], addresses.raw)
		U[ii, ] <- mm.green[jj, ]
		U0[ii, ] <- mm.red[jj, ]
		if (!is.null(beads.M)) {
			beads.U[ii, ] <- mm.beads[jj, ]
		}
		ii <- which(probes.all[, "Design"] == "I" & probes.all[, "Color"] == "Red")
		jj <- findInterval(probes.all[ii, "AddressB"], addresses.raw)
		M[ii, ] <- mm.red[jj, ]
		M0[ii, ] <- mm.green[jj, ]
		if (!is.null(beads.M)) {
			beads.M[ii, ] <- mm.beads[jj, ]
		}
		jj <- findInterval(probes.all[ii, "AddressA"], addresses.raw)
		U[ii, ] <- mm.red[jj, ]
		U0[ii, ] <- mm.green[jj, ]
		if (!is.null(beads.M)) {
			beads.U[ii, ] <- mm.beads[jj, ]
		}

		## TODO: Extract detection p-values if present

		## Extract control probe signals if available
		jj <- findInterval(controls.all, addresses.raw, all.inside = TRUE)
		ii <- which(controls.all == addresses.raw[jj])
		if (length(ii) == 0) {
			qc <- NULL
		} else {
			jj <- jj[ii]
			qc <- matrix(NA_real_, nrow = length(controls.all), ncol = ncol(mm.red))
			rownames(qc) <- as.character(controls.all)
			qc <- list("Cy3" = qc, "Cy5" = qc)
			qc[[1]][ii, ] <- mm.green[jj, ]
			qc[[2]][ii, ] <- mm.red[jj, ]
			if (all(sapply(qc, function(x) { all(is.na(x)) }))) {
				qc <- NULL
			}
		}
    p.data <- as.data.frame(pData(from))
		## Construct the resulting object
		RnBeadRawSet(p.data, rownames(probes.all), M, U, M0, U0, beads.M, beads.U, NULL, qc, platform.name)
	}
)

setAs("RnBeadRawSet", "RGChannelSet", function(from, to){
	assay.name <- from@target
#	probes.all <- rnb.get.annotation(assay.name, "hg19")
#	probes.all <- lapply(probes.all, function(x) {
#			result <- as.data.frame(mcols(x)[, c("Design", "Color", "AddressA", "AddressB")])
#			rownames(result) <- names(x)
#			result
#		}
#	)
#	probes.all <- do.call(rbind, unname(probes.all))
	probes.all <- annotation(from)[,c("Design", "Color", "AddressA", "AddressB")]
	controls.all <- rnb.get.annotation(sub("^probes", "controls", assay.name), ifelse(assay.name=="probesEPICv2", "hg38", "hg19"))
	controls.all <- controls.all[, "ID"]

	# Obtain methylated and unmethylated intensities
	M <- M(from)
	U <- U(from)
#	beads.M <- from@bead.counts.M
#	beads.U <- from@bead.counts.U
	M0 <- from@M0
	U0 <- from@U0

	#mm.green <- matrix(NA_real_, nrow=nrow(probes.all), ncol=ncol(M))
	#mm.red <- matrix(NA_real_, nrow=nrow(probes.all), ncol=ncol(M))

	is.type.II <- which(probes.all[,"Design"]%in%"II")
	addresses <- probes.all[is.type.II,"AddressA"]
	mm.green <- M[is.type.II,]
	mm.red <- U[is.type.II,]
#	mm.beads <- beads.M[is.type.II,,drop=F]
	row.names(mm.green) <- addresses
	row.names(mm.red) <- addresses

	is.type.I.green <- which(probes.all[, "Design"] == "I" & probes.all[, "Color"] == "Grn")
	addresses <- probes.all[is.type.I.green,"AddressB"]
	r.names <- c(row.names(mm.green),addresses)
	mm.green <- rbind(mm.green,M[is.type.I.green,])
	row.names(mm.green) <- r.names
	r.names <- c(row.names(mm.red),addresses)
	mm.red <- rbind(mm.red,M0[is.type.I.green,,drop=F])
	row.names(mm.red) <- r.names
#	mm.beads <- rbind(mm.beads,beads.M[is.type.I.green,,drop=F])
#	mm.beads <- rbind(mm.beads,beads.U[is.type.I.green,,drop=F])

	addresses <- probes.all[is.type.I.green,"AddressA"]
	r.names <- c(row.names(mm.green),addresses)
	mm.green <- rbind(mm.green,U[is.type.I.green,])
	row.names(mm.green) <- r.names
	r.names <- c(row.names(mm.red),addresses)
	mm.red <- rbind(mm.red,U0[is.type.I.green,,drop=F])
	row.names(mm.red) <- r.names

	is.type.I.red <- which(probes.all[, "Design"] == "I" & probes.all[, "Color"] == "Red")
	addresses <- probes.all[is.type.I.red,"AddressB"]
	r.names <- c(row.names(mm.red),addresses)
	mm.red <- rbind(mm.red,M[is.type.I.red,])
	row.names(mm.red) <- r.names
	r.names <- c(row.names(mm.green),addresses)
	mm.green <- rbind(mm.green,M0[is.type.I.red,,drop=F])
	row.names(mm.green) <- r.names
#	mm.beads <- rbind(mm.beads,beads.M[is.type.I.red,,drop=F])
#	mm.beads <- rbind(mm.beads,beads.U[is.type.I.red,,drop=F])

	addresses <- probes.all[is.type.I.red,"AddressA"]
	r.names <- c(row.names(mm.red),addresses)
	mm.red <- rbind(mm.red,U[is.type.I.red,])
	row.names(mm.red) <- r.names
	r.names <- c(row.names(mm.green),addresses)
	mm.green <- rbind(mm.green,U0[is.type.I.red,,drop=F])
	row.names(mm.green) <- r.names

	r.names <- c(row.names(mm.green),controls.all)
	mm.green <- rbind(mm.green,qc(from)[[1]])
	row.names(mm.green) <- r.names
	r.names <- c(row.names(mm.red),controls.all)
	mm.red <- rbind(mm.red,qc(from)[[2]])
	row.names(mm.red) <- r.names
#	mm.beads <- rbind(mm.beads,
#		matrix(NA,nrow=nrow(qc(from)[[1]]),ncol=ncol(mm.beads)))

	if(assay.name %in% "probesEPIC"){
		anno <- "IlluminaHumanMethylationEPIC"
	}else if(assay.name %in% "probesEPICv2"){
		anno <- "IlluminaHumanMethylationEPICv2"
	}else if(assay.name %in% "probes450"){
		anno <- "IlluminaHumanMethylation450k"
	}else if(assay.name %in% "probes27"){
		anno <- "IlluminaHumanMethylation27k"
	}else{
		stop("Unsupported platform")
	}
	p.data <- pheno(from)
	if(!is.null(rnb.getOption("identifiers.column"))){
		row.names(p.data) <- p.data[,rnb.getOption("identifiers.column")]
		colnames(mm.green) <- colnames(mm.red) <- p.data[,rnb.getOption("identifiers.column")]
	}
	to <- RGChannelSet(Green=mm.green,
			Red=mm.red,
	#		NBeads=mm.beads,
	#		GreenSD=matrix(NA,nrow=nrow(mm.green),ncol=ncol(mm.green)),
	#		RedSD=matrix(NA,nrow=nrow(mm.red),ncol=ncol(mm.red)),
			annotation=anno,
			colData=pheno(from))
	return(to)
})

########################################################################################################################
		
## ---------------------------------------------------------------------------------------------------------------------
## ACCESSORS
## ---------------------------------------------------------------------------------------------------------------------

if(!isGeneric("M")) setGeneric('M',
			function(object, ...) standardGeneric('M'))

#' M-methods
#'
#' Extract raw methylated probe intensity from an object of \code{RnBeadRawSet} class.
#'
#' @param object 		Dataset of interest.
#' @param row.names		Flag indicating whether the resulting matrix will be assigned row names
#'  
#' @return \code{matrix} of the methylated probe intensities
#'
#' @rdname M-methods
#' @docType methods
#' @export
#' @aliases M
#' @aliases M,RnBeadRawSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' M.intensity<-M(rnb.set.example)
#' head(M.intensity)
#' } 
#' 
setMethod("M", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@M, object@meth.regions)
		})


if(!isGeneric("U")) setGeneric('U',
			function(object, ...) standardGeneric('U'))

########################################################################################################################

#' U-methods
#'
#' Extract raw unmethylated probe intensity from an object of \code{RnBeadRawSet} class.
#'
#' @param object 		Dataset of interest.
#' @param row.names		Flag indicating whether the resulting matrix will be assigned row names
#'  
#' @return \code{matrix} of the unmethylated probe intensities
#'
#' @rdname U-methods
#' @docType methods
#' @export
#' @aliases U
#' @aliases U,RnBeadRawSet-method
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' U.intensity<-U(rnb.set.example)
#' head(U.intensity)
#' } 
setMethod("U", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@U, object@meth.regions)
		})

########################################################################################################################

setGeneric('M0',
			function(object, ...) standardGeneric('M0'))

setMethod("M0", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@M0, object@meth.regions)
		})

########################################################################################################################

setGeneric('U0',
			function(object, ...) standardGeneric('U0'))

setMethod("U0", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@U0, object@meth.regions)
		})
########################################################################################################################

setGeneric('bead.counts.M',
			function(object, ...) standardGeneric('bead.counts.M'))


setMethod("bead.counts.M", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@bead.counts.M, object@meth.regions)
		})

########################################################################################################################

setGeneric('bead.counts.U',
			function(object, ...) standardGeneric('bead.counts.U'))


setMethod("bead.counts.U", signature(object="RnBeadRawSet"),
		function(object, row.names=FALSE){
			get.dataset.matrix(object, "sites", row.names, object@bead.counts.U, object@meth.regions)
		})

## ---------------------------------------------------------------------------------------------------------------------
## MODIFIERS
## ---------------------------------------------------------------------------------------------------------------------

setGeneric('M<-',
			function(object, value) standardGeneric('M<-'))

setMethod("M<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@M)
				object@M<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@M<-value
			}
			
		})
########################################################################################################################

setGeneric('U<-',
			function(object, value) standardGeneric('U<-'))

setMethod("U<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@U)
				object@U<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@U<-value
			}
		})

########################################################################################################################

setGeneric('M0<-',
			function(object, value) standardGeneric('M0<-'))

setMethod("M0<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@M0)
				object@M0<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@M0<-value
			}
		})
########################################################################################################################

setGeneric('U0<-',
			function(object, value) standardGeneric('U0<-'))

setMethod("U0<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@U0)
				object@U0<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@U0<-value
			}
		})
########################################################################################################################

setGeneric('bead.counts.M<-',
			function(object, value) standardGeneric('bead.counts.M<-'))

setMethod("bead.counts.M<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@bead.counts.M)
				object@bead.counts.M<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@bead.counts.M<-value
			}
		})
########################################################################################################################

setGeneric('bead.counts.U<-',
			function(object, value) standardGeneric('bead.counts.U<-'))

setMethod("bead.counts.U<-", signature(object="RnBeadRawSet", value="matrixOrffOrNULL"),
		function(object, value){
			if(object@status$disk.dump){
				# delete(object@bead.counts.U)
				object@bead.counts.U<-convert.to.ff.matrix.tmp(value)	
			}else{
				object@bead.counts.U<-value
			}
		})
########################################################################################################################

if (!isGeneric("remove.sites")) {
	setGeneric("remove.sites", function(object, probelist, verbose = TRUE) standardGeneric("remove.sites"))
}

#' @rdname remove.sites-methods
#' @aliases remove.sites,RnBeadRawSet-method
#' @docType methods
#' @export
setMethod("remove.sites", signature(object = "RnBeadRawSet"),
		function(object, probelist, verbose = TRUE) {
			inds <- get.i.vector(probelist, rownames(object@sites))
			if (length(inds) != 0) {
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
						if(!is.null(object@status) && object@status$disk.dump){
							new.matrix<-slot(object,sl)[-inds,,drop=FALSE]
							if(isTRUE(object@status$discard.ff.matrices)){
								delete(slot(object,sl))
							}
							slot(object,sl)<-convert.to.ff.matrix.tmp(new.matrix)
							rm(new.matrix); rnb.cleanMem()
						}else{
							slot(object,sl)<-slot(object,sl)[-inds,,drop=FALSE]
						}
					
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
#' @aliases remove.samples,RnBeadRawSet-method
#' @docType methods
#' @export
setMethod("remove.samples", signature(object = "RnBeadRawSet"),
		function(object, samplelist) {
			inds <- get.i.vector(samplelist, samples(object))
			if (length(inds) != 0) {
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
						if(!is.null(object@status) && object@status$disk.dump){
							new.matrix<-slot(object,sl)[,-inds, drop=FALSE]
							if(isTRUE(object@status$discard.ff.matrices)){
								delete(slot(object,sl))
							}
							slot(object,sl)<-convert.to.ff.matrix.tmp(new.matrix)
							rm(new.matrix); rnb.cleanMem()
						}else{
							slot(object,sl)<-slot(object,sl)[,-inds, drop=FALSE]
						}
					}
				}
			}
			callNextMethod()
		}
)

#######################################################################################################################

#if (!isGeneric("update.meth")) {
	setGeneric("update.meth", function(object) standardGeneric("update.meth"))
#}

##  
## update.meth
## 
## Update the methylation calls, after the change of intensity values
##
## param object 		RnBeadRawSet object
##
## return Updated RnBeadRawSet object
##
setMethod("update.meth", signature(object="RnBeadRawSet"),
		function(object){
			
			if(object@status$disk.dump){
				object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))	
			}else{
				object@meth.sites<-beta.value(object@M, object@U)
			}
			return(object)
		})

#######################################################################################################################
## save, load and destroy 

setMethod("save.matrices", signature(object="RnBeadRawSet", path="character"),
		function(object, path){

			if(!is.null(object@status) && object@status$disk.dump){
				
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
					
						if("ff" %in% class(slot(object,sl))){
							ffmatrix<-slot(object,sl)
							ffsave(ffmatrix, file=file.path(path, paste("rnb", sl, sep=".")),
									rootpath=getOption('fftempdir'))
							rm(ffmatrix)
							
						}
					}
				}
			}
			callNextMethod(object, path)
			
		})

#######################################################################################################################

setMethod("load.matrices", signature(object="RnBeadRawSet", path="character"),
		
		function(object, path, temp.dir=tempdir()){
			slot.names <- RNBRAWSET.SLOTNAMES
			for(sl in slot.names){
				if(!is.null(slot(object, sl))){
					
					if(paste("rnb",sl,"RData", sep=".") %in% list.files(path) &&
							paste("rnb",sl,"ffData", sep=".") %in% list.files(path)){
						load_env<-new.env()
						suppressMessages(ffload(file=file.path(path, paste("rnb", sl, sep=".")), 
										envir=load_env,rootpath=getOption("fftempdir")))
						slot(object, sl)<-get("ffmatrix", envir=load_env)
						rm(load_env)
					}
					
				}

			}
			
			callNextMethod(object=object, path=path, temp.dir=temp.dir)
			
		})

#######################################################################################################################

#' @rdname destroy-methods
#' @aliases destroy,RnBeadRawSet-method
#' @docType methods
#' @export
setMethod("destroy", signature(object="RnBeadRawSet"),
		function(object){
			
			if(object@status$disk.dump){
				for(sl in RNBRAWSET.SLOTNAMES){
					if(!is.null(slot(object,sl))){
						delete(slot(object, sl))	
					}
				}
			}
			callNextMethod()
			
		}
)

## ---------------------------------------------------------------------------------------------------------------------
## HELPER ROUTINES
## ---------------------------------------------------------------------------------------------------------------------

beta.value<-function(M,U,offset=100){
	M/(M+U+offset)
}

#######################################################################################################################

m.value<-function(M,U,offset=100){
	log2((M+offset)/(U+offset))
}

#######################################################################################################################

#' intensities.by.color
#' 
#' Rearranges information from "M" and "U" slots of a RnBeadsRawSet object by color channel.
#'
#' @param raw.set          Methylation dataset as an instance of \code{RnBeadRawSet} object.
#' @param address.rownames  if \code{TRUE} the rows of the returned matrices are named with the with the corresponding Illumina probe addresses 
#' @param add.oob			 if \code{TRUE} the "out-of-band" intensities are included
#' @param add.controls if \code{TRUE} the control probe intensities are included
#' @param add.missing  if \code{TRUE} the rows for the probes missing in \code{raw.set} is imputed with \code{NA} values
#' @param re.separate  if \code{TRUE} the type I and type II intensities, as well as the out-of-band and control probe intensities
#'                     (if set to \code{TRUE}), will be returned as separate elements per channel and not as concatenated rows.  
#' @return A \code{list} with elements \code{Cy3} and \code{Cy5} containing average bead intensities measured for each
#'         each probe in the green and red channels, respectively. Exception, if \code{re.separate} is \code{TRUE} a \code{list}
#'         with elements \code{Cy3.I}, \code{Cy5.I}, and \code{II} will be returned. The elements \code{Cy3.I.oob}, 
#'         \code{Cy5.I.oob} and also \code{Cy3.ctl}, \code{Cy5.ctl} will be returned 
#'         if the respective parameters (\code{add.oob} and \code{add.ctl}) are set to true. 
#' 
#' @author Pavlo Lutsik, Nathan Steenbuck

intensities.by.color<-function(raw.set,
                                   address.rownames = TRUE, 
                                   add.oob = all(!is.null(M0(raw.set)), !is.null(U0(raw.set))), 
                                   add.controls = !is.null(qc(raw.set)),
                                   add.missing = TRUE,
                                   re.separate = FALSE) {
  if (raw.set@target == "probesEPIC") { ## FIXME: No minfi support for EPIC v2!
    rnb.require("IlluminaHumanMethylationEPICmanifest")
    manifest.object <- IlluminaHumanMethylationEPICmanifest

#   } else if (raw.set@target == "probesEPICv2") {
#     rnb.require("IlluminaHumanMethylationEPICv2manifest") ## TODO: Doesn't exist for EPIC v2!
#     manifest.object <- IlluminaHumanMethylationEPICv2manifest

  } else if (raw.set@target == "probes450") {
    rnb.require("IlluminaHumanMethylation450kmanifest")
    manifest.object <- IlluminaHumanMethylation450kmanifest
  } else {
    stop("Unsupported platform")
  }
  
  Mmatrix<-M(raw.set, row.names=TRUE)
  Umatrix<-U(raw.set, row.names=TRUE)
  if(add.oob){
      if(!is.null(M0(raw.set))){
        M0matrix<-M0(raw.set, row.names=TRUE)
        U0matrix<-U0(raw.set, row.names=TRUE)
    }else{
        M0matrix<-matrix(NA_integer_, ncol=ncol(Mmatrix), nrow=nrow(Umatrix))
        U0matrix<-matrix(NA_integer_, ncol=ncol(Umatrix), nrow=nrow(Umatrix))
        rownames(M0matrix)<-rownames(Mmatrix)
        rownames(U0matrix)<-rownames(Umatrix)
    }
  }
  pinfos <- annotation(raw.set, add.names=TRUE)
  
  if(add.missing){
    full.ann<-rnb.annotation2data.frame(rnb.get.annotation(raw.set@target))
    ann.missing<-full.ann[setdiff(rownames(full.ann), rownames(pinfos)),]
    pinfos<-rbind(pinfos, ann.missing[,colnames(full.ann)])
    filler<-matrix(NA_real_, nrow=nrow(ann.missing), ncol=length(samples(raw.set)))
    rownames(filler)<-rownames(ann.missing)
    
    Mmatrix<-rbind(Mmatrix, filler)
    Umatrix<-rbind(Umatrix, filler)
    
    if(add.oob){
      M0matrix<-rbind(M0matrix, filler)
      U0matrix<-rbind(U0matrix, filler)
    }
    
    rm(full.ann, ann.missing, filler)
  }
  
  rnb.set.probe.ids<-pinfos[["ID"]]
  dII.probes <- rnb.set.probe.ids[pinfos[, "Design"] == "II"]
  
  if (address.rownames) {
    tII<-rbind(as.data.frame(manifest.object@data$TypeII[,c("Name", "AddressA")]),
               as.data.frame(manifest.object@data$TypeSnpII[,c("Name", "AddressA")]))
    tII<-tII[match(dII.probes, tII$Name),]
  }
  dII.grn <- Mmatrix[pinfos[, "Design"] == "II", , drop = FALSE]
  dII.red <- Umatrix[pinfos[, "Design"] == "II", , drop = FALSE]
  if (address.rownames) {
    rownames(dII.grn) <- rownames(dII.red) <- tII$AddressA
  }
  
  dI.red.probes <- rnb.set.probe.ids[pinfos[, "Color"] == "Red"]
  dI.green.probes <- rnb.set.probe.ids[pinfos[, "Color"] == "Grn"]
  #dI.red.probes <- dI.red.probes[!grepl("rs", dI.red.probes)]
  #dI.green.probes <- dI.green.probes[!grepl("rs", dI.green.probes)]
  
  if (address.rownames) {
    tI <- c("Name","Color", "AddressA", "AddressB")
    tI <- lapply(c("TypeI", "TypeSnpI"), function(x) { as.data.frame(manifest.object@data[[x]][, tI]) })
    tI <- do.call(rbind, unname(tI))
    
    tI.red<-tI[tI$Color=="Red",]
    tI.red<-tI.red[match(dI.red.probes, tI.red$Name),]
    
    tI.grn<-tI[tI$Color=="Grn",]
    tI.grn<-tI.grn[match(dI.green.probes, tI.grn$Name),]
    rm(tI)
  }
  
  dI.red.meth <- Mmatrix[pinfos[, "Color"] == "Red", , drop = FALSE]
  dI.red.umeth <- Umatrix[pinfos[, "Color"] == "Red", , drop = FALSE]
  dI.grn.meth <- Mmatrix[pinfos[, "Color"] == "Grn", , drop = FALSE]
  dI.grn.umeth <- Umatrix[pinfos[, "Color"] == "Grn", , drop = FALSE]
  
  if (address.rownames) {
    rownames(dI.red.meth) <- tI.red[, "AddressB"]
    rownames(dI.red.umeth) <- tI.red[, "AddressA"]
    rownames(dI.grn.meth) <- tI.grn[, "AddressB"]
    rownames(dI.grn.umeth) <- tI.grn[, "AddressA"]
  }
  if (add.oob) {
    dI.red.meth.oob <- M0matrix[pinfos[, "Color"] == "Red", , drop = FALSE]
    dI.red.umeth.oob <- U0matrix[pinfos[, "Color"] == "Red", , drop = FALSE]
    dI.grn.meth.oob <- M0matrix[pinfos[, "Color"] == "Grn", , drop = FALSE]
    dI.grn.umeth.oob <- U0matrix[pinfos[, "Color"] == "Grn", , drop = FALSE]
    if (address.rownames) {
      rownames(dI.red.meth.oob) <- tI.red[, "AddressB"]
      rownames(dI.red.umeth.oob) <- tI.red[, "AddressA"]
      rownames(dI.grn.meth.oob) <- tI.grn[, "AddressB"]
      rownames(dI.grn.umeth.oob) <- tI.grn[, "AddressA"]
    }
  }
  
  if(re.separate){
    intensities.by.channel <- list("Cy3.I" = list("M" = dI.grn.meth, "U" = dI.grn.umeth),
                                   "Cy5.I" = list("M" = dI.red.meth, "U" = dI.red.umeth), 
                                   "II" = list("M" = dII.grn, "U" = dII.red))
    
    if(add.oob){
      intensities.by.channel[["Cy3.I.oob"]] = list("M" = dI.red.meth.oob, "U" = dI.red.umeth.oob)
      intensities.by.channel[["Cy5.I.oob"]] = list("M" = dI.grn.meth.oob, "U" = dI.grn.umeth.oob)
    } 
  }else{
    intensities.by.channel <- list(
      Cy3=rbind(dII.grn, dI.grn.meth, dI.grn.umeth, 
                if(add.oob) dI.red.meth.oob else NULL, if(add.oob) dI.red.umeth.oob else NULL),
      Cy5=rbind(dII.red, dI.red.meth, dI.red.umeth, 
                if(add.oob) dI.grn.meth.oob else NULL, if(add.oob) dI.grn.umeth.oob else NULL))
    
    rm(dII.grn, dI.grn.meth, dI.grn.umeth, dI.red.meth.oob, dI.red.umeth.oob, 
       dII.red, dI.red.meth, dI.red.umeth, dI.grn.meth.oob, dI.grn.umeth.oob)
    invisible(gc())
    
    if(add.oob & address.rownames) {
      intensities.by.channel$Cy5<-intensities.by.channel$Cy5[rownames(intensities.by.channel$Cy3),,drop=FALSE]
    }
  }
  
  if (add.controls) {
    ncd <- rnb.get.annotation(ifelse(raw.set@target == "probesEPIC", "controlsEPIC", ifelse(raw.set@target == "probesEPICv2", "controlsEPICv2", "controls450")))
    #ncd <- ncd[ncd[["Target"]] == "NEGATIVE", ]
    ncd$Target <- tolower(ncd$Target)
    controls.by.channel <- qc(raw.set)
    controls.by.channel$Cy3 <- controls.by.channel$Cy3[as.character(ncd$ID), , drop=FALSE]
    controls.by.channel$Cy5 <- controls.by.channel$Cy5[as.character(ncd$ID), , drop=FALSE]
    
    if(re.separate){
      intensities.by.channel[["Cy3.ctl"]] <- controls.by.channel$Cy3
      intensities.by.channel[["Cy5.ctl"]] <- controls.by.channel$Cy5
      
    }
    else{
      intensities.by.channel$Cy3 <- rbind(intensities.by.channel$Cy3, controls.by.channel$Cy3)
      intensities.by.channel$Cy5 <- rbind(intensities.by.channel$Cy5, controls.by.channel$Cy5)
    }
  }
  return(intensities.by.channel)
}

########################################################################################################################
