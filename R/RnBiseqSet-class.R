########################################################################################################################
## RnBiseqSet-class.R
## created: 2012-10-29
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## RnBiseqSet class definition.
########################################################################################################################

#' RnBiseqSet Class
#'
#' A class for storing the DNA methylation and quality information from bisulfite sequencing experiments 
#'
#' @details TBA
#'
#' @section Slots:
#' \describe{
#'   \item{\code{status}}{Normalization status.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link[BiocGenerics]{combine}}}{Combines two datasets.}
#' }
#'
#' @name RnBiseqSet-class
#' @rdname RnBiseqSet-class
#' @author Pavlo Lutsik
#' @exportClass RnBiseqSet
setClass("RnBiseqSet",
		representation(status="listOrNULL"),
		contains="RnBSet",
		prototype(pheno=data.frame(), meth.sites=matrix(nrow=0,ncol=0), covg.sites=NULL, target="CpG", status=NULL),
		package = "RnBeads")

## ---------------------------------------------------------------------------------------------------------------------
## CONSTRUCTORS
## ---------------------------------------------------------------------------------------------------------------------

setMethod("initialize", "RnBiseqSet",
	function(.Object,
			pheno=data.frame(),
			sites=matrix(ncol=0,nrow=0),
			meth.sites=matrix(ncol=0,nrow=0),
			covg.sites=matrix(ncol=0,nrow=0),
			regions=list(),
			meth.regions=list(),
			assembly="hg19",
			target="CpG",
			status=list(),
			inferred.covariates=list()
			){
		
		.Object@pheno<-pheno	
		.Object@sites<-sites
		.Object@meth.sites<-meth.sites
		.Object@covg.sites<-covg.sites
		.Object@regions<-regions
		.Object@meth.regions<-meth.regions
		.Object@assembly<-assembly
		.Object@target<-target
		.Object@status<-status
		.Object@inferred.covariates <- inferred.covariates
		
		.Object
		
	}
)
		

#' Wrapper function RnBiseqSet
#' 
#' @param pheno       	phenotypic data.
#' @param sites       	CpG site definition, as a \code{data.frame}
#' 						with 3 variables: chromosome (of type
#' 						\code{character}), position (\code{integer}) and
#' 						strand (\code{character}, one of \code{"+"}, \code{"-"} or \code{"*"}
#' @param meth  		summarized methylation calls as a \code{matrix} or \code{ff_matrix}
#' @param covg 			read coverage information as a \code{matrix} or \code{ff_matrix}
#' @param assembly		the genome assembly
#' @param target		target DNA methylation features (CpG sites)
#' @param summarize.regions ...
#' @param region.types	region annotations for which the methylation data should be summarized 
#' @param useff			flag specifying whether the ff functionality should be used
#' @param usebigff		flag specifying whether the extended ff functionality should be used (large matrix support for ff)
#' @param verbose		flag specifying whether the diagnostic messages should be written to the 
#' 						console or to the RnBeads logger, if the latter is initialized
#' 
#' @return an object of class RnBiseqSet
#' 
#' @name RnBiseqSet
#' @rdname RnBiseqSet-class
#' @aliases initialize,RnBiseqSet-method
#' @export
RnBiseqSet<-function(
		pheno,
		sites,
		meth,
		covg = NULL,
		assembly="hg19",
		target="CpG",
		summarize.regions=TRUE,
		region.types=rnb.region.types.for.analysis(assembly),
		useff=rnb.getOption("disk.dump.big.matrices"),
		usebigff=rnb.getOption("disk.dump.bigff"),
		verbose=FALSE){
	
	if(missing(pheno)){
		stop("argument pheno should be supplied")
	}
	
	if(missing(sites)){
		#warning("Valid RnBeadSet object was not created: pheno and/or betas missing")
		stop("argument meth should be supplied")
	}
	
	if(missing(meth)){
		#warning("Valid RnBeadSet object was not created: pheno and/or betas missing")
		stop("argument meth should be supplied")
	}
	
	if(!is.data.frame(pheno)){
		stop("invalid value for pheno: should be a data frame")
	}
	
	if(!(any(c('data.frame', 'matrix') %in% class(sites)))){
		stop("invalid value for sites: should be a data frame or a matrix of type")
	}
	
	if(!any(c('matrix', 'ff_matrix', 'BigFfMat') %in% class(meth))){
		stop("invalid value for meth: should be a matrix or an ff_matrix")
	}
	
	if(!is.null(covg) && !any(c('matrix', 'ff_matrix', 'BigFfMat') %in% class(covg))){
		stop("invalid value for covg: should be a matrix or an ff_matrix")
	}
	
	if(!is.character(assembly) || length(assembly)!=1L){
		stop("invalid value for assembly: should be a character of length one")
	}
	
	if(!is.character(target) || length(target)!=1L){
		stop("invalid value for target: should be a character of length one")
	}
	
	if(!is.logical(summarize.regions) || length(summarize.regions)!=1L){
		stop("invalid value for summarize.regions: should be a logical of length one")
	}
	
	if(!is.logical(useff) || length(useff)!=1L){
		stop("invalid value for useff: should be a logical of length one")
	}
	
	if (nrow(sites)!= nrow(meth)){
		msg<-c("Inconsistent values for sites and meth")
		rnb.error(msg)
	}
	if(usebigff){
		bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")
	}
	
	if(!is.null(rnb.getOption("identifiers.column"))){
		sample.names<-pheno[[rnb.getOption("identifiers.column")]]
	}else if(!is.null(colnames(meth))){
		sample.names<-colnames(meth)
	}else{
		sample.names<-as.character(1:ncol(meth))
	}
	#
	# prepare sites definition
	#
	chroms<-names(rnb.get.chromosomes(assembly=assembly))
	
	if(is.data.frame(sites)){
		sites.df<-sites
		sites<-matrix(nrow=nrow(sites.df), ncol=3, dimnames=list(NULL, c("chr", "coord", "strand")))
		
		sites[,1L]<-as.integer(match(sites.df[,1L], chroms))
		sites[,2L]<-as.integer(sites.df[,2L])
		sites[,3L]<-as.integer(factor(sites.df[,3L], levels=c("+","-","*")))
		
		rm(sites.df)
	}
	
	num.all.sites<-nrow(sites)
	
	#remove sites with illegal chromosomes
	legal.sites <- rowSums(is.na(sites))==0
	num.illegal <- sum(!legal.sites)
	
	sites <- sites[legal.sites,,drop=FALSE]
	#meth <- meth[legal.sites,,drop=FALSE]
	#covg <- covg[legal.sites,,drop=FALSE]
	
	if(verbose && num.illegal > 0) {
		rnb.info(c("Removed", num.illegal, "sites with unknown chromosomes"))
	}
	
	if(nrow(sites)<1L){
		rnb.warning("All sites have been removed, returning NULL")
		return(invisible(NULL))
	}
		
	#
	# match to annotation
	#
	guess.strand <- all(sites[,3L]==3)
	if(verbose && guess.strand) {
		rnb.info("Inferring strand information from annotation enabled")
	}
	# Chromosome ordering fix. 
	#
	#chrdata<-sapply(sort(unique(sites[,1L])), function(chr){
	chrdata<-lapply(sort(unique(sites[,1L])), function(chr){
				#dump<-sapply(unique(sites[,1L]), function(chr){
				chr.subset<-which(sites[,1L]==chr)
				# version with proper ordering
				#
				# chr.subset<-sites[,1L]==chr
				if(guess.strand){
					site.ranges <- GRanges(rep(chroms[chr], length(chr.subset)), 
							IRanges(start=sites[chr.subset,2], 
									width=rep(1,length(chr.subset))), 
							rep("*",length(chr.subset)))
					c.coords.ref <- GenomicRanges::flank(rnb.get.annotation(assembly=assembly)[[chr]],-1)
					
					match <- GenomicRanges::findOverlaps(c.coords.ref, site.ranges, type="start")
				} else {
					site.ranges<-GRanges(rep(chroms[chr], length(chr.subset)), 
							IRanges(start=sites[chr.subset,2L], 
									width=rep(1,length(chr.subset))), 
							c("+","-","*")[sites[chr.subset,3L]])
					
					#				match<-GenomicRanges::findOverlaps(site.ranges, rnb.get.annotation(assembly=assembly)[[chr]], type="start")
					#				sites[chr.subset[queryHits(match)],2]<<-subjectHits(match)
					#				sites[chr.subset[-queryHits(match)],2]<<-NA
					#				return(NULL)
					
					# version with proper ordering
					#			
					match<-GenomicRanges::findOverlaps(rnb.get.annotation(assembly=assembly)[[chr]], site.ranges, type="start")
				}
				bed.h <- subjectHits(match)
				ref.h <- queryHits(match)
				multi.hit.ref <- duplicated(ref.h)
				if(any(multi.hit.ref)){
					rnb.warning(paste("Multiple entries detected in the bed file that match to the same reference C. --> Picking the first one"))
					bed.h <- bed.h[!multi.hit.ref]
					ref.h <- ref.h[!multi.hit.ref]
				}
				multi.hit.bed <- duplicated(bed.h)
				if(any(multi.hit.bed)){
					rnb.warning(paste("Some entries in the bed file match multiple reference Cs. --> Picking the first one"))
					bed.h <- bed.h[!multi.hit.bed]
					ref.h <- ref.h[!multi.hit.bed]
				}
				
				n.match <- length(match)
				
				#chrmeths<-meth[chr.subset,,drop=FALSE][bed.h,,drop=FALSE]
				#chrcovgs<-covg[chr.subset,,drop=FALSE][bed.h,,drop=FALSE]
				#indices<-chr.subset[bed.h]
				chrsites<-cbind(rep(1L,n.match),rep(chr, n.match), ref.h, chr.subset[bed.h])
				#return(list(chrsites,chrmeths,chrcovgs))
				chrsites
				
			})
	
	rnb.cleanMem()
	#version with proper ordering
	#sites<-do.call("rbind", chrdata[1,])
	sites<-do.call("rbind", chrdata)
	#meth<-do.call("rbind", chrdata[2,])
	#covg<-do.call("rbind", chrdata[3,])
	valid<-rowSums(is.na(sites))==0
	
	
	if(nrow(sites)<1L){
		rnb.warning("All sites have been removed, returning NULL")
		return(invisible(NULL))
	}else{
	
		if(verbose){
			rnb.status(c("Matched", sum(valid),"of", num.all.sites,"methylation sites to the annotation"))
		}
	}
	
	if(!is.null(covg)){
		if(verbose) rnb.status(c("Checking site coverage")) #TODO: debug info: remove me
		zero.covg.sites <- rowSums(covg[sites[,4L],,drop=FALSE], na.rm=TRUE)==0 #TODO: could this be too memory wasteful?
		if(verbose && any(zero.covg.sites)) {
			rnb.status(c("Removed",length(intersect(which(zero.covg.sites), which(valid))),"of",
							sum(valid),"methylation sites because they were not covered in any sample"))
		}
		valid<-which(valid | zero.covg.sites)
		rnb.cleanMem()
		if(!useff){
			covg[is.na(covg)]<-0L
		}else{
			# if(verbose) rnb.status(c("Setting zero coverages")) #TODO: debug info: remove me
			for(ci in 1:ncol(covg)){
				covg[is.na(covg[,ci]),ci]<-0L
			}
		}
	}
	rnb.cleanMem()
	
	sites <- sites[valid,]
	if(!useff){
		meth <- meth[which(legal.sites)[sites[,4L]],,drop=FALSE]
		if(!is.null(covg)){
			covg <- covg[which(legal.sites)[sites[,4L]],,drop=FALSE]
		}
	}else{
		meth.old <- meth
		if(verbose) rnb.status(c("Creating methylation matrix")) #TODO: debug info: remove me
		if (usebigff) {
			meth <- BigFfMat(row.n=nrow(sites), col.n=ncol(meth.old), row.names=NULL, col.names=sample.names, finalizer=bff.finalizer)
		} else {
			meth <- ff(NA, dim=c(nrow(sites), ncol(meth.old)), dimnames=list(NULL, sample.names), vmode="double")
		}
		# if(verbose) rnb.status(c("Filling methylation matrix")) #TODO: debug info: remove me
		# meth[,] <- meth.old[which(legal.sites)[sites[,4L]],]
		site.inds <- which(legal.sites)[sites[,4L]]
		for (j in 1:ncol(meth.old)) meth[,j] <- meth.old[site.inds,j]
		rm(meth.old)
		rnb.cleanMem()
		if(!is.null(covg)){
			covg.old <- covg
			if(verbose) rnb.status(c("Creating coverage matrix")) #TODO: debug info: remove me
			if (usebigff) {
				covg <- BigFfMat(row.n=nrow(sites), col.n=ncol(covg.old), row.names=NULL, col.names=sample.names, na.prototype=as.integer(NA), finalizer=bff.finalizer)
			} else {
				covg <- ff(NA_integer_, dim=c(nrow(sites), ncol(covg.old)), dimnames=list(NULL, sample.names))
			}
			# if(verbose) rnb.status(c("Filling coverage matrix")) #TODO: debug info: remove me
			# covg[,] <- covg.old[which(legal.sites)[sites[,4L]],]
			for (j in 1:ncol(covg.old)) covg[,j] <- covg.old[site.inds,j]
			rm(covg.old)
			rnb.cleanMem()
		}
	}
	sites <- sites[,-4L]
	colnames(sites)<-c("context", "chr", "index")
	
	meth.sample.n.valid<-integer()
	for(ci in 1:ncol(meth)){
		meth.sample.n.valid[ci] <- sum(!is.na(meth[,ci]))
	}
	
	if(any(meth.sample.n.valid<1)){
		rnb.error(c("The following samples have no valid methylation values:",colnames(meth)[meth.sample.n.valid<1]))
		#stop("invalid sample methylation values")
	}
	
#	sites<-cbind(rep(1, nrow(sites)), sites)
#	sites<-sites[,2:4]
#	colnames(sites)<-c("context", "chr", "index")
	
	#be memory efficient by cleaning up
	# rownames take a lot of space
	
#	sites <- sites[valid,,drop=FALSE]
#	meth  <- meth[valid,,drop=FALSE]
#	covg  <- covg[valid,,drop=FALSE]
	rnb.cleanMem()
#	if(!is.null(covg)){
#		zero.covg.sites <- rowSums(covg, na.rm=TRUE)==0
#		if(any(zero.covg.sites)){
#			valid <- !zero.covg.sites
#			sites <- sites[valid,,drop=FALSE]
#			meth  <- meth[valid,,drop=FALSE]
#			covg  <- covg[valid,,drop=FALSE]
#			rnb.cleanMem()
#
#			if(verbose) {
#				rnb.status(c("Removed",sum(zero.covg.sites),"of",
#								length(valid),"methylation sites because they were not covered in any sample"))
#			}
#		}
#	}
	
	status <- list(disk.dump=FALSE)
	if(useff){
		#subsampling for large datasets (ff currently supports only objects of size .Machine$integer.max)
		if (prod(dim(meth))>.Machine$integer.max & !usebigff){
			sites.allowed <- as.integer(.Machine$integer.max/ncol(meth))
			sample.site.inds <- sort(sample.int(nrow(meth),sites.allowed))
			msg<-c("Full dataset is too large to be supported by ff. --> downsampling to",sites.allowed,"( of",nrow(meth),") sites")
			rnb.warning(msg)
			meth <- meth[sample.site.inds,]
			if(!is.null(covg)){
				covg <- covg[sample.site.inds,]
			}
			sites <- sites[sample.site.inds,]
		}
		
		if("matrix" %in% class(meth)){
			if (usebigff){
				meth <- BigFfMat(meth, finalizer=bff.finalizer)
			} else {
				meth <- convert.to.ff.matrix.tmp(meth)
			}
		}
		
		if(!is.null(covg)){
			if("matrix" %in% class(covg)){
				if (usebigff){
					covg <- BigFfMat(covg, finalizer=bff.finalizer)
				} else {
					covg <- convert.to.ff.matrix.tmp(covg)
				}
			}
		}
		status <- list(disk.dump=TRUE, disk.dump.bigff=usebigff)
	}

	if(verbose) rnb.status(c("Creating object")) #TODO: debug info: remove me
	object<-new("RnBiseqSet",
			pheno=pheno,
			sites=sites,
			meth.sites=meth,
			covg.sites=covg,
			status=status,
			assembly=assembly,
			target=target)
	
	if(!rnb.getOption("strand.specific")){
		msg<-c("Summarizing","strand","methylation")
		if(verbose){
			rnb.status(msg)
		}
		if(is.null(covg)){
			aggr<-"mean"
		}else{
			aggr<-"coverage.weighted"
		}
		object <- summarize.regions(object, "strands", aggregation=aggr)
		rnb.cleanMem()
	}
	
	for (region.type in region.types) {
		if (region.type %in% rnb.region.types(assembly)) {
			msg<-c("Summarizing",region.type,"methylation")
			if(verbose){
				rnb.status(msg)
			}
			object <- summarize.regions(object, region.type)
			rnb.cleanMem()
		}
	}
	object
}


## ---------------------------------------------------------------------------------------------------------------------
## VALIDITY
## ---------------------------------------------------------------------------------------------------------------------
#
# check.rnb.biseq.set
#
# A routine that checks the validity of a freshly loaded RnBiseqSet object
# 
# #param object RnBiseqSet object to test
# #param verbose A flag specifying whether the logger is initialized
#
# #return TRUE if all checks are successful
# 
# Pavlo Lutsik
#
check.rnb.biseq.set<-function(object, verbose=TRUE){
	
	
	if(!inherits(object, "RnBiseqSet")){
		if(verbose){ 
			rnb.info("The supplied object is not of class RnBiseqSet. Breaking the check...")
		}
		return(FALSE)
	}else{
		if(verbose){
			rnb.info("Checking the supplied RnBiseqSet object")
		}
		valid<-TRUE
	}
	
	#check the methylation data
	
	mm <- meth(object)
	
	if(nrow(mm)<1){
		if(verbose){
			rnb.info("The object contains information for 0 methylation sites")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info(sprintf("The object contains information for %d methylation sites",nrow(mm)))
		}
	}
	
	if(ncol(mm)<1){
		if(verbose){
			rnb.info("The object contains information for 0 samples")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info(sprintf("The object contains information for %d samples",ncol(mm)))
		}
	}
	
	if(sum(!is.na(mm))<1){
		if(verbose){
			rnb.info("All methylation values are missing")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info(sprintf("The object contains %d missing methylation values",sum(is.na(mm))))
		}
	}
	
	if(sum(mm[!is.na(mm)]<0 & mm[!is.na(mm)]>1)>0){
		if(verbose){
			rnb.info("The object contains incorrect methylation values (less than 0 or greater than 1)")
		}
		valid<-FALSE
	}else{
		if(verbose){
			rnb.info("Methylation values are within the expected range")
		}
	}
	
	#check the coverage data
	
	if(is.null(object@covg.sites)){
		if(verbose){
			rnb.info("No coverage information found")
		}
		valid<-FALSE
		return(valid)
	}else{
		if(verbose){
			rnb.info("The object contains coverage information")
		}
		cvg<-covg(object)
	}
	
	if(sum(cvg[!is.na(cvg)]<0)>0){
		if(verbose){
			rnb.info("The object contains incorrect coverage values (negative values)")
		}
		valid<-FALSE
	}else if(sum(cvg[!is.na(cvg)]>1000)>0.5*length(as.numeric(cvg))){
		if(verbose){
			rnb.info("The object contains incorrect coverage values (half of the sites have coverage > 1000 )")
		}
	}else{
		if(verbose){
			rnb.info("Coverage values are within the expected range")
		}
	}
	
	return(valid)
	
}

########################################################################################################################

setMethod("show", "RnBiseqSet",
	function(object) {
		cat("Object of class RnBiseqSet\n")
		cat(sprintf("%8d samples\n", nrow(object@pheno)))
		cat(sprintf("%8d methylation sites\n", nrow(object@sites)))
		if (!is.null(object@regions)){
			cat("Region types:\n")
			for (rn in names(object@regions)) {
				cat(sprintf("\t%8d regions of type %s\n", nrow(object@regions[[rn]]), rn))
			}
		}
		cat(sprintf("Coverage information is %s\n", ifelse(is.null(object@covg.sites), "absent", "present")))
	}
)

########################################################################################################################
