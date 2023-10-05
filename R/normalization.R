########################################################################################################################
## normalization.R
## created: 2012-05-31
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the normalization step.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.execute.normalization
#'
#' Performs normalization of the provided HumanMethylation450 data set.
#'
#' @param object Methylation dataset as an object of type \code{\linkS4class{MethyLumiSet}} or
#'               \code{\linkS4class{RnBSet}}.
#' @param method Normalization method, must be one of \code{"none"}, \code{"illumina"}, \code{"swan"},
#'               \code{"minfi.funnorm"}, \code{"bmiq"}, or \code{wm.*} where \code{*} stands for one of the methods
#'               implemented in \pkg{wateRmelon} package.
#'               Note that the execution of methods SWAN and minfi.funnorm requires packages \pkg{minfi} and
#'               \pkg{IlluminaHumanMethylation450kmanifest}. The BMIQ method requires the package \pkg{RPMM}. The
#'               \code{wm.*} methods naturally require \pkg{wateRmelon}.
#' @param bgcorr.method Character singleton specifying which background subtraction should be used. Only methods impemented
#'               in the \pkg{methylumi} package are supported at the moment, namely \code{methylumi.noob}, \code{methylumi.goob}
#'               and \code{methylumi.doob}. See Triche et al. for detailed description of the methods.
#' @param verbose flag specifying whether diagnostic output should be written to the console or to the RnBeads logger
#' 				 in case the latter is initialized
#'
#' @return Normalized dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#'
#' @references 1. Triche, Timothy J., Jr., Weisenberger, Daniel J., Van Den Berg, David, Laird, Peter W. and Siegmund, Kimberly D. (2013)
#'             Low-level processing of Illumina Infinium DNA Methylation BeadArrays.
#'             Nucleic Acids Research 41(7):e90-e90.
#'
#' @author Pavlo Lutsik
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' rnb.set.norm<-rnb.execute.normalization(rnb.set.example, method="illumina", bgcorr.method="none")
#' }
#' @export
rnb.execute.normalization<-function(
		object,
		method=rnb.getOption("normalization.method"),
		bgcorr.method=rnb.getOption("normalization.background.method"),
		verbose=TRUE){

	if(!(inherits(object,"MethyLumiSet") || inherits(object,"RnBSet"))) {
		stop("invalid value for object; expected MethyLumiSet, RnBeadSet or RnBiseqSet")
	}
	accepted <- .rnb.options[["accepted"]][["normalization.method"]]
	accepted <- c(accepted, paste0(setdiff(accepted, c("none", "bmiq")), "+bmiq"))
	if (!(is.character(method) && length(method) == 1 && isTRUE(method %in% accepted))) {
		msg <- paste0('"', accepted, '"', collapse = ", ")
		stop(paste("invalid value for method; expected one of", msg))
	}
	accepted <- .rnb.options[["accepted"]][["normalization.background.method"]]
	if (!(is.character(bgcorr.method) && length(bgcorr.method) == 1 && isTRUE(bgcorr.method %in% accepted))) {
		stop("invalid value for bgcorr.method")
	}

	## Split the double normalization
	if (grepl("\\+bmiq$", method)) {
		method <- substr(method, 1, nchar(method) - 5)
		secondary.bmiq <- TRUE
	} else {
		secondary.bmiq <- FALSE
	}
	method.to.set <- method
	bgcorr.method.to.set <- bgcorr.method
	disable.method <- function(is.method, txt = '') {
		txt <- paste0('Incompatible dataset and ', ifelse(is.method, 'normalization', 'background correction'),
					  ' method', ifelse(txt == '', '', ': '), txt, '. Changed the method to "none"')
		rnb.warning(txt)
		if (is.method) {
		  rnb.options(normalization.method = "none")
			method <<- "none"
			method.to.set <<- object@status$normalized
			secondary.bmiq <<- FALSE
		} else {
			rnb.options(normalization.background.method = "none")
			bgcorr.method <<- "none"
			bgcorr.method.to.set <<- object@status$background
		}
	}

	## Ignore normalization for bisulfite sequencing datasets
	if (inherits(object, "RnBiseqSet")) {
		if (method != "none") {
			disable.method(TRUE)
		}
		if (bgcorr.method != "none") {
			disable.method(FALSE)
		}
		return(object)
	}

	## Validate there are enough samples present
	if (inherits(object, "MethyLumiSet")) {
		nsamples <- ncol(exprs(object))
	} else {
		nsamples <- ncol(object@meth.sites)
	}
	accepted <- c("none", "bmiq")
	if (nsamples < 2 && !(method %in% accepted)) {
		disable.method(TRUE, 'too few samples')
	}
	rm(nsamples)

	if (inherits(object, "RnBeadSet")) {
		if (object@status$normalized != "none") {
			## Allow only BMIQ as a secondary normalization
			if (method == "bmiq" && (!grepl("bmiq$", object@status$normalized))) {
				method.to.set <- paste0(object@status$normalized, "+bmiq")
			} else if (method != "none") {
				disable.method(TRUE, 'dataset already normalized')
			}
		} else if (!(method %in% c("none", "bmiq") || inherits(object, "RnBeadRawSet"))) {
			disable.method(TRUE, 'missing intensity data')
		} else if (method == "illumina" && is.null(qc(object))) {
			disable.method(TRUE, 'missing data on quality control probes')
		}
	}

	## Perform background subtraction
	if (bgcorr.method != "none") {
		if (inherits(object, "MethyLumiSet")) {

			if (annotation(object) == "IlluminaHumanMethylation27") {
				disable.method(FALSE, 'not supported for Infinium 27k')
			} else if (grepl("enmix", bgcorr.method)) {
				disable.method(FALSE, 'not supported for MethyLumiSet')
			} else if (bgcorr.method == "methylumi.noob") {
				if (!all(c("methylated.OOB", "unmethylated.OOB") %in% ls(object@assayData))) {
					disable.method(FALSE, 'missing out-of-band intensities')
				}
			}

			if (grepl("methylumi", bgcorr.method)) {
				bgcorr.methylumi<-gsub("methylumi\\.", "", bgcorr.method)
				pheno.columns<-colnames(phenoData(object)@data)
				suppressMessages({
					sinkfile<-ifelse("Windows" %in% Sys.info(),"NUL", "/dev/null")
					sink(sinkfile)
					object<-methylumi.bgcorr(object, method=bgcorr.methylumi)
					sink()
				})
				#removing the introduced columns
				phenoData(object)<-phenoData(object)[,pheno.columns]
			}

		} else if (inherits(object, "RnBeadRawSet")) {

			if (object@target == "probes27") {
				disable.method(FALSE, 'not supported for Infinium 27k')
			} else if (bgcorr.method == "methylumi.noob" && object@target == "probesEPICv2") {
				disable.method(FALSE, 'methylumi.noob is not supported for MethylationEPICv2')
            } else if (bgcorr.method == "methylumi.noob" && object@target == "probesEPIC") {
				disable.method(FALSE, 'methylumi.noob is not supported for MethylationEPIC')
            } else if (bgcorr.method == "methylumi.noob" && object@target == "probesMMBC") {
                disable.method(FALSE, 'methylumi.noob is not supported for MethylationMMBC')
			} else if (grepl("oob$", bgcorr.method) && (is.null(M0(object)) || is.null(U0(object)))) {
				disable.method(FALSE, 'missing out-of-band intensities')
			} else if (object@status$normalized == "swan") {
				disable.method(FALSE, 'dataset already normalized using SWAN')
			}

			if (bgcorr.method=="subtraction"){
                object<-rnb.bgcorr.subtr(object)
            }else if (any(grepl("sesame.noob", bgcorr.method))){
                if(bgcorr.method=="sesame.noobsb"){
                    object<-noob.sesame(object, background="bleed-through")
                }else{
                    object<-noob.sesame(object, background="internal")
                }
            }else if (any(grepl("methylumi", bgcorr.method))) {
				bgcorr.methylumi<-gsub("methylumi\\.", "", bgcorr.method)
				pheno.columns<-colnames(pheno(object))
				inferred.covariates <- object@inferred.covariates
				old.obj<-object
				object<-as(object, "MethyLumiSet")
				if(isTRUE(old.obj@status$discard.ff.matrices)){
					rnb.call.destructor(old.obj)
					rm(old.obj)
				}
				rnb.cleanMem()
				suppressMessages({
					sinkfile <- ifelse("Windows" %in% Sys.info(), "NUL", "/dev/null")
					sink(sinkfile)
					object<-methylumi.bgcorr(object, method=bgcorr.methylumi)
					sink()
				})
				#removing the introduced columns
				phenoData(object)<-phenoData(object)[,pheno.columns]
				object<-as(object, "RnBeadRawSet")
				if(object@status$disk.dump && rnb.getOption("enforce.destroy.disk.dumps")){
					object@status$discard.ff.matrices<-TRUE
				}
			} else if (any(grepl("enmix", bgcorr.method))) {
				bgcorr.enmix<-gsub("enmix\\.", "", bgcorr.method)
				object<-rnb.enmix.oob(object)
			}

			object@status$background <- bgcorr.method
		}
		if (verbose) {
			rnb.status(c("Performed background subtraction with method", bgcorr.method))
		}
		rnb.cleanMem()
	}

	## Validate the normalization method is supported for the given data type
	accepted <- c("none", "illumina", "methylumi.illumina")
	if (inherits(object, "MethyLumiSet") && annotation(object) == "IlluminaHumanMethylation27" &&
		!(method %in% accepted)) {
		disable.method(TRUE, 'not supported for Infinium 27k MethyLumiSet')
	}
	if (inherits(object, "RnBeadSet") && object@target == "probes27" && !(method %in% accepted)) {
		disable.method(TRUE, 'not supported for Infinium 27k')
	}
	accepted <- c("none", "scaling", "scaling.internal", "scaling.reference", "bmiq", "swan", "minfi.funnorm", "wm.dasen")
	## TODO: No support for EPICv2?
	if (inherits(object, "RnBeadSet") && object@target == "probesEPIC" && !(method %in% accepted)) {
		disable.method(TRUE, 'not supported for HumanMethylationEPIC')
	}
    accepted<-setdiff(.rnb.options[["accepted"]][["normalization.method"]], c("illumina", "minfi.funnorm", "swan", "wm.tost","wm.fuks", "wm.swan"))
    if (inherits(object, "RnBeadSet") && object@target == "probesMMBC" && !(method %in% accepted)) {
        disable.method(TRUE, 'not supported for Mouse Methylation Bead Chip')
    }
    
	## Perform normalization
	if (any(grepl("scaling", method))) {
        
        if(method=="scaling") {
            sc.method="internal"
        }else{
            sc.method<-gsub("scaling.", "", method)
        }
        
        object<-rnb.norm.scaling(object, method=sc.method)
        
        object@status$normalized<-"scaling"
        object@status$background<-bgcorr.method
        rnb.cleanMem()
        
    } else if (method=="illumina" || method == "methylumi.illumina") {
		if (inherits(object, "RnBeadRawSet")) {
			inferred.covariates <- object@inferred.covariates
			object <- as(object, "MethyLumiSet")
		}
		object <- suppressMessages(as(normalizeMethyLumiSet(object), "RnBeadRawSet"))
		object@status$normalized<-"illumina"
		object@status$background<-bgcorr.method
		rnb.cleanMem()

	}else if (method=="swan"){

		rnb.require("minfi")
		rnb.require("IlluminaHumanMethylation450kmanifest") 
		if(inherits(object,"MethyLumiSet") && (is.null(methylated(object))||is.null(unmethylated(object)))) {
			rnb.error("Invalid value for object; missing intensity information")
		}
		## TODO: No support for EPIC v2?
		rga <- c("IlluminaHumanMethylationEPIC", "ilm10b2.hg19", "IlluminaHumanMethylation450k", "ilmn12.hg19")
		rga <- matrix(rga, 2, 2, TRUE, list(c("EPIC", "450"), c("array", "annotation")))
		if(inherits(object,"MethyLumiSet")){
			intensities.by.channel<-methylumi.intensities.by.color(object)
			## FIXME: Update this MethyLumiSet can contain EPIC data as well
			rga <- rga["450", ]
		}else if(inherits(object,"RnBeadRawSet")){
			intensities.by.channel<-intensities.by.color(object)
			rga <- rga[gsub("^probes", "", object@target), ]
		}
		if (grepl("EPIC", rga[1])) { ## FIXME: No EPICv2 manifest from minfi! This condition passes for EPICv2
			rnb.require("IlluminaHumanMethylationEPICmanifest")
		} else {
			rnb.require("IlluminaHumanMethylation450kmanifest")
		}

		rg.set<-RGChannelSet(intensities.by.channel$Cy3, intensities.by.channel$Cy5)
		annotation(rg.set) <- rga
		suppressMessages({
				sinkfile<-ifelse("Windows" %in% Sys.info(), "NUL", "/dev/null")
				sink(sinkfile); methyl.set<-preprocessSWAN(rg.set); sink()
		})

		meth.minfi<-getMeth(methyl.set)
		umeth.minfi<-getUnmeth(methyl.set)

		if(inherits(object, "MethyLumiSet")){
			methylated(object)<-meth.minfi[match(rownames(meth.minfi), featureNames(object)),]#+methylated(object)[setdiff(featureNames(object), rownames(meth.minfi)),]
			unmethylated(object)<-umeth.minfi[match(rownames(umeth.minfi), featureNames(object)),]#+unmethylated(object)[setdiff(featureNames(object), rownames(umeth.minfi)),]
			rm(rg.set,methyl.set, meth.minfi, umeth.minfi)
			betas(object)<-rbind(methylated(object)/(methylated(object)+unmethylated(object)),
				betas(object)[setdiff(featureNames(object), rownames(methylated(object))),])
			object<-as(object, "RnBeadSet")
		}else if(inherits(object, "RnBeadRawSet")){
			probe.ids<-rownames(annotation(object))
			## rs probes are "lost" during SWAN

			meth.minfi<-meth.minfi[rownames(meth.minfi) %in% probe.ids,]
			umeth.minfi<-umeth.minfi[rownames(umeth.minfi) %in% probe.ids,]

			object@M[,][match(rownames(meth.minfi),probe.ids),]<-meth.minfi
			object@U[,][match(rownames(umeth.minfi),probe.ids),]<-umeth.minfi

			#update.meth(object)
			if(object@status$disk.dump){
				object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))
			}else{
				object@meth.sites<-beta.value(object@M[,,drop=FALSE], object@U[,,drop=FALSE])
			}
			rm(rg.set,methyl.set, meth.minfi, umeth.minfi)
		}

	} else if (method == "bmiq") {

		object <- rnb.execute.normalization.bmiq(object)

	}else if(grepl("wm\\.",method)[1]){

		rnb.require("wateRmelon")
		wm.method<-gsub("wm\\.","", method)

		if(inherits(object, "MethyLumiSet")){

			object<-do.call(wm.method, list(object))

			object<-as(object, "RnBeadSet")

		}else if(inherits(object, "RnBeadRawSet")){

			if(wm.method %in% c("betaqn", "fuks")){

				ann<-annotation(object, add.names=TRUE)

				if(wm.method=="betaqn"){
					betas.norm<-do.call(wm.method, list(meth(object)))
				}else{
					colnames(ann)[5]<-"DESIGN"
					betas.norm<-do.call(wm.method, list(meth(object), anno=ann))
				}
				object<-as(object, "RnBeadSet")
				if(object@status$disk.dump){
					object@meth.sites<-convert.to.ff.matrix.tmp(betas.norm)
				}else{
					object@meth.sites<-betas.norm
				}
				rm(betas.norm)

			}else{

				Mmatrix<-M(object, row.names=TRUE)
				Umatrix<-U(object, row.names=TRUE)

				ann.full<-rnb.annotation2data.frame(rnb.get.annotation(object@target, object@assembly))
				ann<-annotation(object, add.names=TRUE)
				probe.names<-ann[["ID"]]

				if(nrow(Mmatrix)<nrow(ann.full)){
					filler<-matrix(NA_real_, nrow=nrow(ann.full)-nrow(ann), ncol=length(samples(object)))
					rownames(filler)<-rownames(ann.full[!rownames(ann.full) %in% rownames(ann),])
					ann<-rbind(ann, ann.full[!rownames(ann.full) %in% rownames(ann),colnames(ann.full)])
					Mmatrix<-rbind(Mmatrix, filler)
					Umatrix<-rbind(Umatrix, filler)
					rm(filler)
				}
				rm(ann.full)
				rownames(Mmatrix)<-ann[["ID"]]
				rownames(Umatrix)<-ann[["ID"]]

				if(wm.method %in% c("swan")){

					if(wm.method %in% c("swan")){
						betas.norm<-suppressMessages({do.call(wm.method, list(Mmatrix, Umatrix,qc=qc(object)))})
						rnb.warning("Methylation calls at SNP-tagging probes were not normalized")
					}

#					rownames(betas.norm)<-probe.names
#					qcl<-qc(object)
#					object<-new("RnBeadSet",
#							pheno(object),
#							betas.norm[1:nrow(meth(object)),],
#							dpval(object, row.names=TRUE),
#							covg(object, row.names=TRUE))
#					qc(object)<-qcl
					object<-as(object, "RnBeadSet")
					betas.new<-meth(object)
					betas.new[match(rownames(betas.norm), probe.names),]<-betas.norm
					if(object@status$disk.dump){
						object@meth.sites<-convert.to.ff.matrix.tmp(betas.new)
					}else{
						object@meth.sites<-betas.new
					}
					rm(betas.norm, betas.new)

				}else if(wm.method %in% c("tost")){
					cgi.relation<-rep("Sea",nrow(ann))
					cgi.relation[unlist(regionMapping(object, "cpgislands"))]<-"Island"
					ann.tost<-data.frame(
							TargetID=ann[["ID"]],
							#TargetID=as.character(1:nrow(ann)),
							COLOR_CHANNEL=ann[["Color"]],
							INFINIUM_DESIGN_TYPE=ann[["Design"]],
							RELATION_TO_UCSC_CPG_ISLAND=cgi.relation,
							CHROMOSOME=ann[["Chromosome"]],
							POSITION=ann[["Start"]],
							row.names=rownames(ann))
					betas.norm<-do.call(wm.method, list(Mmatrix, Umatrix,
									da=ann.tost, pn=dpval(object)))
					object<-as(object, "RnBeadSet")
					betas.new<-meth(object)
					betas.new[match(rownames(betas.norm), probe.names),]<-betas.norm
					if(object@status$disk.dump){
						object@meth.sites<-convert.to.ff.matrix.tmp(betas.new)
					}else{
						object@meth.sites<-betas.new
					}
					rm(betas.norm, betas.new)

				}else{

					if(wm.method %in% c("nasen", "naten", "nanet")){
						list.norm<-do.call(wm.method, list(Mmatrix, Umatrix,
										ann[,"Design"], ret2=TRUE))
					}else{
						list.norm<-do.call(wm.method, list(Mmatrix, Umatrix,
										ann[,"Design"], ret2=TRUE, roco=pheno(object)[,grep("Sentrix[ |_]Position", colnames(pheno(object)))]))
					}
					object@M[,]<-list.norm[[1]][1:nrow(meth(object)),]
					object@U[,]<-list.norm[[2]][1:nrow(meth(object)),]
					rm(list.norm)
					#update.meth(object)
					if(object@status$disk.dump){
						object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))
					}else{
						object@meth.sites<-beta.value(object@M[,], object@U[,])
					}
				}
			}
		}

	}else if(method == "minfi.funnorm"){

		rnb.require("minfi")
		rnb.require("IlluminaHumanMethylation450kmanifest")
		rnb.require("IlluminaHumanMethylation450kanno.ilmn12.hg19")
		if(inherits(object,"MethyLumiSet") && (is.null(methylated(object))||is.null(unmethylated(object)))) {
			rnb.error("Invalid value for object; missing intensity information")
		}

		if(inherits(object,"MethyLumiSet")){

			if(nrow(methylated(object))<sum(elementNROWS(rnb.get.annotation("probes450"))) ||
					sum(is.na(methylated(object)))>0 || sum(is.na(unmethylated(object)))>0){
				rnb.warning("Funtional normalization is only supported for unfiltered data sets where intensity values are present for all probes. Skipping normalization")
				object<-as(object, "RnBeadRawSet")
				object@status$normalized<-"none"
				object@status$background<-bgcorr.method
				rnb.cleanMem()
				return(object)
			}
			intensities.by.channel<-methylumi.intensities.by.color(object)

		}else if(inherits(object,"RnBeadRawSet")){

			if(nrow(sites(object))<sum(elementNROWS(rnb.get.annotation("probes450"))) ||
					sum(is.na(M(object)))>0 || sum(is.na(U(object)))>0){
				rnb.warning("Funtional normalization is only supported for unfiltered data sets where intensity values are present for all probes. Skipping normalization")
				object@status$normalized<-"none"
				object@status$background<-bgcorr.method
				rnb.cleanMem()
				return(object)
			}
			intensities.by.channel<-intensities.by.color(object)
		}

		rg.set<-RGChannelSet(intensities.by.channel$Cy3, intensities.by.channel$Cy5)
		annotation(rg.set)<-c(array="IlluminaHumanMethylation450k", annotation=minfi:::.default.450k.annotation)
		suppressMessages({
					#sinkfile<-ifelse("Windows" %in% Sys.info(),"NUL", "/dev/null")
					#sink(sinkfile);
					#methyl.set<-preprocessFunnorm(rg.set);
					rg.set <- updateObject(rg.set)
					gmSet <- mapToGenome(rg.set)
					extractedData <- minfi:::.extractFromRGSet450k(rg.set)

					gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
					sex <- rep(1L, length(gmSet$predictedSex))
					sex[gmSet$predictedSex == "F"] <- 2L

					rm(rg.set)
					CN <- getCN(gmSet)
					methyl.set <- minfi:::.normalizeFunnorm450k(object = gmSet,
							extractedData = extractedData,
							sex = NULL, nPCs = 2, verbose = 0)

					#sink()
				})

		meth.minfi<-getMeth(methyl.set)
		umeth.minfi<-getUnmeth(methyl.set)

		meth.minfi<-pmax(meth.minfi, 0)
		umeth.minfi<-pmax(umeth.minfi, 0)

		if(inherits(object, "MethyLumiSet")){

			methylated(object)<-meth.minfi[match(rownames(meth.minfi), featureNames(object)),]#+methylated(object)[setdiff(featureNames(object), rownames(meth.minfi)),]
			umeth.minfi<-getUnmeth(methyl.set)
			unmethylated(object)<-umeth.minfi[match(rownames(umeth.minfi), featureNames(object)),]#+unmethylated(object)[setdiff(featureNames(object), rownames(umeth.minfi)),]
			rm(methyl.set, meth.minfi, umeth.minfi)
			betas(object)<-rbind(methylated(object)/(methylated(object)+unmethylated(object)),
					betas(object)[setdiff(featureNames(object), rownames(methylated(object))),])
			object<-as(object, "RnBeadSet")

		}else if(inherits(object, "RnBeadRawSet")){

			probe.ids<-rownames(annotation(object, add.names=TRUE))
			## rs probes are "lost" during SWAN

			meth.minfi<-meth.minfi[rownames(meth.minfi) %in% probe.ids,]
			umeth.minfi<-umeth.minfi[rownames(umeth.minfi) %in% probe.ids,]

			object@M[,][match(rownames(meth.minfi),probe.ids),]<-meth.minfi
			object@U[,][match(rownames(umeth.minfi),probe.ids),]<-umeth.minfi

			#update.meth(object)
			if(object@status$disk.dump){
				object@meth.sites<-convert.to.ff.matrix.tmp(beta.value(object@M[,], object@U[,]))
			}else{
				object@meth.sites<-beta.value(object@M[,], object@U[,])
			}
			rm(methyl.set, meth.minfi, umeth.minfi)
		}

	} else { # method == "none"
		if (inherits(object, "MethyLumiSet")) {
			object <- as(object, "RnBeadSet")
		}
	}

	## Display
	if (method != "none") {
		if (verbose) {
			rnb.status(c("Performed normalization with method", method))
		}
		rnb.cleanMem()
	}

	## Apply BMIQ as a secondary normalization method
	if (secondary.bmiq) {
		object <- rnb.execute.normalization.bmiq(object)
		if (verbose) {
			rnb.status("Performed normalization with method bmiq")
		}
		rnb.cleanMem()
	}

	## Update the dataset instance
	object@status$normalized <- method.to.set
	object@status$background <- bgcorr.method.to.set
	if (method != "none" || bgcorr.method !="none") {
		object <- updateRegionSummaries(object)
	}
	if (base::exists("inferred.covariates", inherits = FALSE)) {
		object@inferred.covariates <- inferred.covariates
	}

	object
}

########################################################################################################################

#' rnb.execute.normalization.bmiq
#'
#' Performs BMIQ normalization on the given dataset.
#' @param object Methylation dataset as an object of type \code{\linkS4class{MethyLumiSet}} or
#'               \code{\linkS4class{RnBSet}}.
#' @return The normalized dataset.
#'
#' @author Yassen Assenov
#' @noRd
rnb.execute.normalization.bmiq <- function(object) {

	## Extract methylation value matrix and probe design information
	if (inherits(object, "MethyLumiSet")) {
		m.data <- MethyLumiSet2RnBeadSet(object)
		beta.vals <- m.data$betas
		probe.design <- rnb.annotation2data.frame(rnb.get.annotation(object@target), add.names = TRUE)
		probe.design <- as.integer(probe.design[rownames(beta.vals), "Design"])
	} else {
		beta.vals <- object@meth.sites[, , drop = FALSE]
		probe.design <- as.integer(annotation(object)[, "Design"])
	}

	## Perform BMIQ
	samples.skipped <- integer()
	if (parallel.isEnabled() && ncol(beta.vals) > 1) {
		beta.names <- dimnames(beta.vals)
		beta.vals <- foreach(beta.v = as.data.frame(beta.vals), .combine = cbind, .packages = "RPMM",
							 .export = c("BMIQ", "betaEst2", "blc2"),
							 .noexport = c("bgcorr.method", "beta.names", "method", "object")) %dopar% {
			i <- which(!is.na(beta.v))
			p.design <- probe.design[i]
			if (all(tabulate(p.design, 2L) >= 50000L)) {
				x <- BMIQ(beta.v[i], p.design)
				if (!is.null(x)) {
					beta.v[i] <- x$all
				}
			}
			beta.v
		}
		dimnames(beta.vals) <- beta.names
		rm(beta.names)
	} else {
		for (j in 1:ncol(beta.vals)) {
			i <- which(!is.na(beta.vals[, j]))
			p.design <- probe.design[i]
			if (any(tabulate(p.design, 2L) < 50000L)) {
				## There are not enough probes of types I and/or II
				samples.skipped <- c(samples.skipped, j)
				rnb.status(c("Skipped sample", j))
			} else {
				x <- BMIQ(beta.vals[i, j], p.design)
				if (is.null(x)) {
					samples.skipped <- c(samples.skipped, j)
					rnb.status(c("Could not normalize sample", j))
				} else {
					beta.vals[i, j] <- x$all
					rnb.status(c("Normalized sample", j))
				}
			}
		}
		suppressWarnings(rm(j, i, p.design, x))
	}
	if (length(samples.skipped) != 0) {
		## Some samples were skipped due to not enough observations or unsuccessful fits
		rnb.warning(c("The following samples were not normalized:", paste(samples.skipped, collapse = ", ")))
	}
	rm(probe.design, samples.skipped)

	## Construct the resulting dataset
	if (inherits(object, "MethyLumiSet")) {
		object <- new("RnBeadSet", pheno=m.data$pheno, betas=beta.vals, p.values=m.data$p.values,
					  bead.counts=m.data$bead.counts)
		if ("qc" %in% names(m.data)) {
			qc(object) <- m.data[["qc"]]
		}
	} else {
		if (rnb.getOption("disk.dump.big.matrices")) {
			object@meth.sites <- convert.to.ff.matrix.tmp(beta.vals)
		} else {
			object@meth.sites <- beta.vals
		}
		for (region.type in rnb.region.types.for.analysis(object@assembly)) {
			object <- summarize.regions(object, region.type)
		}
	}
	object
}

########################################################################################################################
########################################################################################################################

#' rnb.section.normalization.shifts
#'
#' Generates a report sub-section on methylation beta values before and after correction.
#'
#' @param report       The report to contain the generated section.
#' @param betas.before Vector of raw beta values (before correction).
#' @param shifts       Vector of methylation value shifts (corrections). The \code{i}-th element in this vector must
#'                     correspond to the shift of the \code{i}-th element in \code{betas.before}.
#' @param bgcorr       Flag indicating if background subtraction was performed.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.normalization.shifts <- function(report, betas.before, shifts, bgcorr) {

	## Calculate frequencies of observed shifts and beta values
	shift.break.width <- 0.01
	shift.max <- ceiling(max(abs(shifts)) * 10) / 10
	shift.break.max <- shift.max + shift.break.width / 2
	shift.breaks <- seq(-shift.break.max, shift.break.max, by = shift.break.width)
	shift.mids <- seq(-shift.max, shift.max, length.out = length(shift.breaks) - 1)
	beta.break.width <- 0.025
	beta.breaks <- seq(0, 1, by = beta.break.width)
	beta.mids <- (beta.breaks[-length(beta.breaks)] + beta.breaks[-1]) / 2
	i.sb <- ceiling(betas.before / beta.break.width)
	i.sb[i.sb == 0] <- 1L
	i.sb <- as.integer(ceiling((shifts - shift.breaks[1]) / shift.break.width) - 1) * length(beta.mids) + i.sb
	i.sb <- tabulate(i.sb, length(shift.mids) * length(beta.mids))
	i.sb <- data.frame(
		xmin = rep(beta.breaks[-length(beta.breaks)], length(shift.mids)),
		xmax = rep(beta.breaks[-1], length(shift.mids)),
		ymin = rep(shift.breaks[-length(shift.breaks)], each = length(beta.mids)),
		ymax = rep(shift.breaks[-1], each = length(beta.mids)),
		freq = i.sb)

	## Plot histogram of shifts (differences)
	dframe <- data.frame(x = shift.mids, y = tapply(i.sb$freq, i.sb$ymin, sum))
	pp <- ggplot(dframe, aes_string(x = "x", y = "y")) + ggplot2::geom_bar(stat = "identity", width = 0.01) +
		scale_x_continuous(limits = c(-shift.break.max, shift.break.max), expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0)) + labs(x = 'Shift', y = 'Frequency')
	rplot <- createReportPlot("correction_shifts", report, width = 7, height = 5.2)
	print(pp)
	off(rplot)
	txt <- c("The next figure gives an idea of the magnitude of the correction by showing the distribution of shifts, ",
		"i.e. degrees of modification of the raw methylation values.")
	rnb.add.paragraph(report, txt)
	txt <- "Histogram of observed magnitude of &beta; value correction."
	report <- rnb.add.figure(report, txt, rplot)
	logger.status("Added histogram of observed beta shifts (magnitude of correction)")
	rm(shift.break.width, shift.max, shift.breaks, shift.mids, beta.break.width, beta.breaks, beta.mids)
	rm(dframe, pp)

	## Plot 2D histogram of beta value (before normalization) and shift
	pp <- ggplot(i.sb, aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax", fill = "freq")) +
		labs(x = expression(beta), y = "Shift", fill = "frequency") + geom_rect() +
		scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
		scale_fill_gradient(low = "#FFFFFF", high = "#000000") +
		scale_y_continuous(limits = c(-shift.break.max, shift.break.max), expand = c(0, 0)) +
		theme(legend.justification = c(0, 0.5), legend.position = c(1, 0.5)) +
		theme(panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank()) +
		theme(plot.margin = unit(c(0.1, 1.1, 0.1, 0.1), "in"))
	rplot <- createReportPlot("betas_shifts", report, width = 7.2, height = 6.2)
	print(pp)
	off(rplot)
	txt <- c("The figure below gives a more detailed view. This color-coded 2D histogram shows the uncorrected &beta; ",
		"values and their respective shifts after performing the normalization procedure.")
	rnb.add.paragraph(report, txt)
	txt <- "2D histogram showing the raw &beta; values and the magnitude of the corrections."
	report <- rnb.add.figure(report, txt, rplot)
	logger.status("Added 2D histogram of observed beta values and shifts")

	return(report)
}

########################################################################################################################

#' rnb.section.normalize.regions
#'
#' Adds a section summarizing the requested regions.
#'
#' @param report  Report on loading to contain the newly constructed section.
#' @param rnb.set Methylation dataset to be analyzed.
#' @param regions Non-empty \code{character} vector of region names.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.normalize.regions <- function(report, rnb.set, regions) {
	msg <- function(txt, e.class = "disabled") {
		paste0('<span class="', e.class, '">', txt, '</span>')
	}
	reg.counts <- sapply(regions, function(region) {
			tryCatch(as.character(nrow(meth(rnb.set, region))), error = function(e) { msg("not supported") })
		}
	)
	get.r.description <- function(reg) {
		result <- attr(rnb.get.annotation(reg, assembly = assembly(rnb.set)), "description")
		if (is.null(result) || identical("", result)) {
			result <- msg("n.a.")
		} else {
			result <- paste0('<p style="font-weight:normal;text-align:left">', result, '</p>')
		}
		result
	}
	reg.descriptions <- sapply(regions, function(region) {
			tryCatch(get.r.description(region), error = function(e) { msg("not imported", "outdated") })
		}
	)
	table.statistics <- data.frame(
		"Annotation" = regions, "Description" = reg.descriptions, "Regions in the Dataset" = reg.counts,
		check.names = FALSE, stringsAsFactors = FALSE)
	txt <- c(ifelse(rnb.getOption("analyze.sites"), "In addition to CpG sites, there", "There"),
		ifelse(length(regions) == 1, " is one set", paste(" are", length(regions), "sets")),
		" of genomic regions to be covered in the analysis. The table below gives a summary of these annotations.")
	report <- rnb.add.section(report, "Region Annotations", txt, level = 2)
	table.header <- c("<colgroup>", paste0('\t<col width="', c(210, 420, 150), 'px" />'), "</colgroup>")
	rnb.add.table(report, table.statistics, row.names = FALSE, thead = table.header)

	return(report)
}

#######################################################################################################################

#' rnb.section.normalization
#'
#' Performs normalization of the loaded methylation dataset.
#'
#' @param report    Report to contain the normalization section. This must be an object of type
#'                  \code{\linkS4class{Report}}.
#' @param rnb.set   Normalized methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param betas.raw Matrix of methylation beta values for CpG sites, as they would appear without (prior to)
#'                  normalization. If this parameter is \code{NULL}, the effect of normalization on beta values is
#'                  not examined and visualized.
#' @return The modified report.
#'
#' @details
#' If specified \code{betas.raw}, must be a matrix of methylation values with dimensions identical to the one
#' encapsulated in \code{rnb.set}. The former matrix is expected to store beta values before normalization,
#' whereas the latter one - after normalization. If the provided matrices contain sufficient amounts of non-missing
#' values, this function creates figures that compare the two distributions of values and examines the magnitute of
#' modifications.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.section.normalization <- function(report, rnb.set, betas.raw = NULL) {
	if (!inherits(report, "Report")) {
		stop("Invalid value for report")
	}
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!is.null(betas.raw)) {
		if (!(is.matrix(betas.raw) && (is.double(betas.raw) || is.integer(betas.raw)))) {
			stop("invalid value for betas.raw")
		}
		betas.after <- meth(rnb.set, row.names = TRUE)
		if (!(nrow(betas.raw) == nrow(betas.after) && ncol(betas.raw) == ncol(betas.after) &&
			setequal(rownames(betas.raw), rownames(betas.after)) &&
			setequal(colnames(betas.raw), colnames(betas.after)))) {
			stop("incompatible values for rnb.set and betas.raw")
		}
		betas.after <- betas.after[rownames(betas.raw), colnames(betas.raw)]
	}

	txt.methylumi <- "<a href=\"http://www.bioconductor.org/packages/release/bioc/html/methylumi.html\">methylumi</a>"

	bgcorr.method <- rnb.set@status$background
	bgcorr <- (bgcorr.method != "none")
	bgcorr.method <- gsub("methylumi\\.","", bgcorr.method)
	if (bgcorr) {
		refText <- c("Triche, T.J. Jr, Weisenberger, D.J., Van Den Berg, D., Laird, P.W., and Siegmund, K.D. (2013)",
			"Low-level processing of Illumina Infinium DNA Methylation BeadArrays. ",
			"<i>Nucleic Acids Research</i> <b>41</b>(7), e90")
		report <- rnb.add.reference(report, refText)
		txt <- c("The background was subtracted using the ", txt.methylumi, " package (method \"",
			bgcorr.method, "\") ", rnb.get.reference(report, refText), ".")
		txt.methylumi <- "methylumi"
	} else {
		txt <- ""
	}
	method <- rnb.set@status$normalized
	if (method == "illumina") {
		txt <- c(txt, "The signal intensity values were normalized using Illumina's recommended normalization method, ",
			"as implemented in the ", txt.methylumi, " package.")
	} else if (method == "swan") {
		txt <- c(txt, "The signal intensity values were normalized using the SWAN normalization method, as ",
			"implemented in the ",
			"<a href=\"http://www.bioconductor.org/packages/release/bioc/html/minfi.html\">minfi</a> package.")
	} else if (method == "bmiq") {
		refText <- c("Teschendorff, A. E., Marabita, F., Lechner, M., Bartlett, T., Tegner, J., Gomez-Cabrero, D. and ",
			"Beck, S. (2013) A beta-mixture quantile normalization method for correcting probe design bias in ",
			"Illumina Infinium 450 k DNA methylation data. <i>Bioinformatics</i> <b>29</b>(2), 189-196")
		report <- rnb.add.reference(report, refText)
		txt <- c(txt, "The methylation &beta; values were normalized using the BMIQ normalization method ",
			rnb.get.reference(report, refText), ".")
	} else if (grepl("wm\\.", method)){
		refText <- c("Pidsley, R., Wong, C.,  Volta, M., Lunnon, K., Mill, J., and Schalkwyk, L. (2013)",
				"A data-driven approach to preprocessing Illumina 450K methylation array data."," <i>BMC Genomics</i>, <b>14<b>(1), 293")
		report <- rnb.add.reference(report, refText)
		wm.method<-gsub("wm\\.", "", method)
		txt<-c(txt, sprintf("The data was normalized using method %s from %s.", wm.method, rnb.get.reference(report,refText)))
	}else if (method == "minfi.funnorm"){
		refText<-c("Fortin, J., Labbe, A., Lemire, M., Zanke, B. W., Hudson, T. J., Fertig, E. J., Greenwood, C.M.T., ",
			"Hansen, K. D. (2014) Functional normalization of 450k methylation array data improves replication in ",
			"large cancer studies. <i>BioRxiv</i>.")
		report <- rnb.add.reference(report, refText)
		txt<-c(txt, sprintf("The data was normalized using the functional normalization method from %s.", rnb.get.reference(report,refText)))

	}else{ # method == "none"
		txt <- c(txt, "The measurements in this dataset were not normalized after ",
			ifelse(txt == "", "loading", "background subtraction"), ".")
		if(rnb.getOption("normalization.method")=="minfi.funnorm"){
			if(nrow(sites(rnb.set))<sum(elementNROWS(rnb.get.annotation("probes450")))){
				txt<-c(txt, "The desired normalization method \"minfi.funnorm\" was not applied since the data set did not
					contain information for all probes. This is most likely because of the preceding quality filtering steps.",
					"In order to apply minfi.funnorm in the full RnBeads analysis you should disable SNP filtering and Greedycut.",
					"Otherwise, you apply rnb.execute.normalization() function to an unfiltered RnBead(Raw)Set after loading and start
					the full pipeline with the returned object as input.")
			}else if(sum(is.na(M(rnb.set)))>0 || sum(is.na(U(rnb.set)))>0){
				txt<-c(txt, "The desired normalization method \"minfi.funnorm\" was not applied because intensity information contains
						missing values.", "Try to apply rnb.execute.normalization() function to an unfiltered RnBead(Raw)Set after loading and start
								the full pipeline with the returned object as input.")
			}
		}
	}
	report <- rnb.add.section(report, "Normalization", txt)

	## Add comparison of betas before and after correction
	if ((!is.null(betas.raw)) && (bgcorr || method != "none")) {
		betas.raw <- as.vector(betas.raw)
		betas.after <- as.vector(betas.after)
		min.observations <- 501L

		i.before <- !is.na(betas.raw); i.after <- !is.na(betas.after)
		i.both <- which(i.before & i.after)
		i.before <- which(i.before); i.after <- which(i.after)

		section.title <- "Effect of Correction"
		if (length(i.before) >= min.observations && length(i.after) >= min.observations) {
			txt <- c("This section shows the influence of the applied normalization procedure on CpG methylation ",
				"values. The following figure compares the distributions of the &beta; values before and after ",
				"performing normalization.")
			report <- rnb.add.section(report, section.title, txt, level = 2)

			## Distributions before and after normalization
			beta.values <- list("Before correction" = betas.raw[i.before], "After correction" = betas.after[i.after])
			report.plots <- rnb.plot.beta.comparison(beta.values, "correction_comparison", report, min.observations)
			txt <- "Comparison of &beta; values before and after correction."
			txt <- c(txt, add.text.subsampling(attr(report.plots, "subsampled"), paste("betas", names(beta.values))))
			setting.names <- list("Plot type" =
				c("density" = "density estimation", "histogram" = "histograms", "qq" = "quantile-quantile plot"))
			report <- rnb.add.figure(report, txt, report.plots, setting.names)
			rm(i.before, i.after, beta.values, report.plots, setting.names)
			rnb.cleanMem()
			logger.status("Added comparison between non-normalized and normalized beta values")

			if (rnb.getOption("normalization.plot.shifts")) {
				if (length(i.both) != 0) {
					betas.raw <- betas.raw[i.both]
					betas.after <- betas.after[i.both] - betas.raw
					report <- rnb.section.normalization.shifts(report, betas.raw, betas.after, bgcorr)
				} else {
					txt <- c("No plots on &beta; value shifts were created because there are not enough data to ",
						"visualize.")
					rnb.add.paragraph(report, txt)
				}
			}
		} else {
			txt <- c("The effects of the applied procedure on CpG methylation values is not summarized because there ",
				"are not enough &beta; values to accurately present the observed distributions before and after ",
				"normalization.")
			report <- rnb.add.section(report, section.title, txt)
		}
		rm(betas.raw, betas.after)
	}

	if (!inherits(rnb.set, "RnBeadSet")) {
		return(report)
	}

	## Add point-and-whisker plots of mean methylation for every slide
	report.plots <- rnb.plot.sentrix.distributions(rnb.set, report = report, width = 6, height = 6)
	if (is.null(report.plots)) {
		txt <- "Sample average methylation cannot be visualized because "
		if (length(samples(rnb.set)) == 0) {
			txt <- c(txt, "the normalized dataset is empty.")
		} else if (all(is.na(meth(rnb.set)))) {
			txt <- c(txt, "the dataset contains no valid methylation beta values.")
		} else {
			txt <- c(txt, "no valid Sentrix ID and Sentrix Position information could be extracted from the ",
				"sample annotation table.")
		}
	} else {
		txt <- "The following figure visualizes the average methylation per sample. Samples are grouped by slide."
	}
	report <- rnb.add.section(report, "Sample Mean Methylations", txt, level = 2)
	if (!is.null(report.plots)) {
		description <- "Point-and-whisker plot showing mean and standard deviation among all beta values in a sample."
		if (inherits(report.plots, "ReportPlot")) {
			setting.names <- list()
		} else { # is.list(report.plots)
			setting.names <- list("Slide number" = names(report.plots))
			names(setting.names[[1]]) <- 1:length(report.plots)
		}
		report <- rnb.add.figure(report, description, report.plots, setting.names)
	}

	r.types <- rnb.getOption("region.types")
	if (is.null(r.types)) {
		r.types <- summarized.regions(rnb.set)
	}
	if (length(r.types) != 0) {
		report <- rnb.section.normalize.regions(report, rnb.set, r.types)
	}

	return(report)
}

#######################################################################################################################

#' rnb.step.normalization
#'
#' Performs normalization of the loaded HumanMethylation450 dataset and adds a corresponding section to the given
#' report.
#'
#' @param object Methylation dataset as an object of type inheriting \code{\linkS4class{MethyLumiSet}} or
#'               \code{\linkS4class{RnBSet}}.
#' @param report Report to contain the normalization section. This must be an object of type
#'               \code{\linkS4class{Report}}.
#' @param method Normalization method, must be one of \code{"illumina"}, \code{"swan"} or \code{"none"}. Note that
#'               the execution of method SWAN requires packages \pkg{minfi} and
#'               \pkg{IlluminaHumanMethylation450kmanifest}.
#'
#' @return List with two elements:
#'         \describe{
#'           \item{\code{dataset}}{Normalized dataset as an object of type \code{\linkS4class{RnBeadSet}}.}
#'           \item{\code{report}}{the modified report.}
#'         }
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.step.normalization<-function(object, report, method = rnb.getOption("normalization.method")){

	if(!(inherits(object,"MethyLumiSet") || inherits(object,"RnBSet"))) {
		stop("Invalid value for object; expected MethyLumiSet, RnBeadSet or RnBiseqSet")
	}
	if (!inherits(report, "Report")) {
		stop("Invalid value for report")
	}
	accepted <- .rnb.options[["accepted"]][["normalization.method"]]
	if (!(is.character(method) && length(method) == 1 && isTRUE(method %in% accepted))) {
		msg <- paste0('"', accepted, '"', collapse = ", ")
		stop(paste("invalid value for method; expected one of", msg))
	}

	logger.start("Normalization Procedure")


	if (method == "none") {
		betas.raw <- NULL
	} else if (inherits(object, "MethyLumiSet")) {
		betas.raw <- betas(object)
		pheno.table <- phenoData(object)@data
		id.column <- rnb.getOption("identifiers.column")
		if (!(is.null(pheno.table) || is.null(id.column))) {
			ids <- NA
			if (is.character(id.column) && (id.column %in% colnames(pheno.table))) {
				ids <- pheno.table[, id.column]
			} else if (1 <= id.column && id.column <= ncol(pheno.table)) {
				ids <- pheno.table[, id.column]
			}
			if (any(is.na(ids)) == FALSE && anyDuplicated(ids) == 0) {
				colnames(betas.raw) <- as.character(ids)
			}
			rm(ids)
		}
		rm(pheno.table, id.column)
		rnb.cleanMem()
	} else { # inherits(object, "RnBSet")
		betas.raw <- meth(object, row.names = TRUE)
	}
	#only destroy an object if the new object is not the same as the old object,
	#i.e. if normalization or background subtraction took place
	#if ((rnb.set@status$normalized != "none") || (rnb.set@status$background != "none")){
	#	rnb.call.destructor(object)
	#}

	if(object@status$disk.dump && rnb.getOption("enforce.destroy.disk.dumps")){
		object@status$discard.ff.matrices<-TRUE
	}
	object <- rnb.execute.normalization(object, method)
	rnb.cleanMem()
	if(isTRUE(object@status$discard.ff.matrices)){
		object@status$discard.ff.matrices<-NULL
	}
#	if(object@status$background!="none"){
#		logger.status(sprintf("Performed background correction with method \"%s\"", object@status$background))
#	}
	if(object@status$normalized!="none"){
		logger.status(sprintf("Performed normalization with method \"%s\"", object@status$normalized))
	}

	report<-rnb.section.normalization(report, object, betas.raw)
	logger.status("Added normalization section")

	logger.completed()
	return(list(dataset=object, report=report))
}

#######################################################################################################################

#' mean.imputation
#'
#' Performs mean imputation either for all samples (way=1) or for all CpGs (way=2).
#'
#'@param rnb.set Object containing the methylation information to be changed. Has to be of
#'               type \code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBiseqSet}} or a data matrix containing
#'               the methylation values.
#'@param way Should the sample-wise mean (1) or CpG-wise mean (2) be used to replace the missing value.
#'@return Modified rnb.set object.
#'
#'@author Michael Scherer
#'@noRd
mean.imputation <- function(rnb.set,way=1){
  if(inherits(rnb.set,"RnBSet")){
    methData <- meth(rnb.set)
  }else{
    methData <- as.matrix(rnb.set)
  }
  nas <- is.na(methData)
  if(any(apply(nas,way,all))){
    logger.warning(paste("There are", ifelse(way==1,"CpG sites","sampels"),"that have missing values in all",
                         ifelse(way==1,"samples","CpG sites"),", imputation not performed."))
    return(NULL)
  }
  means <- apply(methData,way,mean,na.rm=TRUE)
  has.nas <- which(apply(nas,way,any))
  if(way==1){
    for(i in has.nas){
      methData[i,nas[i,]] <- means[i]
    }
  }else if(way==2){
    for(i in has.nas){
        methData[nas[,i],i] <- means[i]
    }
  }
  return(methData)
}

#######################################################################################################################

#' imputation.low.memory.samples
#'
#' Performs (sample-wise) imputation in a low memory-footprint way with the specified method
#'
#'@param rnb.set Object containing the methylation information to be changed. Has to be of
#'               type \code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBiseqSet}} or a data matrix containing
#'               the methylation values.
#'@param method Method to be used, should be either 'mean' or 'median'
#'@return Modified rnb.set object.
#'
#'@author Michael Scherer
#'@noRd
imputation.low.memory.samples <- function(rnb.set,method=mean){
 for(i in 1:nsites(rnb.set)){
    cpg <- meth(rnb.set,i=i)
    if(any(is.na(cpg))){
      rnb.set@meth.sites[i,which(is.na(cpg))] <- method(cpg,na.rm=T)
    }
  }
  return(rnb.set)
}

#######################################################################################################################

#' imputation.low.memory.cpgs
#'
#' Performs (cpg-wise) imputation in a low memory-footprint way with the specified method
#'
#'@param rnb.set Object containing the methylation information to be changed. Has to be of
#'               type \code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBiseqSet}} or a data matrix containing
#'               the methylation values.
#'@param method Method to be used, should be either 'mean' or 'median'
#'@return Modified rnb.set object.
#'
#'@author Michael Scherer
#'@noRd
imputation.low.memory.cpgs <- function(rnb.set,method=mean){
  if(any(sampleMethApply(rnb.set,function(x){all(is.na(x))}))){
    logger.warning("There are samples that have missing values at all CpG sites, imputation not performed.")
    return(rnb.set)
  }
  for(i in 1:length(samples(rnb.set))){
    ss <- meth(rnb.set,j=i)
    if(any(is.na(ss))){
      rnb.set@meth.sites[which(is.na(ss)),i] <- method(ss,na.rm=T)
    }
  }
  return(rnb.set)
}


#######################################################################################################################

#' random.imputation
#'
#' Performs random imputation by replacing missing values with randomly selecting (with replacement) from the same CpG site in the samples
#' that do not contain missing values.
#'
#'@param rnb.set Object containing the methylation information to be changed. Has to be of
#'               type \code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBiseqSet}} or a data matrix containing
#'               the methylation values.
#'@return Modified rnb.set object.
#'
#'@author Michael Scherer
#'@noRd
random.imputation <- function(rnb.set){
  if(inherits(rnb.set,"RnBSet")){
    methData <- meth(rnb.set)
  }else{
    methData <- as.matrix(rnb.set)
  }
  nas <- is.na(methData)
  if(any(apply(nas,1,all))){
    logger.warning("There are CpG sites that have missing values in all samples, imputation not performed.")
    return(NULL)
  }
  has.nas <- which(apply(nas,1,any))
  for(i in has.nas){
    row <- methData[i,]
    without_nas <- row[!nas[i,]]
    replacement <- sample(without_nas,sum(nas[i,]),replace=TRUE)
    methData[i,nas[i,]] <- replacement
  }
  return(methData)
}
#######################################################################################################################

#' knn.imputation
#'
#' This function performns k nearest neighbors imputation by calling the \code{impute.knn} function from the
#' \pkg{impute} package.
#'
#'@param rnb.set Object containing the methylation information to be changed. Has to be of
#'               type \code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBiseqSet}} or a data matrix containing
#'               the methylation values.
#'@param k parameter defining the number of nearest neighbors from which the missing value should be inferred
#'@return Modified rnb.set object.
#'
#'@author Michael Scherer
#'@noRd
knn.imputation <- function(rnb.set,k=10){
  rnb.require('impute')
  if(inherits(rnb.set,"RnBSet")){
    methData <- meth(rnb.set)
  }else{
    methData <- as.matrix(rnb.set)
  }
  dummy <- capture.output(methData <- (impute.knn(methData,colmax=1,k=k))$data)
  rm(dummy)
  return(methData)
}

#######################################################################################################################

#' median.imputation
#'
#' Performs median imputation either for all samples (way=1) or for all CpGs (way=2).
#'
#'@param rnb.set Object containing the methylation information to be changed. Has to be of
#'               type \code{\linkS4class{RnBeadSet}} or \code{\linkS4class{RnBiseqSet}} or a data matrix containing
#'               the methylation values.
#'@param way Should the sample-wise median (1) or CpG-wise median (2) be used to replace the missing value.
#'@return Modified rnb.set object.
#'
#'@author Michael Scherer
#'@noRd
median.imputation <- function(rnb.set,way=1){
  if(inherits(rnb.set,"RnBSet")){
    methData <- meth(rnb.set)
  }else{
    methData <- as.matrix(rnb.set)
  }
  nas <- is.na(methData)
  if(any(apply(nas,way,all))){
    logger.warning(paste("There are", ifelse(way==1,"CpG sites","sampels"),"that have missing values in all",
                ifelse(way==1,"samples","CpG sites"),", imputation not performed."))
    return(NULL)
  }
  medians <- apply(methData,way,median,na.rm=TRUE)
  has.nas <- which(apply(nas,way,any))
  if(way==1){
    for(i in has.nas){
      methData[i,nas[i,]] <- medians[i]
    }
  }else if(way==2){
    for(i in has.nas){
      methData[nas[,i],i] <- medians[i]
    }
  }
  return(methData)
}

#######################################################################################################################

#' rnb.execute.imputation
#'
#' Removes missing methylation values in the methylation matrix of the given object
#'
#' @param rnb.set Dataset object inheriting from \code{\linkS4class{RnBSet}}.
#' @param method Imputation method to be used, must be one of \code{"mean.cpgs"}, \code{"mean.samples"},
#'                \code{"random"}, \code{"knn"}, \code{"median.cpgs"}, \code{"median.samples"}, or \code{"none"}.
#' @param update.ff flag indicating if the disk based matrices should be updated. Should be set to FALSE, if methylation
#'                  matrix should only temporarly be changed. If this value is FALSE, the region level methylation values
#'                  are not updated and only the site-wise matrix is changed temporarly.
#' @param ... Optional arguments passed to knn.imputation
#' @return The modified rnb.set object without missing methylation values.
#' 
#' @details Imputes missing values by applying on the following methods:
#'          \describe{
#'            \item{mean.cpgs:}{missing values are inferred as the average methylation value from all other 
#'            (non-mising) CpGs in this sample}
#'            \item{mean.samples:}{missing values are inferred as the average methylation value from all other 
#'            (non-mising) values at this CpG sites in all other samples}
#'            \item{random:}{missing values are inferred by randomly selecting a (non-missing) methylation value
#'            from any other sample at this CpG site}
#'            \item{knn:}{missing values are inferred by k-nearest neighbors imputation (see \pkg{impute})}
#'            \item{median.cpgs:}{missing values are inferred as the median methylation value from all other 
#'            (non-mising) CpGs in this sample}
#'            \item{median.samples:}{missing values are inferred as the median methylation value from all other 
#'            (non-mising) values at this CpG sites in all other samples}
#'            \item{none:}{imputation should not be performed}
#'          }
#'
#' @author Michael Scherer
#'
#' @export
rnb.execute.imputation <- function(rnb.set,method=rnb.getOption("imputation.method"),update.ff=TRUE,...){
  if(!(inherits(rnb.set,"RnBSet")||is.matrix(rnb.set)||is.data.frame(rnb.set))){
    stop("Invalid value for input object, has to be of type RnBeadSet or RnBiseqSet")
  }
  if(inherits(rnb.set,"RnBSet")&&isImputed(rnb.set)){
    logger.info("RnBSet already imputed, imputation skipped.")
    return(rnb.set)
  }
  if(inherits(rnb.set,"RnBSet") && (length(samples(rnb.set))==1) && !(method%in%c("median.cpgs","mean.cpgs"))){
    logger.info(paste("Imputation cannot be employed with method",method,"median.cpgs used"))
    method <- "median.cpgs"
  }
  if(!(method%in%c('mean.cpgs','mean.samples','random','knn','median.cpgs','median.samples'))){
    if(method=='none'){
      if(inherits(rnb.set,"RnBeadSet")){
        logger.info("No imputation method selected, 'knn' method used.")
        rnb.options(imputation.method='knn')
        method = 'knn'
      }else if (inherits(rnb.set,"RnBiseqSet")||is.matrix(rnb.set)||is.data.frame(rnb.set)){
        logger.info("No imputation method selected, 'mean.samples' method used.")
        rnb.options(imputation.method='mean.samples')
        method = 'mean.samples'
      }
    }else{
      stop("Invalid option for imputation method, has to be one of 'mean.cpgs','mean.samples','random','knn'")
    }
  }
  if(inherits(rnb.set,"RnBiseqSet") && rnb.getOption("imputation.method")=="knn"){
    rnb.options("imputation.method"="mean.samples")
    method = "mean.samples"
    logger.info("Knn imputation not applicable to sequencing data sets, switched to 'mean.samples' method")
    rnb.options("imputation.method"=method)
  }
  if(rnb.getOption("enforce.memory.management")){
    if(!method%in%c("mean.cpgs","median.cpgs")){
      logger.info(sprintf("Low memory imputation not compatible with method %s, switched to mean.cpgs",method))
      method <- "mean.cpgs"
      rnb.options("imputation.method"=method)
    }
    if(method=="mean.cpgs"){
      if(is.data.frame(rnb.set)||is.matrix(rnb.set)){
        logger.info("No low memory footprint imputation available for matrix")
        meth.data <- mean.imputation(rnb.set,way=2)
        return(meth.data)
      }
      logger.start("Low memory footprint version of imputation")
      rnb.set <- imputation.low.memory.cpgs(rnb.set,mean)
      logger.completed()
    }
    if(method=="median.cpgs"){
      if(is.data.frame(rnb.set)||is.matrix(rnb.set)){
        logger.info("No low memory footprint imputation available for matrix")
        meth.data <- median.imputation(rnb.set,way=2)
        return(meth.data)
      }
      logger.start("Low memory footprint version of imputation")
      rnb.set <- imputation.low.memory.cpgs(rnb.set,median)
      logger.completed()
    }
    if(inherits(rnb.set,"RnBSet")){
      rnb.set <- updateRegionSummaries(rnb.set)
      rnb.set@imputed <- TRUE
    }
    return(rnb.set)
  }
  logger.start(sprintf("Imputation procedure %s ",method))
  if(method=='mean.cpgs'){
    meth.data <- mean.imputation(rnb.set,2)
  }
  if(method=='mean.samples'){
    meth.data <- mean.imputation(rnb.set,1)
  }
  if(method=='random'){
    meth.data <- random.imputation(rnb.set)
  }
  if(method=='knn'){
    meth.data <- knn.imputation(rnb.set,...)
  }
  if(method=='median.cpgs'){
    meth.data <- median.imputation(rnb.set,2)
  }
  if(method=='median.samples'){
    meth.data <- median.imputation(rnb.set,1)
  }
  if(!is.null(meth.data)){
    if(inherits(rnb.set,"RnBSet")){
      if(update.ff){
        rnb.set <- updateMethylationSites(rnb.set,meth.data)
        rnb.set@imputed <- TRUE
      }else{
        rnb.set@meth.sites <- meth.data
      }
    }else{
      rnb.set <- meth.data
    }
  }
  logger.completed()
  return(rnb.set)
}

#######################################################################################################################

#' rnb.section.imputation
#'
#' Adds information and plots describing the effect of imputation on the dataset
#'
#' @param report Report to which the information should be added (\code{\linkS4class{Report}}).
#' @param rnb.set Dataset on which imputation was performed (\code{\linkS4class{RnBSet}}).
#' @param old.values Methylation matrix before imputation procedure.
#'
#' @return The modified report object
#' @author Michael Scherer
#' @noRd
rnb.section.imputation <- function(report,rnb.set,old.values){
  if(rnb.getOption('imputation.method')=='none'){
    return(report)
  }
  txt <- c("Imputation was performed ")
  if(rnb.getOption('imputation.method')=='mean.samples'){
    txt <- c(txt,"by calculating the mean methylation level for each CpG site across all samples",
             " and replacing all missing values for this CpG site in individual samples with the mean across all samples.\n")
  }
  if(rnb.getOption('imputation.method')=='mean.cpgs'){
    txt <- c(txt,"by calculating the mean methylation level for each sample across all CpG sites",
             " and replacing all missing values for this sample at an individual CpG site with the mean across all CpGs in the sample.\n")
  }
  if(rnb.getOption('imputation.method')=='random'){
    txt <- c(txt,"by randomly selecting values from non-missing CpGs in other samples.",
             " Sampling was done with replacement such that the number of missing values can be larger as the number of existing.\n")
  }
  if(rnb.getOption('imputation.method')=='knn'){
    txt <- c(txt,"by k nearest neighbors imputation from the <a href=https://bioconductor.org/packages/release/bioc/html/impute.html>impute</a> package. Briefly, missing values were replaced",
             " by the average of the methylation values at that CpG in the closest samples. Closeness was defined with",
             " the euclidean distance.\n")
  }
  if(rnb.getOption('imputation.method')=='median.samples'){
    txt <- c(txt,"by calculating the median methylation level for each CpG site across all samples",
             " and replacing all missing values for this CpG site in individual samples with the median across all samples.\n")
  }
  if(rnb.getOption('imputation.method')=='median.cpgs'){
    txt <- c(txt,"by calculating the median methylation level for each sample across all CpG sites",
             " and replacing all missing values for this sample at an individual CpG site with the median across all CpGs in the sample.\n")
  }
  new.values <- meth(rnb.set)
  missing.values <- apply(old.values,2,function(x)sum(is.na(x)))
  median <- median(missing.values)
  txt <- c(txt," Imputation replaced a median of ",median," missing values per sample by estimations.")
  report <- rnb.add.section(report,'Imputation',txt)
  if(sum(is.na(new.values))==sum(is.na(old.values))){
    txt <- "Imputation failed, see logger output for more information."
    report <- rnb.add.paragraph(report,txt)
    return(report)
  }
  old.values <- melt(old.values)
  old.values <- old.values$value
  old.values <- na.omit(old.values)
  new.values <- melt(new.values)
  new.values <- new.values$value
  min.observations <- 501L
  beta.values <- list("NAs removed" = old.values, "After imputation" = new.values)
  report.plot <- rnb.plot.beta.comparison(beta.values, "imputation_comparison", report, min.observations)
  setting.names <- list("Plot type" =
                          c("density" = "density estimation", "histogram" = "histograms", "qq" = "quantile-quantile plot"))
  descr <- 'Visulisation of the effect of Imputation on the methylation distribution.'
  report <- rnb.add.figure(report,descr,report.plot,setting.names)
  return(report)
}
#######################################################################################################################

#' rnb.step.imputation
#'
#' Imputes missing values from the given dataset with the specified method
#'
#' @param rnb.set Dataset object inheriting from \code{\linkS4class{RnBSet}}.
#' @param report Report object to contain imputation information. Has to be of class
#'                \code{\linkS4class{Report}}.
#' @param method Imputation method to be used, must be one of \code{"mean.cpgs"}, \code{"mean.samples"},
#'                \code{"random"}, \code{"knn"} or \code{"none"}.
#'
#' @return List with two elements:
#'         \describe{
#'           \item{\code{dataset}}{Dataset as an object of type \code{\linkS4class{RnBSet}}
#'                 containing the imputed methylation matrix.}
#'           \item{\code{report}}{the modified report.}
#'         }
#'
#' @author Michael Scherer
#' @noRd
rnb.step.imputation <- function(rnb.set, report, method=rnb.getOption("imputation.method")){
  if(!(inherits(rnb.set,"RnBSet"))) {
    stop("Invalid value for object; RnBeadSet or RnBiseqSet required")
  }
  if (!inherits(report, "Report")) {
    stop("Invalid value for report")
  }
  old.data <- meth(rnb.set)
  if(sum(is.na(old.data))==0){
    logger.info("No missing values present, imputation skipped")
    return(list(dataset=rnb.set,report=report))
  }
  if(!(method%in%c('none'))){
    rnb.set <- rnb.execute.imputation(rnb.set)
  }
  report <- rnb.section.imputation(report,rnb.set,old.data)
  return(list(dataset=rnb.set,report=report))
}

#######################################################################################################################
