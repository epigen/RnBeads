########################################################################################################################
## CNV.R
## created: 2013-06-11
## creator: Pavlo Lutsik, edited by Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## Copy number variation based on HumanMethylation450 data.
########################################################################################################################

## utility functions

#' rnb.get.cv.annotations
#' 
#' This function returns reference profiles for CNV computation. Only 450k data is
#' supported in the moment.
#' 
#' @param platform The platform of the refernce data set to be returned. Acceptible
#'                  values are \code{"probes450}.
#' @return A list with three elements stored in the \code{\link{RnBeads.hg19}} package.
#'          \itemize{
#'             \item{\code{seq.info.hg19}} A Seqinfo object containing general sequence information.
#'             \item{\code{ref.intensity}} Reference intensity values used to compute CNVs
#'             \item{\code{bac.array.cpg}} CpG identifieres of the entries in \code{ref.intensity}
#'          }
#' @details The reference profiles have been obtained as the median signal intensities of twins from two twin studies:
#'    GSE85647 and GSE100940. 
#' @noRd
#' @author Pavlo Lutsik, with modifications by Michael Scherer         
rnb.get.cnv.annotations<-function(platform="probes450"){
  if(!platform %in% c("probes450")){
    logger.warning("Reference dataset is only available for HumanMethylation450 chips.")
    return(NULL)
  }
  rnb.require("RnBeads.hg19")
  ref.loc <- system.file("extdata/cnv.reference.RDS", package="RnBeads.hg19")
  if(file.exists(ref.loc)){
    cnv.reference.data <- readRDS(ref.loc)
  }else{
    logger.warning("Reference dataset for CNV estimation not available. Please update RnBeads.hg19.")
    cnv.reference.data <- NULL
  }
	return(cnv.reference.data)
}

#######################################################################################################################

#' getGLADProfiles
#' 
#' This function computes and returns GLAD profiles.
#' 
#' @param rnb.set An object of type \code{\link{RnBeadRawSet}} containing intensity information for
#'                 computing the GLAD profiles.
#' @return GLAD profiles as produced by \code{\link{daglad}}, with names corresponding
#'          to the samples of \code{rnb.set}
#' @noRd
#' @author Pavlo Lutsik, with modifications by Michael Scherer         
getGLADProfiles<-function(rnb.set,refbased=TRUE){
	
	if(!inherits(rnb.set,"RnBeadRawSet")){
		stop("RnBeadRawSet object expected")
	}
	
  annot <- annotation(rnb.set)
  target <- rnb.set@target
  if(target %in% c("probes27","probes450","probesEPIC")){
    target <- "probes450"
  }
  cnv.reference.data<-rnb.get.cnv.annotations(target)
  if(rnb.getOption("qc.cnv.refbased") && !is.null(cnv.reference.data)){
	  matchi <- rownames(annot) %in% cnv.reference.data$bac.array.cps
	  annot<-annot[matchi,]
	  ref.int <- cnv.reference.data$ref.intensity[cnv.reference.data$bac.array.cps %in% row.names(annot)]
  
  	I<-M(rnb.set)+U(rnb.set)
  	I<-I[matchi,]
  }else{
    I<-M(rnb.set)+U(rnb.set)
    ref.int <- rowMeans(I)
  }
	
	sample.names<-samples(rnb.set)
	
	primary.data<-lapply(sample.names, function(cl){
	  log2(I[,cl]/ref.int)
	})
	names(primary.data)<-sample.names
	
	if(parallel.isEnabled()){	
  	ncores <- parallel.getNumWorkers()
  	cgh.profiles <- mclapply(primary.data, function(df){
  	  df.complete <- data.frame(
  	    Clone=row.names(annot),
  	    PosOrder=1:nrow(annot), 
  	    LogRatio=df,
  	    PosBase=annot[,"Start"],
  	    Chromosome=annot[,"Chromosome"],
  	    BAC=row.names(annot)
  	  )
  	  dummy <- capture.output(res <- as.profileCGH(df.complete, value=3))
  	  return(res)
  	}
  	                         ,mc.cores = ncores)
  	names(cgh.profiles)<-sample.names
  	glad.profiles<-mclapply(cgh.profiles, function(cl){
  	  dummy <- capture.output(res <- daglad(cl,mediancenter=TRUE,param=c(d=9)))
  	  return(res)
  	},mc.cores = ncores)
  	names(glad.profiles)<-sample.names
	}else{
	  cgh.profiles <- lapply(primary.data, function(df){
	    df.complete <- data.frame(
	      Clone=row.names(annot),
	      PosOrder=1:nrow(annot), 
	      LogRatio=df,
	      PosBase=annot[,"Start"],
	      Chromosome=annot[,"Chromosome"],
	      BAC=row.names(annot)
	    )
	    dummy <- capture.output(res <- as.profileCGH(df.complete, value=3))
	    return(res)
	  })
	  names(cgh.profiles)<-sample.names
	  glad.profiles<-lapply(cgh.profiles, function(cl){
	    dummy <- capture.output(res <- daglad(cl,mediancenter=TRUE,param=c(d=9)))
	    return(res)
	  })
	  names(glad.profiles)<-sample.names
	}
	
	glad.profiles
	
}

#######################################################################################################################

#' rnb.plot.GLAD.profile
#' 
#' This functions returns a ReportPlot object with the GLAD profile figures added.
#' 
#' @param glad.profiles GLAD profiles as produced by \code{\link{getGLADprofiles}}.
#' @param label Name of the sample to be plotted.
#' @param sample.names Character vector containing all sample names.
#' @param numeric.names Something is replaced.
#' @param ... Further arguments passed to \code{createReportPlot}
#' @return A ReportPlot object with the GLAD profile plots added.
#' @author Pavlo Lutsik, with modifications by Michael Scherer
#' @noRd
rnb.plot.GLAD.profile<-function(glad.profile, label, sample.names = NA, numeric.names = FALSE, ...){
	
	cytoband_env<-new.env()
	data(cytoband, package="GLAD", envir = cytoband_env)
	
	plot.file<-createReportPlot(paste('GLADProfilePlot',  
					ifelse(!is.na(sample.names) && numeric.names, 
							match(label, sample.names), gsub("[ |_]", ".", label)) , sep="_"), ...)
		
	GLAD::plotProfile(glad.profile,cytoband=get("cytoband", envir = cytoband_env), Bkp=TRUE, Smoothing="Smoothing", main=paste(label),
	                  labels=F)
		
	off(plot.file)
	return(plot.file)
}

#######################################################################################################################

#' getCGCounts
#' 
#' This function returns the CpG-wise copy number alterations for each sample.
#' 
#' @param cnv.profiles CNV (GLAD) profiles as produced by \link{getGLADProfiles}.
#' @param rnb.set An object of type \code{\link{RnBeadRawSet}}.
#' @return CNV counts for each CpG in all samples
#' @author Pavlo Lutsik, with modifications by Michael Scherer
#' @noRd
getCGCounts<-function(cnv.profiles, rnb.set){
	target <- rnb.set@target
    assembly<-rnb.set@assembly
	cnv.reference.data<-rnb.get.cnv.annotations(target)
	chrom.data <- seqlengths(rnb.get.annotation(target, assembly))
	
	si<-cnv.reference.data$seq.info.hg19
	annot<-annotation(rnb.set)
	
	results<-sapply(cnv.profiles, function(profile){
		
		chr.basic<-GRanges(seqnames=ChrNumeric(names(chrom.data)), 
				IRanges(start=rep(1,length(chrom.data)), end=chrom.data), 
				strand=(rep("*", length(chrom.data))), regid=rep("native", length(chrom.data)))#,
				#seqinfo=si)		
				
		sample.cgh<-as.data.frame(profile)
		sample.bkps<-profile$BkpInfo
		
		if(!any(is.na(sample.bkps))){			
	  	sample.ranges<-GRanges(seqnames=sample.bkps$Chromosome, 
		  		IRanges(start=sample.bkps$PosBase,width=1), strand=rep("*", nrow(sample.bkps)),
			  	#seqinfo=si,
				  regid=sample.bkps[,1])
		}else{
		  sample.ranges <- GRanges()
		}
		sample.ranges<-c(chr.basic, sample.ranges)
		sample.ranges<-disjoin(sample.ranges)
		start(sample.ranges[which(width(sample.ranges)==1)+1L])<-start(sample.ranges[which(width(sample.ranges)==1)+1L])-1
		sample.ranges<-sample.ranges[width(sample.ranges)>1]
					
		sample.prof.ranges<-GRanges(seqnames=ChrNumeric(sample.cgh$Chromosome), 
				IRanges(start=sample.cgh$PosBase, width=1), strand="*")#, seqinfo=si)
		sample.olaps<-findOverlaps(sample.prof.ranges, sample.ranges)
		
		sample.ranges <- sample.ranges[unique(subjectHits(sample.olaps))]
		mean.zonegnl<-tapply(sample.cgh[queryHits(sample.olaps),"ZoneGNL"], subjectHits(sample.olaps), mean)
		Level<-round(mean.zonegnl) # expecting 2n, can be very wrong if the sample is on average kn, k>2
		sample.counts<-cbind(as.data.frame(sample.ranges),Level)	
		
		sample.counts$chromosome <- seqnames(sample.ranges)
		sample.counts$seqnames<-NULL
		sample.counts<-sample.counts[,c(6,1,2,4,5)]			
		
		cg.ranges<-GRanges(seqnames=ChrNumeric(annot$Chromosome), 
				IRanges(start=annot$Start, width=1), strand="*")#, seqinfo=si)
		cg.olaps<-findOverlaps(cg.ranges, sample.ranges)

		cnv.vector<-sample.counts[subjectHits(cg.olaps),"Level"]
	})
	
	colnames(results)<-names(cnv.profiles)
	row.names(results)<-row.names(annot)
	results
	
}

#######################################################################################################################
#' rnb.execute.cnv
#'
#' Copy number variation calling using GLAD
#'
#' @param object \code{\linkS4class{RnBeadRawSet}} object
#' @return a \code{list} with elements \code{glad.profiles} and \code{cg.counts}
#'
#' @author Pavlo Lutsik
#' @noRd

rnb.execute.cnv<-function(object){
	
    logger.warning("CNV calling currently not supported due to dependency issues")
	#glad.profiles <- getGLADProfiles(object)
	#cg.counts <- getCGCounts(glad.profiles, object)
    glad.profiles <- NULL
    cg.counts <- NULL
	
	list(cnv.profiles=glad.profiles, cg.counts=cg.counts)
	
}

#######################################################################################################################

#' rnb.section.cnv
#'
#' Adds CNV section to quality control report
#'
#' @param report analysis report
#' @param cnv.data a \code{list} output by \code{rnb.execute.cnv}
#'
#' @return The modified report object
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.section.cnv<-function(report, cnv.data){
	
	logger.start("Copy Number Variation Section")
	
	report <- rnb.add.section(report, "Copy number variation analysis", 
			"Visualization and analysis of copy number variations (CNVs) based on Infinium data.")
	
	txt <- 'CNV profile plots visualizes the results of CNV analysis using the <a href="http://bioconductor.org/packages/release/bioc/html/GLAD.html">GLAD</a> package.'
	if(rnb.getOption("qc.cnv.refbased")){
	  txt <- c(txt," \n A reference dataset was used to compute copy number gains and losses. The reference intensity values were obtained ",
	           "from a twin dataset (<a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100940'>GSE100940</a>), as the median intensity value for each position in the twins.")
	}else{
	  txt <- c(txt," \n The mean intensity value for each CpG in all samples was used as a reference to compute copy number alterations.")
	}
	rnb.add.paragraph(report, c("<b>CNV profiles:</b> ", txt))
	report<-add.profile.plots(report, cnv.data$cnv.profiles)
	logger.status("Added CNV profile plots")
	
	fname <- "cnv_per_cpg.csv"
	fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), fname)
	utils::write.csv(cnv.data$cg.counts, file = fname.full, na = "")
	
	txt<-c("Based on the segments with the identical copy number, one can calculate copy number for each CpG. ",
			"The table of computed copy numbers is available as a <a href=\"",
			rnb.get.directory(report, "data"), "/", fname, "\">comma-separated file</a> accompanying this report.",
			"The meanings of each entry are listed below.")
	rnb.add.paragraph(report, c("<b>Copy gains per CpG:</b> ", txt))
	exp.list <- list("<b>-1:</b> represents a copy number loss",
	                 "<b>1:</b> represents a copy number gain",
	                 "<b>2:</b> represents an amplicon",
	                 "<b>-10:</b> represents a deletion",
	                 "<b>0:</b> represents the normal state")
	rnb.add.list(report,exp.list)
	
	logger.status("Added CpG counts")
	
	logger.completed()
	return(report)
	
}
#######################################################################################################################

#' rnb.step.cnv
#'
#' Performs copy number calling from the Infinium intenstity data and adds the results to the report
#'
#' @param rnb.set	An object of type \code{\linkS4class{RnBeadRawSet}}
#' @param report  Report on quality control to contain the generated sections. This must be an object of type
#'                		\code{\linkS4class{Report}}.
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @export
rnb.step.cnv<-function(rnb.set, report){

    logger.warning("CNV calling currently not supported due to dependency issues")
  
	#rnb.require("GLAD")

	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	
	logger.start("Copy Number Variation")
	
	cnv.data<-rnb.execute.cnv(rnb.set)
	
	#report<-rnb.section.cnv(report, cnv.data)
	
	logger.completed()
	
	report
	
}

#######################################################################################################################
#######################################################################################################################

#' add.profile.plots
#' 
#' Produces GLAD profile plots for each sample and adds them to the report.
#' 
#' @param report The \code{\link{Report}} object, to which the plots should be added.
#' @param cnv.profiles CNV (GLAD) profiles as produced by \link{getGLADProfiles}.
#' @return The modified report object.
#' @author Pavlo Lutsik, with modifications by Pavlo Lutsik
#' @noRd
add.profile.plots<-function(report, cnv.profiles){
	
	descr <- c("Profiles visualize copy number variation across the genome of Infinium profiled samples. Each point represents a CpG on ",
"the chip. Yellow is the normal state, black is a deletion, blue an amplification, red a gain and green a loss. The CpGs are annotated ",
"to their position in the genome, with the centromers of the chromosome depicted in red. Red dashed lines represent breakpoints of regions ",
"with the identical amount of DNA.")
	
	ids<-names(cnv.profiles)
	
	if(parallel.isEnabled()){
	  cplots <- mclapply(ids, function(id) {
	    rnb.plot.GLAD.profile(glad.profile=cnv.profiles[[id]], label=id, sample.names=ids, 
	                          report=report, numeric.names=TRUE, width=8, height=7, low.png=100, high.png=300, create.pdf=F)
	  },mc.cores=parallel.getNumWorkers())
	}else{
	  cplots<-lapply(ids, function(id) {
				rnb.plot.GLAD.profile(glad.profile=cnv.profiles[[id]], label=id, sample.names=ids, 
						report=report, numeric.names=TRUE, width=8, height=7, low.png=100, high.png=300, create.pdf=F)
			})
	}
	
	names(cplots)<-1:length(ids)
	
	sn<-list("Sample labels" = ids)
	names(sn[[1]])<-1:length(ids)
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report
}
#######################################################################################################################
