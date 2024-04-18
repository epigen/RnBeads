########################################################################################################################
## qualityControl.R
## created: 2012-05-31
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
##
########################################################################################################################

#' rnb.step.quality
#'
#' Cteates quality control plots section of the quality control report.
#'
#' @param rnb.set Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report on quality control to contain the generated sections. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.step.quality<-function(rnb.set, report){

	if (!inherits(rnb.set, "RnBSet")){
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}

	if(inherits(rnb.set,"RnBeadSet")){
		covg.lists <- NULL
	}else{ # inherits(rnb.set,"RnBiseqSet")
		covg.lists<-rnb.execute.quality(rnb.set)
	}
	report<-rnb.section.quality(report, rnb.set, covg.lists=covg.lists)

	return(report)
}

########################################################################################################################

#' rnb.execute.quality
#'
#' Performs quality control calculations on the loaded DNA methylation data set.
#'
#' @param object            Methylation dataset as an object of class \code{\linkS4class{RnBeadSet}},
#'                          \code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBiseqSet}}.
#' @param type              \code{character} vector of length \code{1} giving the type of genomic regions for which the
#' 				            quality control information is summarized.
#' @param qc.coverage.plots Flag indicating if sequencing coverage information is summarized and returned. This
#'                          parameter is considered only when \code{object} is of type \code{\linkS4class{RnBiseqSet}}.
#' @param verbose	        Flag specifying whether diagnostic output should be written to the console or to the
#'                          RnBeads logger in case the latter is initialized.
#'
#' @details Currently, summarizing coverage for \code{\linkS4class{RnBiseqSet}} object is the only available function.
#'
#' @return \code{\linkS4class{RnBeadSet}} object with imputed quality control information
#'
#' @author Pavlo Lutsik
#' @export
rnb.execute.quality<-function(
		object,
		type="sites",
		qc.coverage.plots=rnb.getOption("qc.coverage.plots"),
		verbose=TRUE){

  if(verbose){
  	rnb.logger.start("Preparing Quality Control Information")
  }

  if(inherits(object, "RnBiseqSet")){
	  if(qc.coverage.plots){
		  covg.rnbs<-covg(object, type)
		  covg.rnbs[is.na(covg.rnbs)]<-0

		  if(type=="sites"){
			  type<-"CpG"
			  map<-sites(object)
		  }else{
			  map<-regions(object, type)
		  }

		  covg.lists<-lapply(samples(object), function(sample){

		  sample.id<-match(sample, samples(object))
		  assembly <- assembly(object)
		  covg.list<-lapply(names(rnb.get.chromosomes(assembly=assembly)), function(chr) {
					  chr.id<-match(chr, names(rnb.get.chromosomes(assembly=assembly)))
					  covg.weight<-rep(0, rnb.annotation.size(assembly=assembly)[chr.id])
					  covg.weight[map[map[,2]==chr.id,3]]<-covg.rnbs[map[,2]==chr.id,sample.id]
					  covg.weight
				  })
	  	  })
		  names(covg.lists)<-samples(object)
		  result<-covg.lists
  	  }else{
		  result<-NULL
	  }
  }else{
	  result<-object
  }
  if(verbose){
  	logger.completed()
  }
  return(result)
}

#######################################################################################################################

#' rnb.section.quality
#'
#' Adds quality control probe or coverage section to the quality control report.
#'
#' @param report             Analysis report to contain the newly generated section. This must be an object of type
#'                           \code{\linkS4class{Report}}.
#' @param rnb.set            Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param qc.boxplots        Flag indicating if quality control box plots are to be generated. This option has effect
#'                           only when \code{rnb.set} is a HumanMethylation450 dataset.
#' @param qc.barplots        Flag indicating if quality control bar plots are to be generated. This option has effect
#'                           only when \code{rnb.set} is a HumanMethylation450 dataset.
#' @param qc.negativeboxplot Flag indicating if plots showing the negative control probes for each sample are to be
#'                           generated. This option has effect only when \code{rnb.set} is a
#'                           HumanMethylation450 dataset.
#' @param covg.lists         ...
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.section.quality<-function(report, rnb.set, qc.boxplots=rnb.getOption("qc.boxplots"),
	qc.barplots=rnb.getOption("qc.barplots"), qc.negative.boxplot=rnb.getOption("qc.negative.boxplot"),
	qc.coverage.plots=rnb.getOption("qc.coverage.plots"), qc.coverage.histograms=rnb.getOption("qc.coverage.histograms"),
	qc.coverage.violins=rnb.getOption("qc.coverage.violins"),
	qc.coverage.threshold.plot=rnb.getOption("qc.coverage.threshold.plot"),covg.lists=NULL){

	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!parameter.is.flag(qc.boxplots)) {
		stop("invalid value for qc.boxplots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.barplots)) {
		stop("invalid value for qc.barplots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.negative.boxplot)) {
		stop("invalid value for qc.negative.boxplots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.coverage.plots)) {
		stop("invalid value for qc.coverage.plots; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.coverage.histograms)) {
		stop("invalid value for qc.coverage.histograms; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(qc.coverage.violins)) {
		stop("invalid value for qc.coverage.violins; expected TRUE or FALSE")
	}
	## TODO: Add error handling for qc.coverage.threshold.plot
	## TODO: Add error handling for covg.lists

	logger.start("Quality Control Section")
	report <- rnb.add.section(report, "Quality Control", NULL)
	txt <- NULL

	if(inherits(rnb.set, "RnBeadSet")){
		if(is.null(qc(rnb.set))){

			txt <- c("The supplied dataset contains no quality information, therefore, quality control graphics could ",
				"not be generated.")
			rnb.add.paragraph(report, txt)
			logger.info("No quality information present in the dataset")

		}else{

			txt <- "This section contains quality control plots and statistics for the methylation data."
			rnb.add.paragraph(report, txt)

			if(qc.boxplots) {
				txt <- c("Each box plot below shows the signal distribution of quality control probes across all ",
					"samples. The control box plots are separated by control types. Detailed description of the ",
					"control probes is given in the RnBeads vignette.")
				report <- rnb.add.section(report, "Quality Control Box Plots", txt, level = 2)
				report<-add.qc.boxplots(report, rnb.set)
				logger.status("Added quality control box plots")
			}

			if(qc.barplots) {
				txt <- c("The plots below visualize the exact signal levels at each quality control probe. Note that ",
					"the scale is not standardized. Background signal is usualy at the level of 1000 to 2000.")
				report <- rnb.add.section(report, "Quality Control Bar Plots", txt, level = 2)
				report<-add.qc.barplots(report, rnb.set, sample.batch.size=rnb.getOption("qc.sample.batch.size"))
				logger.status("Added quality control bar plots")
			}

			if(qc.negative.boxplot){
				txt <- c("Negative control box plots visualize background intensity distributions of all analyzed ",
					"samples. Samples with skewed distributions and high medians are likely to be of low quality and ",
					"should be discarded from the analysis.")
				report <- rnb.add.section(report, "Negative Control Box Plots", txt, level = 2)
				report<-add.negative.control.boxplot(report, rnb.set, sample.batch.size=rnb.getOption("qc.sample.batch.size"))
				logger.status("Added negative control boxplots")
			}
		}

	}else{ #inherits(rnb.set, "RnBiseqSet")


		if(qc.coverage.plots){
			txt <- c("The sequencing coverage plots visualized effective read coverage at all CpGs in the sequenced ",
					"genome. In case certain samples seem to have significantly decreased coverage, they should be excluded ",
					"from the analysis.")
			report <- rnb.add.section(report, "Sequencing Coverage Plots", txt, level = 2)
			report<-add.seq.coverage.plot(report, rnb.set, covg.lists)
			logger.status("Added sequencing coverage boxplots")
		}

		if(qc.coverage.histograms){
			txt <- c("The sequencing coverage histograms show distribution of coverage across all chromosomes. In case ",
			"certain samples seem to have significantly decreased coverage, they should be excluded from the analysis.")
			report <- rnb.add.section(report, "Sequencing Coverage Histograms", txt, level = 2)
			report <- add.seq.coverage.histograms(report, rnb.set)
			logger.status("Added sequencing coverage histograms")
		}
		report <- add.seq.coverage.num.sites.covg.tabs(report, rnb.set)
		logger.status("Added sample coverage section")

		if(qc.coverage.violins){
			txt <- c("The plots below show an alternative approach to visualizing the coverage distribution.")
			report <- rnb.add.section(report, "Sequencing Coverage Violin Plots", txt, level = 2)
			report<-add.seq.coverage.violins(report, rnb.set)
			logger.status("Added sequencing coverage violin plots")
		}

		rnb.cleanMem()

		if (length(qc.coverage.threshold.plot) != 0) {
			rplot <- rnb.plot.coverage.thresholds(rnb.set, qc.coverage.threshold.plot, fname = "coverage_interrogated",
				report = report)
			dframe.coverages <- attr(rplot, "data")
			fname <- "coverage_interrogated.csv"
			write.csv(dframe.coverages, file = file.path(rnb.get.directory(report, "data", TRUE), fname), row.names = FALSE)
			txt <- sprintf("%1.1f", range(dframe.coverages[dframe.coverages[, 2] == max(dframe.coverages[, 2]), 3]))
			txt <- c('In total, between ', txt[1], ' and ', txt[2], ' million sites are covered in all samples of the ',
				'dataset. The figure below shows the change in supports for different coverage thresholds. The exact ',
				'values are available in a dedicated <a href="', rnb.get.directory(report, 'data'), '/', fname,
				'">comma-separated file</a> accompanying this report.')
			report <- rnb.add.section(report, "Sequencing Coverage Thresholds", txt, level = 2)
			txt <- c("Line plot showing the number of CpG sites with a given support for different thresholds of ",
				"minimal coverage. The support of a CpG site is the minimal number of samples that interrogate it.")
			report <- rnb.add.figure(report, txt, rplot)
			rm(rplot, dframe.coverages, fname)
		}
	}

	if (is.null(txt)) {
		txt <- "No quality control plots are generated because all respective options are disabled."
		rnb.add.paragraph(report, txt)
	}

	logger.completed()
	return(report)
}

#######################################################################################################################

#' rnb.step.snp.probes
#'
#' Computes statistics on the SNP-based probes in the given dataset and creates a corresponding section in the report.
#'
#' @param object a \code{\linkS4class{RnBeadSet}} object
#' @param report Report on methylation profiles to contain the dimension reduction section. This must be an object of
#'               type \code{\linkS4class{Report}}.
#' @return the modified report object
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.step.snp.probes<-function(object, report){

	if(!inherits(object,"RnBSet")){
		stop("Supplied object is not of the class inheriting from RnBSet")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}

	logger.start("Visualizing SNP Probe Data")
	if (any(unlist(rnb.options("qc.snp.heatmap", "qc.snp.boxplot", "qc.snp.barplot", "qc.snp.purity")))) {
		logger.start("Mixups Visualization Section")
		report <- rnb.section.snp.probes(report,object)
		logger.completed()
	}
	logger.completed()

	return(report)
}

#######################################################################################################################

#' rnb.section.snp.probes
#'
#' Adds a section on the SNP-based probes to the given report.
#'
#' @param report Report to contain the new section. This must be an object of type \code{\linkS4class{Report}}.
#' @param object Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @return The (possibly modified) report.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.section.snp.probes <- function(report, object) {

	mm.snps <- tryCatch(rnb.get.snp.matrix(object), error = function(err) { NULL })
	section.title <- "Visualization of SNP Probe Data"
	if (is.null(mm.snps)) {
		txt <- c("Overview of SNP-based probes and sample comparison based on them cannot be performed ",
			"because the dataset does not contain such probes.")
		report <- rnb.add.section(report, section.title, txt)
		return(report)
	}
	if (ncol(mm.snps) == 1) {
		txt <- c("Overview of SNP-based probes and sample comparison based on them cannot be performed ",
			"because the dataset contains only one sample.")
		report <- rnb.add.section(report, section.title, txt)
		return(report)
	}
	
	txt <- c("Analysis of the values of the SNP-based probes can help identify sample mixups.")
	report <- rnb.add.section(report, section.title, txt)

	add.info <- function(stitle, ffunction, txt) {
		result <- rnb.add.section(report, stitle, txt, level = 2)
		result <- ffunction(result, mm.snps)
		rnb.status(c("Added", stitle))
		result
	}

	if (rnb.getOption("qc.snp.heatmap")) {
		txt <- "SNP heatmap enables the identification of sample mixups in particular for genetically matched designs."
		report <- add.info("SNP Heatmap", rnb.add.snp.heatmap, txt)
	}
	if (rnb.getOption("qc.snp.barplot")) {
		txt <- "SNP bar plots enable the identification of sample mixups in particular for genetically matched designs."
		report <- add.info("SNP Bar Plots", rnb.add.snp.barplot, txt)
	}
	if (rnb.getOption("qc.snp.boxplot")) {
		txt <- c("The SNP box plot is a visualization tool that can show sample mixups, especially in the case the ",
			"set of samples is genetically homogeneous.")
		report <- add.info("SNP Box Plot", rnb.add.snp.boxplot, txt)
	}
	if (rnb.getOption("qc.snp.distances")) {
		txt <- c("If we inspect the dataset in the space defined by the SNP probes only, samples appearing close to ",
			"each other are genetically similar.")
		report <- add.info("SNP-based Distances", rnb.add.snp.distances, txt)
	}
	if (rnb.getOption("qc.snp.purity")) {
		report <- rnb.add.snp.purity(report, mm.snps, object@pheno)
	}

	return(report)
}

########################################################################################################################
## UTILS
########################################################################################################################

add.qc.boxplots<-function(report, object){
  descr<-"Quality control box plots."

  if(object@target=="probesEPIC"){
	  ctypes<-rnb.infinium.control.targets(object@target)[c(14,4,3,15,1:2,12:13,6,11)]
  }else if(object@target=="probesEPICv2"){
	  ctypes<-rnb.infinium.control.targets(object@target)[c(14,4,3,15,1:2,12:13,6,11)]
  }else if(object@target=="probes450"){
  	ctypes<-rnb.infinium.control.targets(object@target)[c(13,4,3,14,1:2,11:12,6)]
  }else if(object@target=="probes27"){
	ctypes<-rnb.infinium.control.targets(object@target)[c(10,3,2,11,1,9,6)]
  }else if(object@target=="probesMMBC"){
      ctypes<-rnb.infinium.control.targets(object@target)[c(14,4,3,15,1:2,12:13,6,11)]
  }

  cplots<-lapply(ctypes, rnb.plot.control.boxplot, rnb.set=object, report=report, writeToFile=TRUE, numeric.names=TRUE, width=8, height=6, low.png=100, high.png=300)
  names(cplots)<-1:length(ctypes)

  sn<-list("Control probe type"=ctypes)
  names(sn[[1]])<-1:length(ctypes)

  report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
  report

}

#######################################################################################################################

add.qc.barplots<-function(report, object, sample.batch.size=50){
  descr="Quality control bar plots."

  if(object@target=="probesEPIC"){
	  cmd <- rnb.get.annotation("controlsEPIC", assembly = rnb.getOption("assembly"))
	  ctypes<-unique(cmd$Target)[unique(cmd$Target) %in% rnb.infinium.control.targets("probesEPIC")[c(14,4,3,15,1:2,12:13,6,11)]]
  }else if(object@target=="probesEPICv2"){
	  cmd <- rnb.get.annotation("controlsEPICv2", assembly = "hg38")
	  ctypes<-unique(cmd$Target)[unique(cmd$Target) %in% rnb.infinium.control.targets("probesEPICv2")[c(14,4,3,15,1:2,12:13,6,11)]]
  }else if(object@target=="probes450"){
  	cmd <- rnb.get.annotation("controls450", assembly = rnb.getOption("assembly"))
  	ctypes<-unique(cmd$Target)[unique(cmd$Target) %in% rnb.infinium.control.targets("probes450")[c(13,4,14,3,1:2,11:12,6)]]
  }else if(object@target=="probes27"){
	cmd <- rnb.get.annotation("controls27")
	ctypes<-unique(cmd$Type)[unique(cmd$Type) %in% rnb.infinium.control.targets("probes27")[c(10,3,2,11,1,9,6)]]
  }else if(object@target=="probesMMBC"){
    cmd <- rnb.get.annotation("controlsMMBC", assembly="mm10")
    ctypes<-unique(cmd$Target)[unique(cmd$Target) %in% rnb.infinium.control.targets("probesMMBC")[c(14,4,3,15,1:2,12:13,6,11)]]
  }
  nsamp<-length(samples(object))

  plot.names<-NULL

  if(nsamp %% sample.batch.size==1){
	  sample.batch.size<-sample.batch.size-5
  }
  portion.starts<-0:(nsamp %/% sample.batch.size)*sample.batch.size+1
  portion.ends<-portion.starts+sample.batch.size-1
  portion.ends[length(portion.ends)]<-nsamp
  portions<-paste(portion.starts, portion.ends, sep="-")

  plots<-lapply(1:length(portions),function(portion.id){

	  cplots<-lapply(ctypes, function(type){

		if(object@target=="probes450" || object@target=="probesEPIC" || object@target=="probesMMBC" || object@target=="probesEPICv2"){
			cmdt <- cmd[cmd[["Target"]] == type, ]
			pn<-paste(type, 1:(dim(cmdt)[1]),  sep=".")
		}else if(object@target=="probes27"){
			cmdt <- cmd[cmd[["Type"]] == type, ]
			pn<-as.character(cmdt$Name)
		}

		if(portion.id==1) plot.names<<-c(plot.names, pn)
	    plots<-lapply(pn, rnb.plot.control.barplot, rnb.set=object,
				sample.subset=portion.starts[portion.id]:portion.ends[portion.id],
				report=report, writeToFile=TRUE, numeric.names=TRUE, width=8, height=6, low.png=100, high.png=300, verbose=TRUE,
				name.prefix=portions[portion.id])

		if(object@target=="probes450" || object@target=="probesEPIC" || object@target=="probesMMBC" || object@target=="probesEPICv2"){
			names(plots)<-paste(type, 1:(dim(cmdt)[1]))
		}else if(object@target=="probes27"){
			names(plots)<-as.character(cmdt$Name)
		}
	    plots

	  })

	  names(cplots)<-NULL

	  cplots<-unlist(cplots)
	  names(cplots)<-1:length(plot.names)
	  cplots
  })

  plots<-unlist(plots)

  sn<-list("Samples #: " = portions, "Control probe ID" = plot.names)


  names(sn[[1]])<-portions
  if(object@target=="probes450" || object@target=="probesEPIC" || object@target == "probesMMBC" || object@target == "probesEPICv2"){
  	names(sn[[2]])<-1:length(plot.names)
  }else if(object@target=="probes27"){
	names(sn[[2]])<-match(plot.names,cmd$Name[cmd$Type %in% rnb.infinium.control.targets("probes27")[c(10,3,2,11,1,9,6)]])
  }

  report<-rnb.add.figure(report, description=descr, report.plots=plots, setting.names=sn)

}

#######################################################################################################################

add.negative.control.boxplot<-function(report, object, sample.batch.size=50){

	descr<-"Box plots of the negative control probes."

#   cplots<-lapply(ctypes, control.boxplot, rnb.set=object, report=report, writeToFile=TRUE)
#   names(cplots)<-gsub(" ", ".", ctypes)
#
#   sn<-list(tolower(ctypes))
#   names(sn[[1]])<-gsub(" ", ".", ctypes)
	nsamp<-length(samples(object))
	if(nsamp %% sample.batch.size==1){
		sample.batch.size<-sample.batch.size-5
	}
	portion.starts<-0:(nsamp %/% sample.batch.size)*sample.batch.size+1
	portion.ends<-portion.starts+sample.batch.size-1
	portion.ends[length(portion.ends)]<-nsamp
	portions<-paste(portion.starts, portion.ends, sep="-")

	cplots<-lapply(1:length(portions), function(portion.id){
				rnb.plot.negative.boxplot(object, sample.subset=portion.starts[portion.id]:portion.ends[portion.id],
						name.prefix=portions[portion.id], writeToFile=TRUE, report=report, width=10, height=6, low.png=75, high.png=100)
	})

	sn<-list("Samples #: " = portions)

	names(sn[[1]])<-portions

	report<-rnb.add.figure(report, description=descr, report.plots=cplots, sn)
	report

}

#######################################################################################################################

rnb.add.snp.boxplot<-function(report, object){

  txt <- "Box plot of the SNP probes."
  rplot <- rnb.plot.snp.boxplot(object, writeToFile=TRUE, report=report, width=9, height=6, low.png=100, high.png=300)
  rnb.add.figure(report, txt, rplot)
}

#######################################################################################################################

rnb.add.snp.barplot<-function(report, object){

	txt <- "Bar plot of all observed beta values for a SNP probe."
	rplots <- list()
	for (probe.id in rownames(object)) {
		rplot <- rnb.plot.snp.barplot(object, probe.id, TRUE, TRUE, report = report, width=9, height=6, high.png=300)
		rplots <- c(rplots, list(rplot))
	}
	setting.names <- list("Probe" = rownames(object))
	names(setting.names[[1]]) <- 1:length(setting.names[[1]])
	rnb.add.figure(report, txt, rplots, setting.names)
}

#######################################################################################################################

rnb.add.snp.heatmap<-function(report, object){
	txt <- paste0("Heatmap of the SNP probes. Euclidean distance and complete linkage are used for constructing the ",
		"dendrograms. Samples with the same genetic background are expected to cluster together.")
	rplot <- rnb.plot.snp.heatmap(object, writeToFile=TRUE, report=report, width=8, height=9, low.png=100, high.png=300)
	rnb.add.figure(report, txt, rplot)
}

#######################################################################################################################

#' rnb.add.snp.distances
#'
#' Adds a section about sample distances based on beta values of SNP probes.
#'
#' @param report Report to contain the section on SNP probe distances.
#' @param object Non-empty \code{matrix} containing the computed beta values on the SNP probes.
#' @return The (possibly modified) report.
#' @author Yassen Assenov
#' @noRd
rnb.add.snp.distances <- function(report, object) {

	## Extract the matrix of beta values on the SNP probes
	snp.betas <- t(object[!apply(is.na(object), 1, any), , drop = FALSE])
	if (ncol(snp.betas) <= 1) {
		snp.betas <- NULL
	}
	if (is.null(snp.betas)) {
		txt <- c("Distances based on SNP probes could not be calculated either because no such probes are found ",
			"in the dataset, or because almost all of them contain missing values.")
		rnb.add.paragraph(report, txt)
		return(report)
	}

	## Calculate Manhattan distances between samples based on the SNP probe intensities
	snp.distances <- as.matrix(stats::dist(snp.betas, method = "manhattan") / ncol(snp.betas))
	rnb.status("Calculated Manhattan distances between samples based on SNP probes")

	report.plots <- list()
	setting.names <- list()
	if (nrow(snp.betas) <= 24) {

		## Create a diagonal heatmap of distances
		i.width <- 4 + nrow(snp.betas) * 0.3
		i.height <- 2.2 + nrow(snp.betas) * 0.3
		txt <- paste("The figure below shows the relative distances between all pairs of samples based on the",
			"&beta; values of the considered SNP probes. The distance metric used is average absolute difference,",
			"which can be considered a scaled version of Manhattan distance.")

		tbl <- symmetric.melt(snp.distances)
		colnames(tbl)[3] <- "distance"
		colors.g <- rnb.getOption("colors.gradient")
		pp <- ggplot(tbl, aes_string(x = "x", y = "y", fill = "distance")) + labs(x = NULL, y = NULL) +
			coord_fixed(ratio = 1) + geom_tile(color = "white") +
			scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
			scale_fill_gradient(na.value = "white", low = colors.g[1], high = colors.g[2]) +
			theme(legend.justification = c(0, 1), legend.position = c(1, 1)) +
			theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
			theme(panel.grid.major = element_blank(), panel.background = element_blank()) +
			theme(panel.border = element_blank(), plot.margin = unit(c(0.1, 1.1, 0.1, 0.1), "in"))
		## Fix the areas for x and y axis labels
#		pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
#		pp$widths[[3]] <- unit(2, "in")
#		pp$heights[[length(pp$heights) - 2L]] <- unit(2, "in")
		rplot <- createReportPlot("snp_low_dimensional", report, width = i.width, height = i.height)
#		grid.newpage()
#		grid.draw(pp)
		print(pp)
		report.plots <- c(report.plots, off(rplot))
		txt <- c(txt, paste("Distances between pairs of samples based on", ncol(snp.betas), "SNP probes."))
		rm(colors.g)

	} else {

		## Create scatter plots
		i.width <- 6.2 + (nrow(snp.betas) >= 60)
		i.height <- i.width
		if (ncol(snp.betas) == 2) {
			tbl <- snp.betas
			alabels <- colnames(snp.betas)
			txt <- "The figure below shows the samples in the space defined by the two SNP probes."
		} else {
			tbl <- prcomp(snp.betas, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
			alabels <- paste("Principal component", 1:ncol(tbl))
			txt <- paste("The figure below shows the samples in the first two principal components of the space",
				"defined by the", ncol(snp.betas), "SNP probes.")
		}

		tbl <- data.frame(x = tbl[, 1], y = tbl[, 2], label = rownames(tbl))
		plot.types <- c("point" = "points", "label" = "identifiers")
		for (ptype in names(plot.types)) {
			pp <- ggplot(tbl, aes_string(x = "x", y = "y", label = "label")) + coord_fixed(ratio = 1)
			if (ncol(snp.betas) == 2) {
				pp <- pp + scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, length.out = 11), expand = c(0, 0)) +
					scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, length.out = 11), expand = c(0, 0))
			}
			if (ptype == "point") {
				pp <- pp + geom_point(size = 3)
			} else {
				pp <- pp + geom_text(size = 2, angle = 45)
			}
			pp <- pp + labs(x = alabels[1], y = alabels[2]) + theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
			fname <- paste0("snp_low_dimensional_", ptype)
			rplot <- createReportPlot(fname, report, width = i.width, height = i.height)
			print(pp)
			report.plots <- c(report.plots, off(rplot))
		}
		txt <- c(txt, "Scatter plot showing the samples in a space defined by the signal of their SNP probes.")
		setting.names <- list("Display samples as" = plot.types)
		rm(alabels, plot.types, ptype, fname)
	}
	rnb.add.paragraph(report, txt[1])
	report <- rnb.add.figure(report, txt[2], report.plots, setting.names)
	rm(snp.betas, report.plots, setting.names, txt, i.width, i.height, tbl, pp, rplot)

	## Export the distances to file
	fname <- "snp_distances.csv"
	write.csv(snp.distances, file = file.path(rnb.get.directory(report, "data", TRUE), fname))
	txt <- c('The full table of all pairwise distances is stored in a dedicated <a href="',
		rnb.get.directory(report, "data"), '/', fname, '">comma-separated file</a> accompanying this report.')
	rnb.add.paragraph(report, txt)

	report
}

#######################################################################################################################

#' rnb.add.snp.purity
#'
#' Adds a section about genetic purity estimated using the beta values of SNP probes.
#'
#' @param report       Report to contain the section on genetic purity.
#' @param mm.snp       Non-empty \code{matrix} containing the computed beta values on the SNP probes.
#' @param s.annotation Optionally, a sample annotation table to be scanned for information on Sentrix IDs.
#' @return The (possibly modified) report.
#' @author Yassen Assenov
#' @noRd
rnb.add.snp.purity <- function(report, mm.snp, s.annotation = data.frame()) {
	tbl <- pmin(mm.snp, abs(mm.snp - 0.5), abs(mm.snp - 1))
	tbl <- data.frame(
		"Genetic Noise" = unname(colMeans(tbl, na.rm = TRUE)),
		"ID" = factor(colnames(mm.snp), levels = colnames(mm.snp)),
		check.names = FALSE, stringsAsFactors = FALSE)
	ids <- s.annotation
	ids <- ids[, intersect(c("Sentrix ID", "Sentrix_ID"), colnames(ids)), drop = FALSE]
	if (ncol(ids) != 0) {
		ids <- as.character(ids[, 1])
		id.levels <- tryCatch(unique(sort(as.double(ids))), warning = function(wr) { unique(ids) })
		ids <- factor(ids, levels = id.levels)
		if (length(id.levels) > 1 && anyDuplicated(na.omit(ids)) != 0) {
			tbl[["Sentrix ID"]] <- ids
		}
		rm(ids, id.levels)
	}
	tbl <- tbl[, intersect(c("ID", "Sentrix ID", "Genetic Noise"), colnames(tbl)), drop = FALSE]
	tbl <- tbl[order(tbl[, "Genetic Noise"]), , drop = FALSE]

	## Export the table to a file
	fname <- "genetic_noise.csv"
	write.csv(tbl, file.path(rnb.get.directory(report, "data", TRUE), fname), row.names = FALSE)

	## Plot the data
	tbl[, "Genetic Noise"] <- tbl[, "Genetic Noise"] * 100
	xstep <- 2
	xmax <- ceiling(max(tbl[, "Genetic Noise"]) / xstep) * xstep
	i.width <- 6.2
	pp <- ggplot2::ggplot(tbl) + ggplot2::aes_string(x = '`ID`', y = '`Genetic Noise`') +
		ggplot2::labs(x = NULL, y = "Genetic noise (%)")
	if ("Sentrix ID" %in% colnames(tbl)) {
		i.width <- i.width + 1.4
		pp <- pp + ggplot2::aes_string(fill = '`Sentrix ID`') +
			ggplot2::theme(legend.position = c(1, 0.5), legend.justification = c(0, 0.5)) +
			ggplot2::theme(plot.margin = grid::unit(c(0, 1.4, 0, 0) + 0.1, "in"))
	}
	pp <- pp + ggplot2::geom_bar(stat = "identity") + ggplot2::coord_flip() +
		ggplot2::scale_y_continuous(breaks = seq(0, xmax, by = xstep), limits = c(0, xmax), expand = c(0, 0)) +
		ggplot2::theme(axis.ticks.y = ggplot2::element_blank())
	rplot <- createReportPlot("genetic_noise", report, width = i.width, height = nrow(tbl) * 0.12 + 1.2)
	print(pp)
	rplot <- off(rplot)
	
	## Add a section to the report
	txt <- c("The relative purity of a sample can be estimated by calculating the mean absolute distance from the ",
		"observed &beta; values of its SNP probes to an idealized profile, consisting of the values 0, 0.5 and 1 ",
		"only. This score is referred to as <em>genetic noise</em> in this report. Samples with large genetic noise ",
		"might contain subpopulations with genomic aberrations (deletions and amplifications), or be contaminated ",
		"with DNA from multiple individuals. It is important to note that the genetic noise could be strongly ",
		"influenced by batch effects and should therefore be considered only in the context of samples with similar ",
		"profiles of intensity values.")
	report <- rnb.add.section(report, "Genetic Purity", txt, 2)
	txt <- c('The figure below lists the genetic noise of all samples in the dataset. The full table with these ',
		'values is available in a dedicated <a href="', rnb.get.directory(report, 'data'), '/', fname,
		'">comma-separated value</a> file accompanying this report.')
	rnb.add.paragraph(report, txt)
	txt <- c("Bar plot showing the genetic noise of each sample. Color coding, if present, denotes slide number.")
	report <- rnb.add.figure(report, txt, rplot)
}

#######################################################################################################################

add.seq.coverage.plot<-function(report, object, covg.lists=NULL){

	txt <- "Effective read coverage of a single sample over all chromosomes of the genome."

	ids <- samples(object)
	names(ids) <- 1:length(ids)

	rplots<-lapply(ids, function(id) rnb.plot.biseq.coverage(rnbs.set=object, sample=id, report=report,
		writeToFile=TRUE, numeric.names=TRUE, create.pdf=FALSE, width=8, height=16, low.png=100, high.png=200,
		covg.lists=covg.lists[[id]]))

	rnb.add.figure(report, txt, rplots, list("Sample" = ids))
}


#######################################################################################################################

add.seq.coverage.histograms <- function(report, object){

	descr<-"Sequencing coverage histogram visualize the bulk distribution of read coverage for each sample."

	ids <- samples(object)

	cplots<-lapply(ids, rnb.plot.biseq.coverage.hist, rnbs.set=object, report=report, writeToFile=TRUE, numeric.names=TRUE, covg.max.percentile=0.99, width=8, height=7, low.png=100, high.png=300)

	names(cplots)<-1:length(ids)

	sn<-list("Sample labels" = ids)
	names(sn[[1]])<-1:length(ids)

	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report

}

#######################################################################################################################

add.seq.coverage.num.sites.covg.tabs <- function(report, object){
	selCols <- c("sampleName", "sites_num", "sites_covgMean", "sites_covgMedian", paste0("sites_covgPerc",c("25","75")), paste0("sites_numCovg",c("5","10","30","60")))
	summary.tab <- rnb.sample.summary.table(object)[,selCols]

	fname <- paste("sample_covg_summary.csv",sep="")
	fname <- rnb.write.table(summary.tab, fname,fpath=rnb.get.directory(report, "data", absolute = TRUE), format="csv", gz=FALSE, row.names = FALSE, quote=FALSE)
	sectionText <- paste("This section contains summary metrics on the number of sites and coverage can be found in the table below.",
		" The summary table below is also available as ",
			"<a href=\"", rnb.get.directory(report, "data"), "/", fname,"\">csv file</a>.",sep="")
	column.description <- c(
		"<ul>",
			"<li>sampleName: Name of the sample</li>",
			"<li>sites_num: Number of sites covered in the sample</li>",
			"<li>sites_covgMean: Mean coverage of sites in the sample</li>",
			"<li>sites_covgMedian: Median coverage of sites in the sample</li>",
			"<li>sites_covgPerc25,75: 25 and 75 percentile of coverage of sites in the sample</li>",
			"<li>sites_numCovg5,10,30,60: Number of sites with coverage greater or equal to 5,10,30,60</li>",
		"</ul>"
	)
	sectionText <- paste(sectionText, paste(column.description, collapse="\n"))

	report <- rnb.add.section(report, "Sample coverage summary", sectionText, level=2)
	rnb.add.table(report, summary.tab)

	pp <- rnb.plot.num.sites.covg(object)
	fname <- paste("num_site_vs_coverage",sep="_")
	overviewPlot <- createReportGgPlot(pp, fname, report, width=12, height=7, create.pdf=TRUE, high.png=as.integer(300))
	overviewPlot <- off(overviewPlot,handle.errors=TRUE)

	desc <- "Covered sites and median coverages for each sample. Vertical bars depict the inter-quartile range."
	report <- rnb.add.figure(report, desc, overviewPlot)
	report
}

#######################################################################################################################

add.seq.coverage.violins<-function(report, object, sample.chunk.size=20){

	descr<-paste("Sequencing coverage histogram visualized as violin plots. Distributions are based on",
				 nsites(object),"methylation sites.")

	ids <- samples(object)
	n <- length(ids)
	num.intervals <- ceiling(n/sample.chunk.size)
	if (num.intervals>1){
		ssets <- cut(1:n,ceiling(n/sample.chunk.size),labels=FALSE)
	} else {
		ssets <- rep(1,n)
	}
	ssets.n <- length(unique(ssets))
	max.covg.2plot<-quantile(covg(object),probs=c(0.99), na.rm=TRUE)[1]

	cplots <- list()
	for (i in 1:ssets.n){
		ssamples <- ids[ssets==i]
		rep.plot <- rnb.plot.biseq.coverage.violin(rnbs.set=object,samples=ssamples,fname=paste("coverageViolin_chunk",i,sep=""),covg.range=c(0,max.covg.2plot),report=report, width=8, height=7, low.png=100, high.png=300)
		cplots<-c(cplots,list(rep.plot))
	}

	names(cplots)<-1:ssets.n

	sn<-list("Sample chunk" = 1:ssets.n)
	names(sn[[1]])<-paste("chunk",1:ssets.n,sep="")

	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report
}
