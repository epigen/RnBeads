########################################################################################################################
## controlPlots.R
## created: 2012-04-04
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Quality control plots for HumanMethylation450 experiments.
########################################################################################################################

## G L O B A L S #######################################################################################################

.control.plot.theme.samples<-ggplot2::theme(
					plot.margin=grid::unit(c(0.1,1,0.1,0.1), "cm"),
					axis.title.x = ggplot2::element_blank(),
					axis.text.x=ggplot2::element_text(size=7, angle=-45, hjust=0))

HM450.SAMPLE.INDEPENDENT<-c("STAINING", "EXTENSION", "HYBRIDIZATION", "TARGET REMOVAL")
HM450.SAMPLE.DEPENDENT<-c("BISULFITE CONVERSION I","BISULFITE CONVERSION II",
		"SPECIFICITY I","SPECIFICITY II", "NON-POLYMORPHIC" )

## F U N C T I O N S ###################################################################################################

#' rnb.plot.control.boxplot
#'
#' Box plots of various control probes
#'
#' @param rnb.set \code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object with valid quality
#' 				  control information.
#' @param type    type of the control probe; must be one of the \code{"BISULFITE CONVERSION I"},
#'                \code{"BISULFITE CONVERSION II"}, \code{"EXTENSION"}, \code{"HYBRIDIZATION"}, \code{"NEGATIVE"},
#'                \code{"NON-POLYMORPHIC"}, \code{"NORM_A"}, \code{"NORM_C"}, \code{"NORM_G"}, \code{"NORM_T"},
#'                \code{"SPECIFICITY I"}, \code{"SPECIFICITY II"}, \code{"STAINING"}, \code{"TARGET REMOVAL"}.
#' @param writeToFile flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param ... other arguments to \code{\link{createReportPlot}}
#'
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class
#' 			\code{\link{ggplot}} otherwise.
#'
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.control.boxplot(rnb.set.example)
#' }
#'
#' @author Pavlo Lutsik
#' @export
rnb.plot.control.boxplot <- function(
		rnb.set,
		type = rnb.infinium.control.targets(rnb.set@target)[1],
		writeToFile = FALSE,
		numeric.names = FALSE,
		...) {

	if(!type %in% rnb.infinium.control.targets(rnb.set@target)) {
		stop(paste("Unrecognized control type:", type))
	}

	## Extract intensities of the control probes
	if(rnb.set@target=="probesEPIC"){
		meta <- rnb.get.annotation("controlsEPIC", assembly = "hg19")
	}else if(rnb.set@target=="probesEPICv2"){
		meta <- rnb.get.annotation("controlsEPICv2", assembly = "hg38")
	}else if(rnb.set@target=="probes450"){
		meta <- rnb.get.annotation("controls450", assembly = "hg19")
	}else if(rnb.set@target=="probes27"){
		meta <- rnb.get.annotation("controls27")
	}else if(rnb.set@target=="probesMMBC"){
        meta <- rnb.get.annotation("controlsMMBC", assembly="mm10")
    }

	if(rnb.set@target=="probesEPIC"){
		types<-rnb.infinium.control.targets(rnb.set@target)[c(14,4,3,15,1:2,12:13,6,11)]
	}else if(rnb.set@target=="probesEPICv2"){ 
		types<-rnb.infinium.control.targets(rnb.set@target)[c(14,4,3,15,1:2,12:13,6,11)]
	}else if(rnb.set@target=="probes450"){
		types<-rnb.infinium.control.targets(rnb.set@target)[c(13,4,3,14,1:2,11:12,6)]
	}else if(rnb.set@target=="probes27"){
		types<-rnb.infinium.control.targets(rnb.set@target)[c(10,3,2,11,1,9,6)]
	}else if(rnb.set@target=="probesMMBC"){
        types<-rnb.infinium.control.targets(rnb.set@target)[c(14,4,3,15,1:2,12:13,6,11)]
    }

	if(!type %in% types){
		warning("Unoptimized probe type, plotting performance may be decreased")
	}

	if(rnb.set@target=="probes450" || rnb.set@target=="probesEPIC" || rnb.set@target=="probesMMBC" || rnb.set@target=="probesEPICv2"){
		rownames(meta)<-meta[["ID"]]
		### TODO: Remove the following passage
		### for testing purposes only!
		if(rnb.set@target=="probesEPIC" || rnb.set@target=="probesEPICv2"){ 
			meta<-rnb.update.controlsEPIC.enrich(meta) ## TODO: Not validated for EPICv2
		}
		meta <- meta[type == meta[["Target"]], ]
		ids<-as.character(meta[["ID"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[ids,]
	}else if(rnb.set@target=="probes27"){
		meta <- meta[type == meta[["Type"]], ]
		ids<-as.character(meta[["Address"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[meta$Address %in% ids,]
	}

	if(nrow(meta)==0){
		plot.obj<-rnb.message.plot(sprintf("Box plot for control type %s not available", type))
		if(writeToFile){
			plot.file<-createReportGgPlot(plot.obj, paste("ControlBoxPlot", ifelse(numeric.names, match(type, types), gsub(" ", ".", type)) , sep="_"), ...)
			
			print(plot.obj)
			off(plot.file)
			return(plot.file)
		}
		return(plot.obj)
	}

	intensities <- lapply(c("green" = "Cy3", "red" = "Cy5"), function(channel) {
		result <- qc(rnb.set)[[channel]][ids, ]
		if (is.vector(result)) {
			return(as.matrix(result))
		}
		return(t(result))
	})

	scales<-lapply(qc(rnb.set), get.unified.scale)

	if(rnb.set@target=="probes450" || rnb.set@target=="probesEPIC" || rnb.set@target=="probesMMBC" || rnb.set@target=="probesEPICv2"){
		## Shorten the words describing probe's expected intensity
		INTENSITIES <- c("Background" = "Bgnd", "High" = "High", "Low" = "Low", "Medium" = "Med")
		levels(meta[, "Expected Intensity"]) <- INTENSITIES[levels(meta[, "Expected Intensity"])]
		meta[, "Expected Intensity"] <- as.character(meta[, "Expected Intensity"])

		## Define bar labels in the plot
		plot.names<-sapply(1:length(ids), sprintf, fmt="Probe%2g")

		exp.intensity <- lapply(c("green" = "Green", "red" = "Red"), function(channel) {
			ifelse(meta[[paste("Evaluate", channel)]] == "+", meta[["Expected Intensity"]], "Bgnd")
		})
		plot.names <- lapply(exp.intensity, function(exps) { paste(plot.names, exps, sep = ":\n") })

		plot.data <- lapply(names(intensities), function(channel) {
			channel.data <- intensities[[channel]]
			data.frame(
				"Intensity" = as.numeric(channel.data),
				"Probe" = as.factor(sapply(plot.names[[channel]], rep, times=nrow(channel.data))))
		})
	}else if(rnb.set@target=="probes27"){
		plot.names<-list()
		plot.names$green <- gsub("conversion |Extension |mismatch ", "", meta$Name)
		plot.names$red <- gsub("conversion |Extension |mismatch ", "", meta$Name)
		plot.data <- lapply(names(intensities), function(channel) {
					channel.data <- intensities[[channel]]
					data.frame(
							"Intensity" = as.numeric(channel.data),
							"Probe" = as.factor(sapply(plot.names[[channel]], rep, times=nrow(channel.data))))
				})
	}

	if(writeToFile){
		plot.file<-createReportPlot(paste("ControlBoxPlot", ifelse(numeric.names, match(type, types), gsub(" ", ".", type)) , sep="_"), ...)
		grid.newpage()
	}

	## Define viewports and assign it to grid layout
	#grid.newpage()
	#pushViewport(viewport(layout = grid.layout(2,1)))

	plots<-list()

	for (i in 1:length(intensities)) {
		channel <- names(intensities)[i]
#		print(qplot(Probe, Intensity, data=na.omit(plot.data[[i]]), geom = "boxplot",
		plots[[i]]<-qplot(Probe, Intensity, data=na.omit(plot.data[[i]]), geom = "boxplot",
				main=paste(type, paste(channel, "channel"), sep=": "),
				ylab="Intensity",fill=I(channel)) + scale_y_continuous(limits=scales[[i]])
		#,vp=viewport(layout.pos.row=i, layout.pos.col=1))
	}

	## FIXME: Rewrite the plotting function without using the function arrangeGrob in gridExtra
	grb<-do.call(arrangeGrob, plots)

	grid.draw(grb)
	if(writeToFile) {
		off(plot.file)
		return(plot.file)
	}else{
		return(invisible(grb))
	}
}

#######################################################################################################################

#' rnb.plot.negative.boxplot
#'
#' Box plots of negative control probes
#'
#' @param rnb.set 		\code{\linkS4class{RnBeadSet}} object with valid quality control information
#' @param sample.subset an integer vector specifying the subset of samples for which the plotting should be performed
#' @param writeToFile 	flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param name.prefix	in case \code{writeToFile} is \code{TRUE}, a \code{character} singleton specifying a prefix to the variable part of the image file names
#' @param ... 			other arguments to \code{\link{createReportPlot}}
#'
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class
#' 			\code{\link{ggplot}} otherwise.
#'
#' @author Pavlo Lutsik
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.negative.boxplot(rnb.set.example)
#' }
rnb.plot.negative.boxplot<- function(
		rnb.set,
		sample.subset=1:length(samples(rnb.set)),
		writeToFile = FALSE,
		name.prefix=NULL,
		...) {

	
	if(rnb.set@target=="probesEPIC" || rnb.set@target=="probesEPICv2" || rnb.set@target=="probesMMBC" ){
		genome.assembly<-rnb.getOption("assembly")
		meta <- rnb.get.annotation(gsub("probes", "controls", rnb.set@target), assembly=c("probesEPIC"=genome.assembly, "probesMMBC"="mm10", "probesEPICv2"="hg38")[rnb.set@target])
		## Extract intensities of the control probes
		### TODO: Remove the following passage
		### for testing purposes only!
		if(rnb.set@target=="probesEPIC" || rnb.set@target=="probesEPICv2"){
			meta<-rnb.update.controlsEPIC.enrich(meta) ## TODO: Not validated for EPICv2 (duplicate)
		}
		meta <- meta["NEGATIVE" == meta[["Target"]], ]
		ids<-as.character(meta[["ID"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[ids,]
	}else if(rnb.set@target=="probes450"){
		meta <- rnb.get.annotation("controls450")
		meta <- meta["NEGATIVE" == meta[["Target"]], ]
		ids<-as.character(meta[["ID"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[ids,]
	}else if(rnb.set@target=="probes27"){
		meta <- rnb.get.annotation("controls27")
		meta <- meta["Negative" == meta[["Type"]], ]
		ids<-as.character(meta[["Address"]])
		ids<-intersect(rownames(qc(rnb.set)[[1]]), ids)
		meta<-meta[meta$Address %in% ids,]
	}
	if(nrow(meta)==0){
		if(writeToFile) plot.file<-createReportPlot(paste("ControlBoxPlot", ifelse(numeric.names, match(type, types), gsub(" ", ".", type)) , sep="_"), ...)
		print(rnb.message.plot("Box plot for negative control probes is not available"))
		if(writeToFile){
			off(plot.file)
			return(plot.file)
		}
		return(invisible())
	}

	intensities <- lapply(c("green" = "Cy3", "red" = "Cy5"), function(channel) {
				result <- qc(rnb.set)[[channel]][ids,sample.subset]
				if (is.vector(result)) {
					return(as.matrix(result))
				}
				return(result)
			})

	## Shorten the words describing probe's expected intensity
#	INTENSITIES <- c("Background" = "Bgnd", "High" = "High", "Low" = "Low", "Medium" = "Med")
#	levels(meta[, "Expected Intensity"]) <- INTENSITIES[levels(meta[, "Expected Intensity"])]
#	meta[, "Expected Intensity"] <- as.character(meta[, "Expected Intensity"])

	## Define bar labels in the plot
#	plot.names<-sapply(1:length(ids), sprintf, fmt="Probe%2g")
	sample.ids<-samples(rnb.set)[sample.subset]

	## FIXME: This might leads to duplicated identifiers and/or reordering of samples
	sample.labels<-abbreviate.names(sample.ids)

#	exp.intensity <- lapply(c("green" = "Green", "red" = "Red"), function(channel) {
#				ifelse(meta[[paste("Evaluate", channel)]] == "+", meta[["Expected Intensity"]], "Bgnd")
#			})

#	plot.names <- lapply(exp.intensity, function(exps) { paste(plot.names, exps, sep = ":\n") })

	plot.data <- lapply(names(intensities), function(channel) {
				channel.data <- intensities[[channel]]
				data.frame(
						"Intensity" = as.numeric(channel.data),
						"Sample" = factor(sapply(sample.ids, rep, times=nrow(channel.data)), as.character(sample.ids)))
			})
	if(is.null(name.prefix)){
		plot.name<-"NegativeControlBoxPlot"
	}else{
		plot.name<-paste("NegativeControlBoxPlot", name.prefix, sep="_")
	}
	if(writeToFile) {
		plot.file<-createReportPlot(plot.name, ...)
		grid.newpage()
	}

	## Define viewports and assign it to grid layout
	#grid.newpage()
	#pushViewport(viewport(layout = grid.layout(2,1)))

	plots<-list()

	for (i in 1:length(intensities)) {
		channel <- names(intensities)[i]
		plots[[i]]<-qplot(Sample, Intensity, data=na.omit(plot.data[[i]]), geom = "boxplot",
						main=paste("Negative control probes ", paste(channel, "channel"), sep=": "),
						ylab="Intensity",fill=I(channel))+
						scale_x_discrete(labels=sample.labels)+
						.control.plot.theme.samples
						#,vp=viewport(layout.pos.row=i, layout.pos.col=1))
	}

	## FIXME: Rewrite the plotting function without using the function arrangeGrob in gridExtra
	grb<-do.call(arrangeGrob, plots)

	grid.draw(grb)
	if(writeToFile) {
		off(plot.file)
		return(plot.file)
	}else{
		return(invisible(grb))
	}
}
#######################################################################################################################
#' rnb.plot.control.barplot
#'
#' Per-sample bar plots of Illumina HumanMethylation control probes
#'
#' @param rnb.set 		\code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}} object with valid
#' 						quality control information
#' @param probe 		exact id of the control probe consisting of the control probe type (see \code{\link{rnb.plot.control.boxplot}})
#' @param sample.subset an integer vector specifying the subset of samples for which the plotting should be performed
#' @param writeToFile 	flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file name with digits
#' @param name.prefix	in case \code{writeToFile} is \code{TRUE}, a \code{character} singleton specifying a prefix to the variable part of the image file names
#' @param verbose		if \code{TRUE} additional diagnostic output is generated
#' @param ... 			other arguments to \code{\link{createReportPlot}}
#'
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class
#' 			\code{\link{ggplot}} otherwise.
#'
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' control.meta.data <- rnb.get.annotation("controls450")
#' ctrl.probe<-paste0(unique(control.meta.data[["Target"]])[4], ".3")
#' print(ctrl.probe) # EXTENSION.3
#' rnb.plot.control.barplot(rnb.set.example, ctrl.probe)
#' }
#'
#' @author Pavlo Lutsik
#' @export
rnb.plot.control.barplot<-function(
		rnb.set,
		probe,
		sample.subset=1:length(samples(rnb.set)),
		writeToFile=FALSE,
		numeric.names=FALSE,
		name.prefix=NULL,
		verbose=FALSE,
		...)
{
	
	
	if(rnb.set@target=="probesEPIC" || rnb.set@target=="probesEPICv2" || rnb.set@target=="probesMMBC"){
		genome.assembly<-rnb.getOption("assembly") ## TODO: Improve genome build selection
		control.meta.data <- rnb.get.annotation(gsub("probes", "controls", rnb.set@target), assembly=c("probesEPIC"=genome.assembly, "probesEPICv2"="hg38", "probesMMBC"="mm10")[rnb.set@target])
		### TODO: Remove the following passage
		### for testing purposes only!
		#if(rnb.set@target=="probesEPIC"){
		#	control.meta.data<-rnb.update.controlsEPIC.enrich(control.meta.data)
		#}
		type=strsplit(probe,".", fixed=TRUE)[[1]][1]
		index=as.integer(strsplit(probe,".", fixed=TRUE)[[1]][2])
		if(!(type %in% rnb.infinium.control.targets(rnb.set@target))){
			rnb.error(c("Unrecognized control probe:", probe))
		}
		id<-subset(control.meta.data, Target==type & Index==index)$ID
		id.col<-"ID"
		target.col<-"Target"
		targets<-c(14,4,3,15,1:2,12:13,6,11)
	}else if(rnb.set@target=="probes450"){
		control.meta.data <- rnb.get.annotation("controls450")
		type=strsplit(probe,".", fixed=TRUE)[[1]][1]
		index=as.integer(strsplit(probe,".", fixed=TRUE)[[1]][2])
		if(!(type %in% rnb.infinium.control.targets(rnb.set@target))){
			rnb.error(c("Unrecognized control probe:", probe))
		}
		id<-subset(control.meta.data, Target==type & Index==index)$ID
		id.col<-"ID"
		target.col<-"Target"
		targets<-c(13,4,3,14,1:2,11:12,6)
	}else if(rnb.set@target=="probes27"){
		control.meta.data <- rnb.get.annotation("controls27")
		if(!probe %in% control.meta.data$Name){
			rnb.error(c("Unrecognized control probe:", probe))
		}
		id<-control.meta.data[probe==control.meta.data[["Name"]],"Address"]
		id.col<-"Address"
		target.col<-"Type"
		targets<-c(10,3,2,11,1,9,6)
	}

	#get intensities
	green<-qc(rnb.set)$Cy3[,sample.subset, drop=FALSE]
	red<-qc(rnb.set)$Cy5[,sample.subset, drop=FALSE]

	#get full set of probes
	full.probe.set<-control.meta.data[[id.col]][control.meta.data[[target.col]] %in% rnb.infinium.control.targets(rnb.set@target)[targets]]
	
	probe.name<-gsub(" ", ".", probe)

	if(! id %in% rownames(green)){
		
		plot.obj<-rnb.message.plot(sprintf("Box plot for control type %s not available", type))
		if(writeToFile){
			plot.name<-paste('ControlBarPlot', name.prefix, ifelse(numeric.names, 
							match(id, full.probe.set), 
							probe.name) , sep="_")
			plot.file<-createReportGgPlot(plot.obj, plot.name, ...)
			print(plot.obj)
			off(plot.file)
			return(plot.file)
		}
		return(plot.obj)
	}

	green<-green[as.character(id),]
	red<-red[as.character(id),]

#   if(is.null(dim(green))) green<-t(as.matrix(green))
#   if(is.null(dim(red))) red<-t(as.matrix(red))

	## get meta information

	if(rnb.set@target=="probes450" || rnb.set@target=="probesEPIC" || rnb.set@target=="probesMMBC" || rnb.set@target=="probesEPICv2"){
		meta<-subset(control.meta.data, ID==id)
	}else if(rnb.set@target=="probes27"){
		meta<-subset(control.meta.data, Address==id)
	}


	plot.name<-paste('ControlBarPlot', name.prefix, ifelse(numeric.names, 
					match(id, full.probe.set), 
					probe.name) , sep="_")
	if(writeToFile){
		plot.file<-createReportPlot(plot.name, ...)
		grid.newpage()
	}

	Samples<-factor(samples(rnb.set)[sample.subset], as.character(samples(rnb.set)[sample.subset]))
	## shorten sample names

	## FIXME: This might leads to duplicated identifiers and/or reordering of samples
	sample.labels<-abbreviate.names(samples(rnb.set)[sample.subset])

	#grid.newpage()
	# define viewports and assign it to grid layout
	#pushViewport(viewport(layout = grid.layout(2,1)))

	### plot green channel

	if(rnb.set@target=="probes450" || rnb.set@target=="probesEPIC" || rnb.set@target=="probesMMBC" || rnb.set@target=="probesEPICv2"){
		main_txt_grn<-paste(probe, meta[,"Description"], "green channel", if(meta[, "Evaluate Green"]=="+") meta[, "Expected Intensity"] else "Background", sep=": ")
		main_txt_red<-paste(probe, meta[,"Description"], "red channel",if(meta[, "Evaluate Red"]=="+") meta[, "Expected Intensity"] else "Background", sep=": ")
	}else{
		main_txt_grn<-paste(probe, sep="")
		main_txt_red<-paste(probe, sep="")
	}

	empty.probe<-sum(sapply(green, is.na))==length(green)

    if(empty.probe) {
		green<-rep(0, length(green))
	}
	
	green.plot<-ggplot(data = data.frame(Samples=factor(Samples, levels=Samples), Intensity=green), aes(x = Samples, y = Intensity)) +
			geom_bar(stat = "identity", position = "stack",fill = "green") + ggtitle(main_txt_grn) +
					scale_x_discrete(labels=sample.labels)+ylab("Intensity")+
					.control.plot.theme.samples

	if(empty.probe) {
		green.plot<-green.plot+annotate("text", label="No intensities found for this quality control probe", x=length(green)%/%2, y=0.5, size=5)
	}
	#print(green.plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))

	### plot red channel
	if(sum(sapply(red, is.na))==length(red)) empty.probe<-TRUE else empty.probe<-FALSE
	if(empty.probe) {
		red<-rep(0, length(red));
	}
	red.plot<-ggplot(data = data.frame(Samples=factor(Samples, levels=Samples), Intensity=red), aes(x = Samples, y = Intensity)) +
			geom_bar(stat = "identity", position = "stack",fill = "red") + ggtitle(main_txt_red) +
			scale_x_discrete(labels=sample.labels)+ylab("Intensity")+
			.control.plot.theme.samples

	if(empty.probe){
		red.plot<-red.plot+annotate("text", label="No intensities found for this quality control probe", x=length(red)%/%2, y=0.5, size=5)
	}
	#print(red.plot, vp=viewport(layout.pos.row=2, layout.pos.col=1))

	## FIXME: Rewrite the plotting function without using the function arrangeGrob in gridExtra
	grb<-arrangeGrob(green.plot, red.plot)
  grid.draw(grb)
  
  if(writeToFile) {
	  off(plot.file)
	  return(plot.file)
  }else{
	  return(invisible(grb))
  }

}

#######################################################################################################################

#' rnb.get.snp.matrix
#' 
#' Gets the matrix with "beta" values of SNP probes in the given dataset.
#' 
#' @param dataset       Microarray-based methylation dataset as an object of type inheriting \code{RnBeadSet}.
#' @param threshold.nas Threshold for removing probes (rows) and samples (columns). If more than this fraction of
#'                      \code{NA}s are present in the SNP matrix, the corresponding row or column is removed.
#' @return Two-dimensional \code{matrix} with the values at the SNP probes.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.get.snp.matrix <- function(dataset, threshold.nas = 1) {
	if (dataset@target %in% c("probes450", "probesEPIC", "probesMMBC", "probesEPICv2")) {
		result <- meth(dataset, row.names=TRUE)
		result <- result[grep("^rs", rownames(result)), , drop = FALSE]
	} else if (dataset@target == "probes27") {
		result <- HM27.snp.betas(qc(dataset))
	} else {
		stop("invalid value for dataset")
	}
	if (length(result) == 0) {
		stop("invalid value for dataset; no SNP probes found")
	}
	
	i <- which(apply(is.na(result), 1, mean) > threshold.nas)
	if (length(i) != 0) {
		if (length(i) == nrow(result)) {
			rnb.error("Not enough SNP data available (too many missing values per probe)")
		}
		result <- result[-i, , drop = FALSE]
		rnb.warning(paste(length(i), "probes ignored because their they contain too many NAs"))
	}
	i <- which(apply(is.na(result), 2, mean) > threshold.nas)
	if (length(i) != 0) {
		if (length(i) == ncol(result)) {
			rnb.error("Not enough SNP data available (too many missing values per sample)")
		}
		result <- result[, -i, drop = FALSE]
		rnb.warning(paste(length(i), "samples ignored because their SNP probes contain too many NAs"))
	}
	return(result)
}

#######################################################################################################################

#' rnb.plot.snp.boxplot
#'
#' Box plots of beta-values from the genotyping probes
#'
#' @param dataset     Dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}, or a matrix of
#'                    methylation beta values.
#' @param writeToFile Flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}.
#' @param ...         Additional named arguments passed to \code{\link{createReportPlot}}.
#' @return If \code{writeToFile} is \code{TRUE}: plot as an object of type \code{\linkS4class{ReportPlot}}. Otherwise:
#'         plot as an object of type \code{\link{ggplot}}.
#'
#' @export
#' @author Pavlo Lutsik
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.snp.boxplot(rnb.set.example)
#' }
rnb.plot.snp.boxplot <- function(dataset, writeToFile=FALSE, ...) {

	if (inherits(dataset, "RnBeadSet")) {
		dataset <- rnb.get.snp.matrix(dataset)
	}
	if (!(is.matrix(dataset) && is.numeric(dataset) && length(dataset) != 0)) {
		stop("invalid value for dataset")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}

	dframe <- data.frame(
		bv = as.numeric(dataset),
		SNP = factor(rep(rownames(dataset), ncol(dataset)), levels = rownames(dataset)))
	rplot <- qplot(SNP, bv, data = dframe, geom = "boxplot", main = "", ylab = expression(beta)) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0))

	if (writeToFile) {
		plot.file <- createReportPlot("SNPBoxplot", ...)
		print(rplot)
		off(plot.file)
		return(plot.file)
	}
	return(invisible(rplot))
}

#######################################################################################################################

#' rnb.plot.snp.barplot
#'
#' Bar plots of beta-values from the genotyping probes
#'
#' @param dataset       Dataset as an instance of \code{\linkS4class{RnBeadRawSet}} or \code{\linkS4class{RnBeadSet}}.
#'                      Alternatively, the dataset can be specified as a non-empty \code{matrix} containing the computed
#'                      beta values on the SNP probes.
#' @param probeID       Probe identifier. This must be one of \code{rownames(meth(dataset))}.
#' @param writeToFile   Flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}.
#' @param numeric.names if \code{TRUE} and \code{writeToFile} is \code{TRUE}substitute the plot options in the plot file
#'                      name with digits.
#' @param ...           Additional named arguments passed to \code{\link{createReportPlot}}.
#'
#' @return plot as an object of type \code{\linkS4class{ReportPlot}} if \code{writeToFile} is \code{TRUE} and of class
#' 			\code{\link{ggplot}} otherwise.
#'
#' @export
#' @author Pavlo Lutsik
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' samp<-samples(rnb.set.example)[1]
#' rnb.plot.snp.barplot(rnb.set.example, samp)
#' }
rnb.plot.snp.barplot <- function(dataset, probeID, writeToFile=FALSE, numeric.names=FALSE, ...) {

	if (inherits(dataset, "RnBeadSet")) {
		dataset <- rnb.get.snp.matrix(dataset)
	}
	if (!(is.matrix(dataset) && is.numeric(dataset) && length(dataset) != 0)) {
		stop("invalid value for dataset")
	}
	ids <- rownames(dataset)
	if (is.null(ids)) {
		stop("invalid value for dataset; missing row names")
	}
	if (is.null(colnames(dataset))) {
		stop("invalid value for dataset; missing column names")
	}
	if (!(is.character(probeID) && length(probeID) == 1 && isTRUE(probeID != ""))) {
		stop("invalid value for probeID")
	}
	if (!(probeID %in% ids)) {
		stop("invalid value for probeID; missing in dataset")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}
	if (!parameter.is.flag(numeric.names)) {
		stop("invalid value for numeric.names; expected TRUE or FALSE")
	}

	dframe <- data.frame(
		bv = dataset[probeID, ],
		sm = factor(colnames(dataset), levels = colnames(dataset)))
	rplot <- ggplot(dframe, aes_string(x = 'sm', y = 'bv')) + geom_bar(stat = "identity", position = "stack") + 
		scale_y_continuous(limits = c(0, 1)) + labs(title = probeID, x = NULL, y = expression(beta)) +
		.control.plot.theme.samples

	if(writeToFile){
		pfile <- ifelse(numeric.names, match(probeID, ids), gsub("[ |_]", ".", probeID))
		pfile <- createReportPlot(paste0('SNPBarPlot_', pfile), ...)
		print(rplot)
		off(pfile)
		return(pfile)
	}
	return(invisible(rplot))
}

#######################################################################################################################

#' rnb.plot.snp.heatmap
#'
#' Heatmap of beta values from genotyping probes.
#'
#' @param dataset     Dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}, or a matrix of
#'                    methylation beta values.
#' @param writeToFile Flag specifying whether the output should be saved as \code{\linkS4class{ReportPlot}}.
#' @param ...         Additional named arguments passed to \code{\link{createReportPlot}}. These are used only if
#'                    \code{writeToFile} is \code{TRUE}.
#' @return If \code{writeToFile} is \code{TRUE}, plot as an object of type \code{\linkS4class{ReportPlot}}. Otherwise,
#'         there is no value returned (invisible \code{NULL}).
#'
#' @export
#' @author Pavlo Lutsik
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' rnb.plot.snp.heatmap(rnb.set.example)
#' }
rnb.plot.snp.heatmap <- function(dataset, writeToFile = FALSE, ...) {

	if (inherits(dataset, "RnBeadSet")) {
		dataset <- rnb.get.snp.matrix(dataset)
	}
	if (!(is.matrix(dataset) && is.numeric(dataset) && length(dataset) != 0)) {
		stop("invalid value for dataset")
	}
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}

	if (writeToFile) {
		plot.file <- createReportPlot('SNPHeatmap', ...)
	}
	if (nrow(dataset) < 2 || ncol(dataset) < 2) {
		rnb.message.plot("Heatmap is not available due to insufficient data.")
	} else {
		sample.ids <- abbreviate.names(colnames(dataset))
		meth.colors <- get.methylation.color.panel()
		meth.breaks <- seq(0, 1, length.out = length(meth.colors) + 1L)
		suppressWarnings(heatmap.2(dataset, scale = "none", na.rm = FALSE,
				breaks = meth.breaks, col = meth.colors, trace = "none",
				margins = c(8, 8), labRow = rownames(dataset), labCol = sample.ids,
				density.info = "density", key.title = NA, key.xlab = expression(beta), key.ylab = "Density"))
	}
	if (writeToFile) {
		off(plot.file)
		return(plot.file)
	}
	return(invisible(NULL))
}

#######################################################################################################################

## point.whisker.ggplot
##
## Creates a ggplot object containing a point-and-whisker plot.
##
## @param dframe Table of values in the form of a \code{data.frame} containig at least the following columns:
##               \code{"Position"}, \code{"Average"}, \code{"Deviation"}.
## @param yrange Requested value range of the y axis.
## @param xlab Title of the x axis.
## @param ylab Title of the y axis.
##
## @author Yassen Assenov
point.whisker.ggplot <- function(dframe, yrange, xlab = "Position on the slide", ylab = "Methylation") {
	ggplot(dframe, aes_string(x = "Position", y = "Average")) +
		coord_cartesian(ylim = yrange) + scale_x_discrete(drop = FALSE) + ggplot2::geom_point(size = 3) +
		geom_errorbar(aes_string(ymin = "Average - Deviation", ymax = "Average + Deviation"), width = 0.4) +
		labs(x = xlab, y = ylab) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

#######################################################################################################################

#' rnb.plot.sentrix.distribution
#'
#' Creates a point-and-whisker plots showing beta value distributions at Sentrix positions for the given slide.
#'
#' @param rnb.set    HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param sentrix.id Slide number (Sentrix ID) as an \code{integer} or \code{character} singleton.
#' @return Generated point-and-whisker plot (an instance of \code{\link{ggplot}}) of mean methylations for the samples
#'         on the specified slide, or \code{FALSE} if the dataset is non-empty but does not contain samples on the
#'         given slide. If the provided dataset does not contain valid Sentrix ID and position information (or is an
#'         empty dataset), this method returns \code{NULL}.
#'
#'
#' @author Yassen Assenov
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' sid<-as.character(pheno(rnb.set.example)[["Sentrix_ID"]][1])
#' rnb.plot.sentrix.distribution(rnb.set.example,sid)
#' }
rnb.plot.sentrix.distribution <- function(rnb.set, sentrix.id) {
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!((is.character(sentrix.id) || is.integer(sentrix.id)) && length(sentrix.id) == 1 && (!is.na(sentrix.id)))) {
		stop("invalid value for sentrix.id")
	}
	dframe <- rnb.get.sentrix.data(rnb.set)
	if (is.null(dframe) || nrow(dframe) == 0) {
		return(NULL)
	}

	## Find the samples in the dataset on the specified slide
	s.indices <- which(dframe[, "Slide"] == sentrix.id) # this works both for integer and character sentrix.ids
	if (length(s.indices) == 0) {
		return(FALSE)
	}

	## Compute mean methylation and standard deviation
	m.methylations <- apply(as.matrix(meth(rnb.set)), 2, function(x) { c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)) })
	if(!is.matrix(m.methylations)){
		m.methylations<-matrix(m.methylations, ncol=1)
	}
	yrange <- range(m.methylations[1, ] - m.methylations[2, ], m.methylations[1, ] + m.methylations[2, ])
	m.methylations <- cbind(dframe[s.indices, ], data.frame(
			"Identifier" = colnames(meth(rnb.set))[s.indices],
			"Average" = m.methylations[1, ],
			"Deviation" = m.methylations[2, ]))
	m.methylations <- m.methylations[order(m.methylations[, "Position"]), ]

	## Create the plot
	point.whisker.ggplot(m.methylations, yrange)
}

#######################################################################################################################

#' rnb.plot.sentrix.distributions
#'
#' Creates one or more point-and-whisker plots showing beta value distributions at Sentrix positions.
#'
#' @param rnb.set HumanMethylation450K dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param fprefix File name prefix to be used in the generated plots. In order to ensure independence of the operating
#'                system, there are strong restrictions on the name of the file. See the documentation of
#'                \code{\link{createReportPlot}} for more information.
#' @param ...     Other arguments passed to \code{\link{createReportPlot}}. These can include the named parameters
#'                \code{report}, \code{width}, \code{height}, and others.
#' @return Point-and-whisker plot (an instance of \code{\linkS4class{ReportPlot}}), or a list of such plots - one per
#'         slide. If the provided dataset does not contain valid Sentrix ID and position information (or is an empty
#'         dataset), this method returns \code{NULL}.
#'
#' @details
#' If no additional parameters are specified, this function creates one PDF and one low-resolution PNG file for every
#' generated plot.
#'
#' @seealso \code{\link{rnb.plot.sentrix.distribution}} for creating a single plot for a specified slide number
#'
#' @author Yassen Assenov
#' @export
rnb.plot.sentrix.distributions <- function(rnb.set, fprefix = "sentrix_whisker", ...) {
	if (!inherits(rnb.set, "RnBeadSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.character(fprefix) && length(fprefix) == 1 && (!is.na(fprefix)))) {
		stop("invalid value for fprefix")
	}
	if (!grepl("^[A-Za-z0-9._-]+$", fprefix)) {
		stop("invalid value for fprefix")
	}
	dframe <- rnb.get.sentrix.data(rnb.set)
	if (is.null(dframe) || nrow(dframe) == 0) {
		return(NULL)
	}

	m.methylations <- apply(meth(rnb.set), 2, function(x) { c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)) })
	if(!is.matrix(m.methylations)){
		m.methylations<-matrix(m.methylations, ncol=1)
	}
	i.valid <- which(!is.na(m.methylations[2, ]))
	if (length(i.valid) == 0) {
		return(NULL)
	}
	m.methylations <- m.methylations[, i.valid, drop=FALSE]
	yrange <- range(m.methylations[1, ] - m.methylations[2, ], m.methylations[1, ] + m.methylations[2, ])
	m.methylations <- cbind(dframe[i.valid, ], data.frame(
			"Identifier" = samples(rnb.set)[i.valid],
			"Average" = m.methylations[1, ],
			"Deviation" = m.methylations[2, ]))
	m.methylations <- m.methylations[with(m.methylations, order(Slide, Position)), ]
	slides <- tapply(1:nrow(m.methylations), m.methylations[, "Slide"], identity)

	report.plots <- list()
	for (si in 1:length(slides)) {
		fname <- ifelse(length(slides) == 1, fprefix, paste(fprefix, si, sep = "_"))
		rplot <- createReportPlot(fname = fname, ...)
		print(point.whisker.ggplot(m.methylations[slides[[si]], ], yrange))
		report.plots[[names(slides)[si]]] <- off(rplot)
	}
	if (length(report.plots) == 1) {
		return(report.plots[[1]])
	}
	return(report.plots)
}

#######################################################################################################################

get.unified.scale<-function(intensities){

	ir<-range(as.numeric(intensities), na.rm=TRUE)

	return(c(0,ir[2]+1000))

}
#######################################################################################################################
HM27.snp.betas<-function(qc.int){
	cmd<-rnb.get.annotation("controls27")
	cmd.snp<-cmd[grep("Genotyping", cmd$Type),]
	snp.ids<-cmd.snp[["Address"]]
	cy3.snp<-qc.int[["Cy3"]][rownames(qc.int[["Cy3"]]) %in% snp.ids,]
	cy5.snp<-qc.int[["Cy5"]][rownames(qc.int[["Cy5"]]) %in% snp.ids,]

	rownames(cy3.snp)<-cmd.snp$Name
	rownames(cy5.snp)<-cmd.snp$Name

	meth.cy3<-cy3.snp[paste(HM27.CY3.SNP.PROBES, "B", sep="_"),]
	umeth.cy3<-cy3.snp[paste(HM27.CY3.SNP.PROBES, "A", sep="_"),]

	meth.cy5<-cy5.snp[paste(HM27.CY5.SNP.PROBES, "B", sep="_"),]
	umeth.cy5<-cy5.snp[paste(HM27.CY5.SNP.PROBES, "A", sep="_"),]

	snp.betas.cy3<-meth.cy3/(meth.cy3+umeth.cy3+1)
	snp.betas.cy5<-meth.cy5/(meth.cy5+umeth.cy5+1)

	snp.betas<-rbind(snp.betas.cy3, snp.betas.cy5)
	rownames(snp.betas)<-gsub("_B", "", rownames(snp.betas))
	return(snp.betas)
}
#######################################################################################################################
rnb.update.controlsEPIC.enrich <- function(control.probe.infos) {
	
	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.COLORS.GREEN <- c("Black", "Blue", "Cyan", "Green", "Lime", "Limegreen", "Skyblue")
	
	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.COLORS.RED <- c("Gold", "Orange", "Purple", "Red", "Tomato", "Yellow")
	
	## Add columns Evaluate Green and Evaluate Red
	control.probe.infos[["Evaluate Green"]] <- "-"
	control.probe.infos[["Evaluate Red"]] <- "-"
	i <- grep("^DNP", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("^Biotin", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Color"]] %in% CONTROL.COLORS.GREEN)
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Color"]] %in% CONTROL.COLORS.RED)
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("^NEGATIVE", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	control.probe.infos[["Evaluate Green"]] <- factor(control.probe.infos[["Evaluate Green"]], levels = c("-", "+"))
	control.probe.infos[["Evaluate Red"]] <- factor(control.probe.infos[["Evaluate Red"]], levels = c("-", "+"))
	
	## Add column Expected Intensity
	control.probe.infos[["Expected Intensity"]] <- as.character(NA)
	i <- control.probe.infos[, "Target"] %in% c("NEGATIVE", "TARGET REMOVAL", "RESTORATION")
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- c("BISULFITE CONVERSION II", "SPECIFICITY II", "EXTENSION", "NON-POLYMORPHIC")
	i <- control.probe.infos[, "Target"] %in% i
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- control.probe.infos[, "Target"] %in% paste("NORM", c("A", "C", "G", "T"), sep = "_")
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((High)|(20K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((Medium)|(5K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Medium"
	i <- grep("\\(Low\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Low"

	#i <- grep("\\((Bkg)|(5K)\\)$", control.probe.infos[, "Description"])
	i <- grep("\\(Bkg\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^BS Conversion I[- ]C", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^BS Conversion I[- ]U", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^GT Mismatch.+\\(PM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^GT Mismatch.+\\(MM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	control.probe.infos[["Expected Intensity"]] <- factor(control.probe.infos[["Expected Intensity"]])
	
	## Add column Sample-dependent
	control.probe.infos[["Sample-dependent"]] <-
			!(control.probe.infos[["Target"]] %in% CONTROL.TARGETS.SAMPLE.INDEPENDENT)
	
	control.probe.infos[["Index"]][order(control.probe.infos$Target)]<-unlist(sapply(sort(unique(control.probe.infos$Target)), function(target) 1:length(which(control.probe.infos$Target==target))))
	
	return(control.probe.infos)
}
