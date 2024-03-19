#######################################################################################
##		agePrediction.R							     ##
##	This file contains all the relevant functions for the		 	     ##
##	age prediction section of RnBeads.					     ##
#######################################################################################

#######################################################################################################################
#' rnb.section.ageprediction
#'
#' Adds a section for the age prediction part to the report. In this function, we check the options
#' \bold{\code{inference.age.prediction.training}} and \bold{\code{inference.age.prediction.predictor}}
#' and accordingly either the predefined predictor (in extdata), a newly trained predictor (in extdata) or a
#' formerly trained predictor that is stored in the path given by \bold{\code{inference.age.prediction.predictor}} is applied
#' to execute epigenetic age prediction.
#'
#' @param object           methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param report           report in which the age prediction section should be included. This must be an object of type \code{\linkS4class{Report}}.
#' @return	(possibly modified) report.
#'
#' @author Michael Scherer
#' @noRd

rnb.section.ageprediction <- function(object,report){
	title <- "Age Prediction"
	descr <- "Plots for the visualization of predicted ages based on the DNA methylation pattern."
	report <- rnb.add.section(report,title,descr)
	predictedAges <- pheno(object)$predicted_ages
	if(is.null(predictedAges)){
		object <- rnb.execute.age.prediction(object)
		predictedAges <- pheno(object)$predicted_ages
		if(!is.null(predictedAges)){
			logger.start("Adding Age Prediction Section to Report")
			report <- rnb.section.ageprediction(object,report)
			logger.completed()
		}else{
			logger.info("Age prediction was not successful")
			prediction_path <- rnb.getOption("inference.age.prediction.predictor")
			if(!file.exists(prediction_path)){
					text <- paste("The specified predictor does not exist at",prediction_path,", please check the option <em>inference.age.prediction.predictor</em> again.")
					report <- rnb.add.paragraph(report,text)
			}
			return(report)
		}
	}
	if(!is.null(predictedAges) && sum(is.na(predictedAges))==length(samples(object))){
		txt <- "The Age prediction was not successful, therefore no age prediction section was added to the report"
		logger.info("No Age Prediction Section Added")
		report <- rnb.add.paragraph(report,txt)
		return(report)
	}

	prediction_path <- rnb.getOption("inference.age.prediction.predictor")
	if(inherits(object,"RnBeadSet")){
		header <- "A new predictor was trained for the specified data set"
		if(is.null(prediction_path) || prediction_path==""){
			if(rnb.getOption("inference.age.prediction.training")){
					txt <- "No age information was available for the column age in the annotation, therefore the predefined predictor is used in the further calculation."
					report <- rnb.add.paragraph(report,txt)
			}
			header <- "A predefined predictor was used"
			if(assembly(object)=="hg38"){
				text <- "Currently there is no predictor available for Genome Build <em>hg38</em>. Prediction will be performed on a predictor trained on a hg19-data set."
				logger.warning(text)
				report <- rnb.add.paragraph(report,text)
			}
			if(dim(annotation(object))[1] < 30000){
			  prediction_path <- system.file(file.path("extdata", "predefined_predictor_27K.csv"), package="RnBeads")
	  		}else{
	  		  prediction_path <- system.file(file.path("extdata", "predefined_predictor_450K.csv"), package="RnBeads")
			}
		}else{
			if(!rnb.getOption("inference.age.prediction.training")){
				header <- "A user-specified predictor was used"
			}
		}
	}else if(inherits(object,"RnBiseqSet")){
		header <- "A new predictor was trained for the specified data set"
		if(is.null(prediction_path) || prediction_path==""){
			if(rnb.getOption("inference.age.prediction.training")){
					txt <- "No age information was available for the column age in the annotation, therefore the predefined predictor is used in the further calculation."
					report <- rnb.add.paragraph(report,txt)
			}
		  header <- "A predefined predictor was used"
		  if(assembly(object)=="hg19"){
		    prediction_path <- system.file(file.path("extdata", "predefined_predictor_Biseq_hg19.csv"), package="RnBeads")
		  }else{
		    prediction_path <- system.file(file.path("extdata", "predefined_predictor_Biseq_hg38.csv"), package="RnBeads")
		  }
		}else{
			if(!rnb.getOption("inference.age.prediction.training")){
				header <- "A user-specified predictor was used"
			}
		}
	}
	txt <- c(header, " and is available as a <a href=\"",prediction_path,"\">","comma-separated file</a>.")
	rnb.add.paragraph(report, txt)

	add.info <- function(stitle, ffunction, txt, actualAges, predictedAges,...) {
		result <- rnb.add.section(report, stitle, txt, level = 2)
		result <- ffunction(report, object,actualAges,predictedAges,...)
		rnb.status(c("Added", stitle))
		result
	}

	ph <- pheno(object)
	age <- rnb.getOption("inference.age.column")
	if(age %in% colnames(ph)){
		actualAges <- ph[,age]
		if(!is.null(actualAges)){
			options(warn=-1)
			if(is.factor(actualAges)){
					actualAges <- as.character(actualAges)
			}
			if(is.character(actualAges)){
					actualAges <- lapply(actualAges,convert.string.ages)
					if(length(actualAges)==0 || is.null(actualAges)){
						predictedAges <- ph$predicted_ages
						report <- add.age.histogram(report,predictedAges)
						logger.warning("Could not read age column")
						text <- "The age column has a format that can not be read by RnBeads. Please check the option <em>inference.age.column</em>."
						report <- rnb.add.paragraph(report,text)
						return(report)
					}
					set.na <- function(z){
						if(length(z)==0){
							NA
						}else{
							z
						}
					}
					actualAges <- unlist(lapply(actualAges,set.na))
			}
			actualAgesNumeric <- as.numeric(actualAges)
			if(sum(is.na(actualAgesNumeric)) != length(actualAges)){
				predictedAges <- ph$predicted_ages
				txt <- "Plotting annotated ages versus predicted ages and indicating different traits with different colors and different shapes."
				report <- add.info("Comparison Plot",add.agecomparison.plot,txt,actualAgesNumeric,predictedAges)
				txt <- "Plotting differences between predicted and annotated age for each sample."
				report <- add.info("Error Plot",add.combination.plot,txt,actualAgesNumeric,predictedAges)
				#txt <- "Plotting quantiles for the difference between predicted and annotated ages."
				#report <- add.info("Quantile Plot",add.quantile.plot,txt,actualAgesNumeric,predictedAges)
				# stratification in sample groups
				txt <- "Plotting predicted ages in different sample groups."
				report <- rnb.add.section(report,"Stratification Plot",txt,level=2)
				report <- add.stratification.plot(report,actualAgesNumeric,predictedAges,sample.groups=rnb.sample.groups(object))
				rnb.status(c("Added","Stratification Plot"))
			}else{
				if(is.factor(actualAges)){
					actualAges <- as.character(actualAges)
				}
				if(is.character(actualAges)){
					actualAges <- lapply(actualAges,convert.string.ages)
					if(length(actualAges)==0 || is.null(actualAges)){
						predictedAges <- ph$predicted_ages
						report <- add.age.histogram(report,predictedAges)
						logger.warning("Could not read age column")
						text <- "The age column has a format that can not be read by RnBeads. Please check the option <em>inference.age.column</em>."
						report <- rnb.add.paragraph(report,text)
						return(report)
					}
					set.na <- function(z){
						if(length(z)==0){
							NA
						}else{
							z
						}
					}
					actualAges <- unlist(lapply(actualAges,set.na))
					actualAgesNumeric <- as.numeric(actualAges)
					predictedAges <- ph$predicted_ages
					txt <- "Plotting annotated ages versus predicted ages and indicating different traits with different colors."
					report <- add.info("Comparison Plot",add.agecomparison.plot,txt,actualAgesNumeric,predictedAges)
					#txt <- "Plotting differences between predicted ages for each sample."
					report <- add.info("Error Plot",add.combination.plot,txt,actualAgesNumeric,predictedAges)
					#txt <- "Plotting quantiles for the difference between predicted and annotated ages."
					#report <- add.info("Quantile Plot",add.quantile.plot,txt,actualAgesNumeric,predictedAges)
				}else{
					txt <- "The age annotation column of the dataset has a format that can not be read."
					predictedAges <- ph$predicted_ages
					report <- add.age.histogram(report,predictedAges)
					report <- rnb.add.paragraph(report,txt)
				}
			}
		}
	}else{
		txt <- "There was no age annotation available for the data set, therefore no inference based on age prediction could be done."
		report <- rnb.add.paragraph(report,txt)
		report <- add.age.histogram(report,predictedAges)
	}
	options(warn=1)
	return(report)
}

#######################################################################################################################
#' rnb.step.ageprediction
#'
#' This function is the starting point for the report generation of the epigenetic age prediction. It call the functions
#' \code{rnb.section.ageprediction} and, if activated, \code{rnb.run.cross.validation}.
#'
#' @param object	a \code{\linkS4class{RnBSet}} object
#' @param report	report in which the age prediction section should be created. This must be an object of
#'              	type \code{\linkS4class{Report}}.
#'
#' @return	modified report object
#'
#' @author Michael Scherer
#' @noRd

rnb.step.ageprediction <- function(object,report){

	if(!inherits(object,"RnBSet")){
		stop("Supplied object is not of the class inheriting from RnBSet")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}

	logger.start("Adding Age Prediction Section to Report")
	report <- rnb.section.ageprediction(object,report)
	logger.completed()

	if(rnb.getOption("inference.age.prediction.cv")){
		if(rnb.getOption("inference.age.prediction.training")){
			report <- run.cross.validation(object,report)
		}else{
			logger.info("Could not run Cross-Validation without Training. Please set inference.age.prediction.training.")
			txt <- "Cross-Validation cannot be performed without training of a new predictor. Please set the option <em>inference.age.prediction.training</em> and run inference again."
			report <- rnb.add.paragraph(report,txt)
		}
	}

	return(report)
}

#######################################################################################################################
#' rnb.execute.age.prediction
#'
#' Performs age prediction by either the specified predictor in the option \code{inference.age.prediction.predictor}
#' or by the corresponding predefined predictor.
#'
#' @param object  a \code{\linkS4class{RnBSet}} object for which age prediction should
#' 			be performed
#'
#' @return	modified  \code{\linkS4class{RnBSet}} object
#'
#' @author Michael Scherer
#' @export

rnb.execute.age.prediction <- function(object){
	if(!inherits(object,"RnBSet")){
		stop("Supplied object is not of the class inheriting from RnBSet")
	}
	prediction_path <- rnb.getOption("inference.age.prediction.predictor")
	if(inherits(object,"RnBSet")){
			logger.start("Performing Age Prediction")
			rnbSet <- agePredictor(object,prediction_path)
			logger.completed()
	}
	return(rnbSet)
}

#######################################################################################################################
#' rnb.execute.training
#'
#' Trains a new age predictor on the specified data set and writes it to the given path. Elastic net regression is
#' to fit the input ages to the methylation values .
#'
#' @param object	a \code{\linkS4class{RnBSet}} object on which a new predictor should be created
#' @param path		path to which the predictor should be written out
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @author Michael Scherer
#' @export

rnb.execute.training <- function(object,path="",alpha=0.8){
	if(path==""){
		stop("A path to store the predictor has to be specified.")
	}else{
		logger.start("Training new age predictor")
		age <- rnb.getOption("inference.age.column")
		if(!(age %in% colnames(pheno(object)))){
			logger.warning("No Age information available, training not successful")
		}else{
			prediction_path <- trainPredictor(object,path,alpha=0.8)
		}
		logger.completed()
	}
}

#######################################################################################################################
#' add.agecomparison.plot
#'
#' Adds a plot that compares predicted ages by the age predictor with the annotated ages of the samples. The individual
#' points are marked in different colors and shapes depending on different samples characteristics as defined by
#' \code{rnb.sample.groups}.
#'
#' @param report	report in which the corresponding plot should be integrated.
#' @param object	methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param actualAges	annotated ages as found in the sample annotation sheet.
#' @param predictedAges	ages as predicted by the age prediction algorithm.
#'
#' @return	modified report with the age plots integrated.
#'
#' @author	Michael Scherer
#' @noRd

add.agecomparison.plot <- function(report, object, actualAges, predictedAges){

	descr <- "Comparing predicted ages by age predictor with annotated ages. Points that are labeled with their identifiers have a larger difference than 15 years between predicted and annotated age."

	ph <- pheno(object)
	Default <- rep(NA,dim(ph)[1])
	ph <- cbind(ph,Default=Default)
	naAges <- is.na(actualAges)
	actualAges <- actualAges[!naAges]
	predictedAges <- predictedAges[!naAges]
	notPredicted <- is.na(predictedAges)
	actualAges <- actualAges[!notPredicted]
	predictedAges <- predictedAges[!notPredicted]
	ph <- ph[!naAges,]
	ph <- ph[!notPredicted,]
	sample_names <- samples(object)
	sample_names <- sample_names[!naAges]
	sample_names <- sample_names[!notPredicted]

	traits <- names(rnb.sample.groups(object))
	# Remove non-ASCII characters
	traits <- traits[!grepl("[^ -~]",traits)]
	traits <- c("Default",traits)
	remove.age <- traits %in% rnb.getOption("inference.age.column")
	traits <- traits[!remove.age]

	report.plots <- list()
	if(length(traits)<=1){
		actualAges <- as.numeric(actualAges)
		report.plot <- createReportPlot("age_comparison_",report)
		data <- data.frame(actualAges,predictedAges,Sample=sample_names)
		diff <- abs(predictedAges-actualAges)
		plot <- ggplot(data,aes(x=actualAges,y=predictedAges),environment=environment())+geom_point()+geom_text(aes(label=ifelse(diff > 15, as.character(Sample),"")),hjust=0,vjust=0,size=3)+xlim(-1,100)+ylim(-1,100)+geom_abline(intercept=0,slope=1)+xlab("Annotated Ages")+ylab("Predicted Age")+theme(legend.title=element_blank())
		print(plot)
		report.plot <- off(report.plot)
		report.plots <- c(report.plots,report.plot)
	}else{
		report.plots.add <- foreach(i=1:length(traits),.combine=c) %dopar%{
	    trait <- traits[i]
		  report.plots.inner <- list()
			if(trait %in% colnames(ph)){
				placeholder <- ph[,trait]
				if(!is.null(placeholder)){
					placeholder <- as.factor(placeholder)
					for(secondTrait in traits){
						if(secondTrait %in% colnames(ph)){
							placeholder2 <- ph[,secondTrait]
							if(!is.null(placeholder2)){
								if(!all(is.na(placeholder)) && !all(is.na(placeholder2))){
									if((mean(nchar(as.character(placeholder)),na.rm=TRUE)<42) && (mean(nchar(as.character(placeholder2)),na.rm=TRUE)<42)){
										placeholder2 <- as.factor(placeholder2)
										trait <- gsub(" ","",trait)
										trait <- gsub("[[:punct:]]","",trait)
										trait <- gsub("[^ -~]","",trait)
										secondTrait <- gsub(" ","",secondTrait)
										secondTrait <- gsub("[[:punct:]]","",secondTrait)
										secondTrait <- gsub("[^ -~]","",secondTrait)
										report.plot <- createReportPlot(paste0("age_comparison_",trait,"_",secondTrait), report)
										if(length(placeholder)==0){
											placeholder <- rep(NA,length(actualAges))
										}
										if(length(placeholder2)==0){
											placeholder2 <- rep(NA,length(actualAges))
										}
										data <- data.frame(actualAges,predictedAges,placeholder,Sample=sample_names)
										diff <- abs(predictedAges-actualAges)
										cvalues <- rep(rnb.getOption("colors.category"),length.out=length(placeholder))
										ptvalues <- rep(rnb.getOption("points.category"),length.out=length(placeholder2))
										plot <- ggplot(data,aes(x=actualAges,y=predictedAges,group=placeholder),environment=environment())+geom_point(aes(colour=placeholder,shape=placeholder2))+geom_text(aes(label=ifelse(diff > 15, as.character(Sample),"")),hjust=0,vjust=0,size=3)+xlim(-1,100)+ylim(-1,100)+geom_abline(intercept=0,slope=1)+xlab("Annotated Ages")+ylab("Predicted Age")+scale_colour_manual(name=trait,na.value = ifelse(trait=="Default","black","#C0C0C0"), values = cvalues)+scale_shape_manual(name=secondTrait,na.value=ifelse(trait=="Default",20,1L),values=ptvalues)
										print(plot)
										report.plot <- off(report.plot)
										report.plots.inner <- c(report.plots.inner,report.plot)
									}
								}else{
placeholder2 <- as.factor(placeholder2)
										trait <- gsub(" ","",trait)
										trait <- gsub("[[:punct:]]","",trait)
										trait <- gsub("[^ -~]","",trait)
										secondTrait <- gsub(" ","",secondTrait)
										secondTrait <- gsub("[[:punct:]]","",secondTrait)
										secondTrait <- gsub("[^ -~]","",secondTrait)
										report.plot <- createReportPlot(paste0("age_comparison_",trait,"_",secondTrait), report)
										if(length(placeholder)==0){
											placeholder <- rep(NA,length(actualAges))
										}
										if(length(placeholder2)==0){
											placeholder2 <- rep(NA,length(actualAges))
										}
										data <- data.frame(actualAges,predictedAges,placeholder,Sample=sample_names)
										diff <- abs(predictedAges-actualAges)
										cvalues <- rep(rnb.getOption("colors.category"),length.out=length(placeholder))
										ptvalues <- rep(rnb.getOption("points.category"),length.out=length(placeholder2))
										plot <- ggplot(data,aes(x=actualAges,y=predictedAges,group=placeholder),environment=environment())+geom_point(aes(colour=placeholder,shape=placeholder2))+geom_text(aes(label=ifelse(diff > 15, as.character(Sample),"")),hjust=0,vjust=0,size=3)+xlim(-1,100)+ylim(-1,100)+geom_abline(intercept=0,slope=1)+xlab("Annotated Ages")+ylab("Predicted Age")+scale_colour_manual(name=trait,na.value = ifelse(trait=="Default","black","#C0C0C0"), values = cvalues)+scale_shape_manual(name=secondTrait,na.value=ifelse(trait=="Default",20,1L),values=ptvalues)
										print(plot)
										report.plot <- off(report.plot)
										report.plots.inner <- c(report.plots.inner,report.plot)
								}
							}
						}
					}
					return(report.plots.inner)
				}
			}
		}
		report.plots <- c(report.plots,report.plots.add)
	}
	traits <- gsub(" ","",traits)
	traits <- gsub("[^ -~]","",traits)
	traits <- gsub("[[:punct:]]","",traits)
	setting_names <- list("Compare first trait"=traits,"Compare second trait"=traits[])
	names(setting_names[[1]]) <- traits
	names(setting_names[[2]]) <- traits
	if(length(traits)==1){
		setting_names = list()
	}
	report <- rnb.add.figure(report, descr, report.plots,setting_names)

	return(report)
}


#######################################################################################################################
#' add.error.plot
#'
#' Adds a plot to the given report that describes the differences between predicted ages by the age prediction
#' algorithm and the ages as annotated in the sample annotation sheet(if available). Possible outliers are marked by
#' their sample names.
#'
#' @param report	report in which the corresponding plot should be integrated.
#' @param object	methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param actualAges	annotated ages as found in the sample annotation sheet.
#' @param predictedAges	ages as predicted by the age prediction algorithm.
#'
#' @return	modified report with the error plot integrated.
#'
#' @author	Michael Scherer
#' @noRd


add.error.plot <- function(report, object ,actualAges ,predictedAges){

	descr <- "Differences between predicted ages and annotated ages. The mean of the difference is shown in blue, +/- two times the Standard Deviation in red. Points that are labeled with their identifiers have a higher deviation from the mean difference than four times the standard deviation."

	ph <- pheno(object)

	diff <- predictedAges-actualAges
	mean <- mean(diff,na.rm=TRUE)
	sd <- sd(diff,na.rm=TRUE)
	deviance <- mean-diff
	deviance <- abs(deviance)
	data <- data.frame(row.names(ph),diff,deviance)
	colnames(data) <- c("Sample","Value","Deviance")
	count <- length(predictedAges)
	hline_mean <- data.frame(yint=mean,Measure="Mean")
	hline_sd <- data.frame(yint=c(mean+4*sd,mean-4*sd),Measure="4*Standard Deviation")
	x_width <- count + 0.1*count
	x_width <- round(x_width,0)

	plot <- ggplot(data,aes(x=order(sort(row.names(data))),y=Value,label=Sample),environment=environment())+geom_text(aes(label=ifelse(Deviance > 4*sd, as.character(Sample),"")),hjust=0,vjust=0,size=3)+geom_point()+geom_hline(data=hline_mean,aes(yintercept=yint,linetype=Measure),color="blue",show_guide=TRUE)+geom_hline(data=hline_sd,aes(yintercept=yint,linetype=Measure),color="red",show_guide=TRUE)+xlab("Sample Number")+ylab("Difference between predicted age and annotated age")+xlim(0,x_width)
	report.plot <- createReportPlot("age_prediction_error", report)
	print(plot)
	report.plot <- off(report.plot)

	report <- rnb.add.figure(report,descr,report.plot)

	return(report)
}

#######################################################################################################################
#' add.quantile.plot
#'
#' Adds a plot to the given report object, that shows the distribution of the differences between
#' predicted and annotated ages. Point inside the 1- and 99-quantile, respectively are marked
#' by their sample identifiers.
#'
#' @param report	report in which the corresponding plot should be integrated.
#' @param object	methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param actualAges	annotated ages as found in the sample annotation sheet.
#' @param predictedAges	ages as predicted by the age prediction algorithm.
#'
#' @return	modified report with the quantile plot integrated.
#'
#' @author	Michael Scherer
#' @noRd

add.quantile.plot <- function(report, object, actualAges, predictedAges){

	descr <- "Shown in red are the 1- and 99- percentiles for the difference between predicted ages and annotated ages. The black line indicates the median of this difference. The point labels are the sample identifiers given by the phenotypic table."

	ph <- pheno(object)

	diff <- predictedAges - actualAges
	data <- data.frame(Difference=diff,Sample=row.names(ph))
	q1 <- quantile(diff,0.01,na.rm=TRUE)
	q99 <- quantile(diff,0.99,na.rm=TRUE)
	density <- density(diff)
	density.frame <- data.frame(Difference=density$x,Density=density$y)
	median <- median(diff)

	plot <- ggplot(data,aes(x=Difference,y=0.0))+geom_density(aes(y=after_stat(density)))+geom_point(aes(label=Sample))+geom_text(aes(label=Sample),hjust=.5,vjust=1.5,size=3,colour='#BE1E2D')+geom_area(data=subset(density.frame,Difference>=q1 & Difference<=q99),aes(x=Difference,y=Density),fill='#75BC1C')+geom_area(data=subset(density.frame,Difference<=q1),aes(x=Difference,y=Density),fill='#BE1E2D')+geom_area(data=subset(density.frame,Difference>=q99),aes(x=Difference,y=Density),fill='#BE1E2D')+geom_vline(xintercept=median)+ylab("Density")
	report.plot <- createReportPlot("age_prediction_quantile", report)
	print(plot)
	report.plot <- off(report.plot)

	report <- rnb.add.figure(report,descr,report.plot)

	return(report)
}

#######################################################################################################################
#' add.combination.plot
#'
#' This function adds a plot that is the combination of two plots in one plot.The upper
#' panel shows the distribution and the quantiles for the differences between predicted and annotated
#' ages and the lower panel shows the individual errors for each sample.
#'
#' @param report	report in which the corresponding plot should be integrated.
#' @param object	methylation dataset as an object of type \code{\linkS4class{RnBSet}}.
#' @param actualAges	annotated ages as found in the sample annotation sheet.
#' @param predictedAges	ages as predicted by the age prediction algorithm.
#'
#' @return	modified report with the combination plot integrated.
#'
#' @author	Michael Scherer
#' @noRd

add.combination.plot <- function(report, object, actualAges,predictedAges){

	descr <- "Upper Panel: Shown is the distribution of the differences between predicted ages and annotated ages. The black line indicates the median of these differences. \n \\ Lower panel: Differences between predicted ages and annotated ages. The mean of the difference is shown as a dashed line, the 1- and 99-perecntile, respectively, as solid lines. Points that are labeled with their identifiers lay outside of the quantiles."

	ph <- pheno(object)

	actualAges <- as.numeric(actualAges)
	naAges <- is.na(actualAges)
	actualAges <- actualAges[!naAges]
	predictedAges <- predictedAges[!naAges]
	notPredicted <- is.na(predictedAges)
	actualAges <- actualAges[!notPredicted]
	ph <- ph[!naAges,]
	ph <- ph[!notPredicted,]
	sample_names <- samples(object)
	sample_names <- sample_names[!naAges]
	sample_names <- sample_names[!notPredicted]

	diff <- predictedAges - actualAges
	q1 <- quantile(diff,0.01,na.rm=TRUE)
	q99 <- quantile(diff,0.99,na.rm=TRUE)
	na <- is.na(diff)
	diff <- diff[!na]

	density <- density(diff)
	Sample <- rep(NA,length(density$x))
	Set <- rep("Density",length(density$y))
	density.frame <- data.frame(Sample,Difference=as.numeric(density$x),Density=as.numeric(density$y),Set)
	median <- median(diff)

	mean <- mean(diff,na.rm=TRUE)
	data <- data.frame(sample_names,as.numeric(diff))
	data <- cbind(data,as.numeric(row.names(data)))
	Set <- rep("Difference",length(diff))
	data <- cbind(data,Set)
	colnames(data) <- c("Sample","Difference","Density","Set")
	hline_mean <- data.frame(yint=mean,Measure="Mean")
	hline_quantiles <- data.frame(yint=c(q1,q99),Measure="1- and 99-quantiles")

	complete_data <- rbind(density.frame,data)
	cvalues <- rep(rnb.getOption("colors.category"))

	report.plots <- list()

	report.plot <- createReportPlot("combination_plot_Points",report)
	plot <- ggplot(complete_data,aes(x=Difference,y=Density,colour=Set))+geom_point()+facet_grid(Set~.,scales="free",drop=TRUE)+geom_text(aes(label=ifelse(Difference <= q1, as.character(Sample),""),hjust=.5),vjust=0,size=3,colour="black")+geom_text(aes(x=Difference,y=Density,label=ifelse(Difference >= q99, as.character(Sample),"")),hjust=.5,vjust=0,size=3,colour="black")+geom_vline(data=hline_mean,aes(xintercept=yint,linetype=Measure),show_guide=TRUE)+geom_vline(data=hline_quantiles,aes(xintercept=yint,linetype=Measure),show_guide=TRUE)+ylab("Sample Number / Density")+xlab("Difference between predicted age and annotated age")+scale_colour_manual(na.value = "#C0C0C0", values = cvalues, guide=FALSE)
	print(plot)
	report.plot <- off(report.plot)
	report.plots <- c(report.plots,report.plot)

	report.plot <- createReportPlot("combination_plot_Identifiers",report)
	plot <- ggplot(complete_data,aes(x=Difference,y=Density,colour=Set))+geom_point()+facet_grid(Set~.,scales="free",drop=TRUE)+geom_text(aes(label=as.character(Sample)),colour="black",size=3)+geom_vline(data=hline_mean,aes(xintercept=yint,linetype=Measure),show_guide=TRUE)+geom_vline(data=hline_quantiles,aes(xintercept=yint,linetype=Measure),show_guide=TRUE)+ylab("Sample Number / Density")+xlab("Difference between predicted age and annotated age")+scale_colour_manual(na.value = "#C0C0C0", values = cvalues, guide=FALSE)
	print(plot)
	report.plot <- off(report.plot)
	report.plots <- c(report.plots,report.plot)

	setting_names <- list("Sample representation"=c("Points","Identifiers"))
	names(setting_names[[1]]) <- c("Points","Identifiers")

	report <- rnb.add.figure(report,descr,report.plots,setting_names)

	return(report)
}

#######################################################################################
#' add.age.histogram
#'
#' This function creates an age histogram for the predicted ages in the case
#' where there is no age annotation available for comparison with predicted ages.
#'
#' @param report	report object to be modified.
#' @param ages		ages for which the histogram should be created.
#'
#' @return		modified report object with the age histogram.
#'
#' @author	Michael Scherer
#' @noRd

add.age.histogram <- function(report,ages){
	data <- data.frame(Age=ages)
	colors <- rnb.getOption("colors.category")
	gradient <- rnb.getOption("colors.gradient")
	plot <- ggplot(data,aes(x=Age))+geom_histogram(aes(fill=after_stat(count),y=after_stat(density)),binwidth=5) +geom_density(color=colors[2])+scale_fill_gradient(low=gradient[1],high=gradient[2],name="Count")+ylab("Density")
	report.plot <- createReportPlot("predicted_ages_histogram",report)
	print(plot)
	report.plot <- off(report.plot)

	descr <- "Age distributions of the predicted ages by the age prediction algorithm."
	report <- rnb.add.figure(report,descr,report.plot)

	return(report)
}

#######################################################################################
#' add.stratification.plot
#'
#' This function creates a plot comparing predicted ages and the difference between annotated and 
#' predicted ages in the different sample groups defined by \code{\link{rnb.sample.groups}}
#'
#' @param report	report object to be modified.
#' @param annotated.ages		annotated ages for the individuals.
#' @param predicted.ages  ages as predicted by the age prediction algorithm
#' @param sample.groups sample groups as defined by \code{\link{rnb.sample.groups}}
#'
#' @return		modified report object with the age stratification plot
#'
#' @author	Michael Scherer
#' @noRd

add.stratification.plot <- function(report,annotated.ages,predicted.ages,sample.groups){
  report.plots <- list()
  nsamples <- length(annotated.ages)
  report.plots <- foreach(c=1:length(sample.groups),.combine = c) %dopar%{
    sample.group <- sample.groups[[c]]
    cvalues <- rep(rnb.getOption("colors.category"),length.out=length(sample.group))
    cat <- rep(NA,nsamples)
    for(i in 1:length(sample.group)){
      cat[sample.group[[i]]] <- names(sample.group)[i]
    }
    trait <- gsub("[[:punct:]]","",names(sample.groups)[c])
    trait <- gsub("[^ -~]","",trait)
    trait <- gsub(" ","",trait)
    to.plot <- data.frame(Increase=predicted.ages-annotated.ages,Group=cat)
    plot <- ggplot(to.plot,aes(x=Group,y=Increase,fill=Group))+geom_boxplot()+scale_fill_manual(values=cvalues)+ylab("Predicted Age - Annotated Age")+
      theme(axis.text.x=element_text(angle=90,hjust=1))
    report.plot.1 <- createReportPlot(paste("age_stratification",trait,"increase",sep="_"),report)
    print(plot)
    report.plot.1 <- off(report.plot.1)
    to.plot <- data.frame(Predicted=predicted.ages,Group=cat)
    plot <- ggplot(to.plot,aes(x=Group,y=Predicted,fill=Group))+geom_boxplot()+scale_fill_manual(values=cvalues)+ylab("Predicted Age")+
      theme(axis.text.x=element_text(angle=90,hjust=1))
    report.plot.2 <- createReportPlot(paste("age_stratification",trait,"predicted",sep="_"),report)
    print(plot)
    report.plot.2 <- off(report.plot.2)
    to.plot <- data.frame(Annotated=annotated.ages,Group=cat)
    plot <- ggplot(to.plot,aes(x=Group,y=Annotated,fill=Group))+geom_boxplot()+scale_fill_manual(values=cvalues)+ylab("Annotated Age")+
      theme(axis.text.x=element_text(angle=90,hjust=1))
    report.plot.3 <- createReportPlot(paste("age_stratification",trait,"annotated",sep="_"),report)
    print(plot)
    report.plot.3 <- off(report.plot.3)
    return(c(report.plots,report.plot.1,report.plot.2,report.plot.3))
  }
  s.groups <- gsub("[[:punct:]]","",names(sample.groups))
  s.groups <- gsub("[^ -~]","",s.groups)
  s.groups <- gsub(" ","",s.groups)
  names(s.groups) <- s.groups
  p.types <- c("Predicted Ages - Annotated Ages","Predicted Ages","Annotates Ages")
  names(p.types) <- c("increase","predicted","annotated")
  s.names <- list(Group=s.groups,Type=p.types)
  descr <- "Predicted Ages are stratified over different sample groups."
  report <- rnb.add.figure(report,descr,unlist(report.plots),s.names)
  return(report)
}

#######################################################################################
#' age.transformation
#'
#' This function is used to transform the input ages, to account for the fact that
#' the human DNA methylation pattern changes more in the first phase of life than in the later.
#'
#' The function is directly taken from Steve Horvath's prediction algorithm
#' DNA methylation age of human tissues and cell types, Steve Horvath, Genome Biology, 2013
#'
#' @param x		input age about to be transformed
#' @param adult.age	threshold for deciding between 'fast' and 'slow' changes in DNA methylation
#'
#' @return	transformed age
#'
#' @author 	Steve Horvath
#' @noRd

age.transformation <- function(x,adult.age=20)
{
	x=(x+1)/(1+adult.age);
	y=ifelse(x<=1, log( x),x-1);
	y
}

#######################################################################################
#' age.anti.transformation
#'
#' This function the ages that were transformed by \code{age.transformation} back to actual ages.
#'
#' The function is directly taken from Steve Horvath's prediction algorithm
#' DNA methylation age of human tissues and cell types, Steve Horvath, Genome Biology, 2013
#'
#' @param x		input age about to be transformed back
#' @param adult.age	threshold for deciding between 'fast' and 'slow' changes in DNA methylation
#'
#' @return	age transformed back
#'
#' @author	Steve Horvath
#' @noRd

age.anti.transformation <- function(x,adult.age=20)
{
	ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age)
}

#######################################################################################
#' agePredictor
#'
#' This function is the core of the age prediction algorithm. It takes a \code{RnBSet} object
#' and loads a given predictor from a csv file to directly predict the age of the samples
#' by this predictor.
#'
#' @param rnbSet	a \code{RnBSet} object containing the relevant methylation infos.
#' @param path		path to a csv file in which a trained predictor should lay.
#'			DEFAULT:	No predictor is specified. In this case the corresponding
#'					predefined predictor is used.
#'
#' @return	modified \code{RnBSet} with two new columns in the phenotypic table (predicted_ages
#'		for the predicted epigenetic ages and age_increase for the difference between predicted
#'		and annnotated ages)
#'
#' @author	Michael Scherer
#' @noRd

agePredictor <- function(rnbSet, path=""){
	if(!is.null(path) && path != ""){
		if(!file.exists(path)){
			logger.info(c("The specified predictor does not exists at ",path))
	  		return(rnbSet)
	  	}
	}else{
		if(inherits(rnbSet,"RnBeadSet")){
			anno <- annotation(rnbSet)
			if(dim(anno)[1] > 30000){
				path <- system.file(file.path("extdata", "predefined_predictor_450K.csv"), package="RnBeads")
			}else{
				path <- system.file(file.path("extdata", "predefined_predictor_27K.csv"), package="RnBeads")
			}
		}else if(inherits(rnbSet,"RnBiseqSet")){
			if(assembly(rnbSet)=="hg19"){
				path <- system.file(file.path("extdata", "predefined_predictor_Biseq_hg19.csv"), package="RnBeads")
			}else{
				path <- system.file(file.path("extdata", "predefined_predictor_Biseq_hg38.csv"), package="RnBeads")
			}
		}
	}
	if(inherits(rnbSet,"RnBeadSet")){
		rnbSet <- agePredictorChip(rnbSet,path)
	}else if(inherits(rnbSet,"RnBiseqSet")){
		rnbSet <- agePredictorBiseq(rnbSet,path)
	}
	return(rnbSet)
}

#######################################################################################
#' agePredictorChip
#'
#' This function is the corresponding instance of the general age prediction algorithm for data
#' created by Illumina Infinium Bead Chips
#'
#' @param rnbSet	a \code{RnBeadSet} object containing the relevant methylation infos.
#' @param path		path to a csv file in which a trained predictor should exist.
#'
#' @return	modified \code{RnBeadSet} with two new columns in the phenotypic table (predicted_ages
#'		for the predicted epigenetic ages and age_increase for the difference between predicted
#'		and annnotated ages)
#'
#' @author	Michael Scherer
#' @noRd
agePredictorChip <- function(rnbSet, path){
	rnb.require("impute")
	options(warn=-1)
	coeffs <- read.csv(path)
	methData <- meth(rnbSet)
	anno <- annotation(rnbSet)
	ph <- pheno(rnbSet)
	row.names(methData) <- row.names(anno)
	intercept <- coeffs[1,2]
	coeffs <- coeffs[-1,]
	usedCpGs <- coeffs[,1]
	match <- match(sort(usedCpGs),usedCpGs)
	usedCpGs <- sort(usedCpGs)
	usedCoeffs <- coeffs[,2]
	usedCoeffs <- usedCoeffs[match]
	selectCpGs <- row.names(methData) %in% usedCpGs
	existingCpGs <- usedCpGs  %in% row.names(methData)
	selected <- methData[selectCpGs,]
	selected <- selected[sort(row.names(selected)),]
	options(warn=-1)
	dummy <- capture.output(selected <- (impute.knn(selected,colmax=1))$data)
	options(warn=1)
	selected <- t(selected)
	if(length(usedCoeffs) > dim(selected)[2]){
		usedCoeffs <- usedCoeffs[existingCpGs]
	}
	predictedAges <- intercept + selected%*%usedCoeffs
	predictedAges <- age.anti.transformation(predictedAges)
	predictedAges <- as.numeric(predictedAges)
	ages <- ph$predicted_ages
	increase <- ph$age_increase
	if(is.null(ages)){
		rnbSet <- addPheno(rnbSet,predictedAges,"predicted_ages")
	}
	if(is.null(increase) || any(is.na(increase))){
		age <- rnb.getOption("inference.age.column")
		if(age %in% colnames(ph)){
			actualAges <- ph[,age]
			nas <- is.na(actualAges)
			difference <- rep(NA,length(actualAges))
			actualAges <- actualAges[!nas]
			predictedAges <- predictedAges[!nas]
			if(!is.null(actualAges)){
				if(is.factor(actualAges)){
					actualAges <- as.character(actualAges)
				}
				if(is.character(actualAges)){
					actualAges <- unlist(lapply(actualAges,convert.string.ages))
				}
				if(length(actualAges)==0 || is.null(actualAges)){
					logger.warning("Could not read age column")
					return(rnbSet)
				}
				actualAges <- as.numeric(actualAges)
				difference[!nas] <- predictedAges-actualAges
				rnbSet <- addPheno(rnbSet,difference,"age_increase")
			}
		}
	}
	options(warn=1)
	return(rnbSet)
}

#######################################################################################
#' agePredictorBiseq
#'
#' This function is the corresponding instance of the general age prediction algorithm for data
#' created by a sequencing-based approach.
#'
#' @param rnbSet	a \code{RnBiseqSet} object containing the relevant methylation infos.
#' @param path		path to a csv file in which a trained predictor should exist.
#'
#' @return	modified \code{RnBiseqSet} with two new columns in the phenotypic table (predicted_ages
#'		for the predicted epigenetic ages and age_increase for the difference between predicted
#'		and annnotated ages)
#'
#' @author	Michael Scherer
#' @noRd

agePredictorBiseq <- function(rnbSet, path){
	rnb.require("impute")
	options(warn=-1)
	coeffs <- read.csv(path)
	methData <- meth(rnbSet)
	anno <- annotation(rnbSet)
	ph <- pheno(rnbSet)
	start <- anno$Start
	end <- anno$End
	chromosome <- anno$Chromosome
	names <- paste0(chromosome,"_",start,"_",end)
	row.names(methData) <- names
	intercept <- coeffs[1,2]
	coeffs <- coeffs[-1,]
	usedCpGs <- coeffs[,1]
	match <- match(sort(usedCpGs),usedCpGs)
	usedCpGs <- sort(usedCpGs)
	usedCoeffs <- coeffs[,2]
	usedCoeffs <- usedCoeffs[match]
	selectCpGs <- row.names(methData) %in% usedCpGs
	existingCpGs <- usedCpGs  %in% row.names(methData)
	selected <- methData[selectCpGs,]
	selected <- selected[sort(row.names(selected)),]
	options(warn=-1)
	dummy <- capture.output(selected <- (impute.knn(selected,colmax=1))$data)
	options(warn=1)
	selected <- t(selected)
	if(dim(selected)[1]==0){
		logger.info("Age prediction could not be performed; NAs introduced")
		predictedAges <- rep(NA,length(samples(rnbSet)))
	}else{
		if(length(usedCoeffs) > dim(selected)[2]){
			usedCoeffs <- usedCoeffs[existingCpGs]
		}
		na_coeffs <- is.na(usedCoeffs)
		usedCoeffs <- usedCoeffs[!na_coeffs]
		selected <- selected[,!na_coeffs]
		predictedAges <- intercept + selected%*%usedCoeffs
		predictedAges <- age.anti.transformation(predictedAges)
		predictedAges <- as.numeric(predictedAges)
	}
	ages <- ph$predicted_ages
	increase <- ph$age_increase
	if(is.null(ages)){
		rnbSet <- addPheno(rnbSet,predictedAges,"predicted_ages")
	}
	if(is.null(increase) || any(is.na(increase))){
		age <- rnb.getOption("inference.age.column")
		if(age %in% colnames(ph)){
			actualAges <- ph[,age]
			nas <- is.na(actualAges)
			difference <- rep(NA,length(actualAges))
			actualAges <- actualAges[!nas]
			predictedAges <- predictedAges[!nas]
			if(!is.null(actualAges)){
				if(is.factor(actualAges)){
					actualAges <- as.character(actualAges)
				}
				if(is.character(actualAges)){
					actualAges <- unlist(lapply(actualAges,convert.string.ages))
				}
				if(length(actualAges)==0 || is.null(actualAges)){
					logger.warning("Could not read age column")
					return(rnbSet)
				}
				actualAges <- as.numeric(actualAges)
				difference[!nas] <- predictedAges-actualAges
				rnbSet <- addPheno(rnbSet,difference,"age_increase")
			}
		}
	}
	options(warn=1)
	return(rnbSet)
}

#######################################################################################
#' trainPredictor
#'
#' This function is the starting point for training a new predictor based on the the input dataset.
#'
#' @param rnbSet	An \code{RnBSet} object containing the methylation info on which
#'			the new predictor should be trained
#' @param data.dir	directory to which the resulting predictor should be written out
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @return	absolute path to the corresponding predictor. The path is also set as the
#'		option \code{inference.age.prediction.predictor}
#'
#' @author	Michael Scherer
#' @noRd

trainPredictor <- function(rnbSet,data.dir,alpha=0.8){
	if(!file.exists(data.dir)){
		stop("The specified directory does not exist!")
	}
	age <- rnb.getOption("inference.age.column")
	if(inherits(rnbSet,"RnBeadSet")){
		if(!(age %in% colnames(pheno(rnbSet)))){
			logger.warning("No Age information available, training not successful.")
			return("")
		}
		path <- simpleGlmnet(rnbSet,file.path(data.dir,"trained_predictor.csv"),alpha=alpha)
		if(!is.null(path)&&path!=""){
			rnb.options(inference.age.prediction.predictor=path)
		}else{
			logger.warning("No new predictor created, deault predictor remains.")
		}
	}else if(inherits(rnbSet,"RnBiseqSet")){
		if(!(age %in% colnames(pheno(rnbSet)))){
			logger.warning("No Age information available, training not successful.")
			return("")
		}
		path <- simpleGlmnetBiseq(rnbSet,file.path(data.dir,"trained_predictor_Biseq.csv"),alpha=alpha)
		if(!is.null(path)&&path!=""){
			rnb.options(inference.age.prediction.predictor=path)
		}else{
			logger.warning("No new predictor created, deault predictor remains.")
		}
	}
	return(path)
}

#######################################################################################
#' simpleGlmnet
#'
#' This function fits the metylation data to the available ages by fitting
#' a generalized linear model (alpha parameter = 0.8). The lambda parameter is chosen by
#' 10-fold cross-validation. This function was created for array-based data.
#'
#' @param trainRnBSet	a \code{RnBeadSet} object containing the methylation info on which
#'			the new predictor should be trained
#' @param filePath	path to which the new predictor should be written out
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @return	absolute path to the corresponding predictor. \code{NULL} if the function was
#'		unable to create such an predictor.
#'
#' @author	Michael Scherer
#' @noRd

simpleGlmnet <- function(trainRnBSet,filePath="",alpha=0.8){
	rnb.require("glmnet")
	rnb.require("impute")
	methData <- meth(trainRnBSet)
	anno <- annotation(trainRnBSet)
	ph <- pheno(trainRnBSet)
	age <- rnb.getOption("inference.age.column")
	ages <- ph[,age]
	naAges <- is.na(ages)
	if(!is.null(ages)){
		if(is.factor(ages)){
			ages <- as.character(ages)
		}
		if(is.character(ages)){
			ages <- unlist(lapply(ages,convert.string.ages))
		}
		if(length(ages)==0){
			return(NULL)
		}
		ages <- as.numeric(ages)
		ages <- age.transformation(ages)
		methData <- methData[,!naAges]
		ages <- ages[!naAges]
		row.names(methData) <- row.names(anno)
		Xchrom <- is.element(anno$Chromosome,"chrX")
		methData <- methData[!Xchrom,]
		anno <- anno[!Xchrom,]
		Ychrom <- is.element(anno$Chromosome,"chrY")
		methData <- methData[!Ychrom,]
		anno <- anno[!Ychrom,]
		#Imputing of missing values is performed by a k-nearest-neighbors approach
		options(warn=-1)
		dummy <- capture.output(methData <- (impute.knn(methData,colmax=1))$data)
		options(warn=1)
		methData <- t(methData)
		missingAges <- is.na(ages)
		methData <- methData[!missingAges,]
		ages <- ages[!missingAges]
		cv <- cv.glmnet(methData,ages,parallel=TRUE)
		model <- glmnet(methData,ages,alpha=alpha,lambda=cv$lambda.min)
		coeffs <- as.matrix(coef(model))
		names <- row.names(coeffs)
		non_zero <- which(coeffs!=0)
		coeffs <- coeffs[non_zero]
		names <- names[non_zero]
		names(coeffs) <- names
		if(writePredictorToCsv(coeffs,filePath)){
			return(filePath)
		}
	}
	return(NULL)
}

#######################################################################################
#' simpleGlmnetBiseq
#'
#' This function fits the metylation data to the available ages by fitting
#' a generalized linear model (alpha parameter = 0.8). The lambda parameter is chosen by
#' 10-fold cross-validation. This function was created for sequencing-based data.
#'
#' @param trainRnBSet	a \code{RnBiseqSet} object containing the methylation info on which
#'			the new predictor should be trained
#' @param filePath	path to which the new predictor should be written out
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @return	absolute path to the corresponding predictor. \code{NULL} if the function was
#'		unable to create such an predictor.
#'
#' @author	Michael Scherer
#' @noRd

simpleGlmnetBiseq <- function(trainRnBSet,filePath="",alpha=0.8){
	rnb.require("glmnet")
	methData <- meth(trainRnBSet)
	anno <- annotation(trainRnBSet)
	ph <- pheno(trainRnBSet)
	age <- rnb.getOption("inference.age.column")
	ages <- ph[,age]
	naAges <- is.na(ages)
	if(!is.null(ages)){
		if(is.factor(ages)){
			ages <- as.character(ages)
		}
		if(is.character(ages)){
			ages <- unlist(lapply(ages,convert.string.ages))
		}
		if(length(ages)==0){
			return(NULL)
		}
		ages <- as.numeric(ages)
		ages <- age.transformation(ages)
		methData <- methData[,!naAges]
		ages <- ages[!naAges]
		Xchrom <- is.element(anno$Chromosome,"chrX")
		methData <- methData[!Xchrom,]
		anno <- anno[!Xchrom,]
		Ychrom <- is.element(anno$Chromosome,"chrY")
		methData <- methData[!Ychrom,]
		anno <- anno[!Ychrom,]
		start <- anno$Start
		end <- anno$End
		chromosome <- anno$Chromosome
		names <- paste0(chromosome,"_",start,"_",end)
		row.names(methData) <- names
		options(warn=-1)
		methData <- imputeBiseq(methData)
		options(warn=1)
		methData <- t(methData)
		missingAges <- is.na(ages)
		methData <- methData[!missingAges,]
		ages <- ages[!missingAges]
		cv <- cv.glmnet(methData,ages)
		model <- glmnet(methData,ages,alpha=alpha,lambda=cv$lambda.min)
		coeffs <- as.matrix(coef(model))
		names <- row.names(coeffs)
		non_zero <- which(coeffs!=0)
		coeffs <- coeffs[non_zero]
		names <- names[non_zero]
		names(coeffs) <- names
		if(writePredictorToCsv(coeffs,filePath)){
			return(filePath)
		}
	}
	return(NULL)
}

#######################################################################################
#' run.cross.validation
#'
#' This function performs 10-fold cross validation to estimate the performance of a
#' newly trained predictor. If \code{parallel.isEnabled()}, the function perfoms cross
#' validation in parallel. The function adds a table to the specified \code{report} containing
#' the result of the 10-fold cross validation.
#'
#' @param rnbSet	a \code{RnBSet} object containing the methylation info and ages on which
#'			the new predictor should be trained
#' @param report	report to which the table should be added
#' @param alpha alpha parameter used in the elastic net regression
#'
#' @return	modified report object
#'
#' @author	Michael Scherer
#' @export

run.cross.validation <- function(rnbSet,report,alpha=0.8){
	logger.start("10-fold Cross Validation")
	if(length(samples(rnbSet))<30){
		txt <- "ATTENTION: Cross-validated correlation result might be misleading, since there are less than 3 samples per fold."
		logger.warning(txt)
		report <- rnb.add.paragraph(report,txt)
	}
	txt <- "Result of the 10-fold cross validation on the given dataset. The corresponding error measures are: Correlation between predicted and annotated ages, Mean absolute deviation and Median absolute deviation"
	cvalues <- rep(rnb.getOption("colors.category"))
	descr <- "Boxplot for the two error measures Mean and Median absolute Error. Each Boxplot consists of 10 different points for each cross-validation fold, respectively."
	if(inherits(rnbSet,"RnBeadSet")){
		result <- cv.array(rnbSet,alpha = alpha)
		if(is.null(result)){
			text <- "Problem when performing age prediction."
			report <- rnb.add.paragraph(report,text)
			return(report)
		}
		report <- rnb.add.table(report,result,tcaption=txt)
		means <- result["Mean",]
		means <- t(means)
		medians <- result["Median",]
		medians <- t(medians)
		toPlot <- data.frame(Mean=means,Median=medians)
		toPlot <- melt(toPlot,id=c())
		colnames(toPlot) <- c("Measure","Error")
		plot <- ggplot(toPlot,aes(x=Measure,y=Error,fill=Measure))+geom_boxplot()+scale_fill_manual(values=cvalues)+ylab("Error [years]")
		report.plot <- createReportPlot("cv_error_boxplot",report)
		print(plot)
		report.plot <- off(report.plot)
		report <- rnb.add.figure(report,descr,report.plot)
	}else if(inherits(rnbSet,"RnBiseqSet")){
		result <- cv.biseq(rnbSet,alpha=alpha)
		if(is.null(result)){
			text <- "Problem when performing age prediction."
			report <- rnb.add.paragraph(report,text)
			return(report)
		}
		report <- rnb.add.table(report,result,tcaption=txt)
		means <- result["Mean",]
		means <- t(means)
		medians <- result["Median",]
		medians <- t(medians)
		toPlot <- data.frame(Mean=means,Median=medians)
		toPlot <- melt(toPlot,id=c())
		colnames(toPlot) <- c("Measure","Error")
		plot <- ggplot(toPlot,aes(x=Measure,y=Error,fill=Measure))+geom_boxplot()+scale_fill_manual(values=cvalues)+ylab("Error [years]")
		report.plot <- createReportPlot("cv_error_boxplot",report)
		print(plot)
		report.plot <- off(report.plot)
		report <- rnb.add.figure(report,descr,report.plot)
	}
	logger.completed()
	return(report)
}

######################################################################################
#' general.cv
#'
#' This functions performs k-fold-cross-validation on the predictor with methylation
#' data as input and ages to be trained on
#'
#' @param	fitFunction	a function that fits a predictor from training data to
#'				predict age from methylation data (e.g \code{simpeGlmnetEvaluate})
#' @param	ages		the ages to be trained on
#' @param	methData	input methylation matrix
#' @param	k		the fold parameter
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @return	a matrix that contains the summarized quality measures for
#'		the predictor which are:
#'			$cor[k+1]:	the mean correlation between predicted
#'					ages and actual ages
#'			$cor[1:k]:	individual correlations between predicted
#'					ages and actual ages for each fold
#'			$mean[k+1]:	the mean of the mean absolute deviation
#'					between predicted ages and actual ages
#'			$mean[1:k]:	the individual mean absolute deviations for
#'					each fold
#'			$median[k+1]:	the mean of the median absolute deviation
#'			$median[1:k]:	the individual median absolute devations
#'					for each fold
#'
#' @author	Michael Scherer
#' @noRd

general.cv <- function(fitFunction,ages,methData,k=10,alpha=0.8){
	nSamples <- length(ages)
	size <- nSamples%/%k
	count <- 1
	ret <- foreach(i=seq(1,k*size,by = size),.combine='cbind',.packages=c("RnBeads","impute","glmnet"),.export=c("age.transformation","age.anti.transformation","simpleGlmnetEvaluate")) %dopar% {
		logger.info(as.character(count))
		choose <- rep(FALSE,nSamples)
		choose[i:(i+size-1)] <- TRUE
		testSet <- methData[,choose]
		testAges <- ages[choose]
		notChosen <- !choose
		trainSet <- methData[,notChosen]
		trainAges <- ages[notChosen]
		predictor <- fitFunction(trainSet,trainAges,alpha)
		predictedAges <- predictor(testSet)
		cor <- cor(predictedAges,testAges)
		cor <- round(cor,2)
		testAges <- age.anti.transformation(testAges)
		mean <- mean(abs(predictedAges-testAges))
		mean <- round(mean,2)
		median <- median(abs(predictedAges-testAges))
		median <- round(median,2)
		column <- c(cor,mean,median)
		count <- count+1
		column
	}
	lastCol <- c(mean(ret[1,1:10]),mean(ret[2,1:10]),mean(ret[3,1:10]))
	ret <- cbind(ret,lastCol)
	colnames(ret) <- c(paste("Fold",seq(1:10)),"Mean")
	row.names(ret) <- c("Correlation","Mean","Median")
	return(ret)
}

###################################################################################
#' cv.array
#'
#' This function randomly changes the order of the samples in the data set and then
#' calls the \code{general.cv} to perform cross-validation for a RnBSet derived from
#' BeadChip arrays
#'
#' @param	rnbSet	\code{RnBeadSet} object on which the cross validation should be perfomed
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @return	the result of the cross validation in a data.frame format
#'
#' @author	Michael Scherer
#' @noRd

cv.array <- function(rnbSet,alpha=0.8){
	rnb.require("impute")
	ph <- pheno(rnbSet)
	age <- rnb.getOption("inference.age.column")
	if(age %in% colnames(ph)){
		ages <- ph[,age]
		if(is.factor(ages)){
			ages <- as.character(ages)
		}
		if(is.character(ages)){
			ages[!is.na(ages)] <- unlist(lapply(ages[!is.na(ages)],convert.string.ages))
		}
		ages <- as.numeric(ages)
		if(length(ages)==0 || is.null(ages) || sum(is.na(ages)) == length(ages)){
			logger.warning("Could not read age column")
			return(NULL)
		}
		ages <- age.transformation(ages)
		missingAges <- is.na(ages)
		ages <- ages[!missingAges]
		rnbSet <- remove.samples(rnbSet,missingAges)
		set.seed(Sys.time())
		samples <- samples(rnbSet)
		sampled <- sample(samples)
		methData <- meth(rnbSet)
		match <- match(sampled,samples)
		methData <- methData[,match]
		ages <- ages[match]
		anno <- annotation(rnbSet)
		row.names(methData) <- row.names(anno)
		Xchrom <- is.element(anno$Chromosome,"chrX")
		methData <- methData[!Xchrom,]
		anno <- anno[!Xchrom,]
		Ychrom <- is.element(anno$Chromosome,"chrY")
		methData <- methData[!Ychrom,]
		anno <- anno[!Ychrom,]
		options(warn=-1)
		dummy <- capture.output(methData <- (impute.knn(methData,colmax=1))$data)
		options(warn=1)
		rm(dummy)
		result <- general.cv(simpleGlmnetEvaluate,ages,methData,alpha=alpha)
		result <- as.data.frame(result)
		return(result)
	}
	return(NULL)
}

###################################################################################
#' cv.biseq
#'
#' This function randomly changes the order of the samples in the data set and then
#' calls the \code{general.cv} to perform cross-validation for a RnBSet derived from
#' bisulfite sequencing
#'
#' @param	rnbSet	\code{RnBiseqSet} object on which the cross validation should be perfomed
#' @param alpha     alpha parameter used in the elastic net regression
#'
#' @return	the result of the cross validation in a data.frame format
#'
#' @author	Michael Scherer
#' @noRd

cv.biseq <- function(rnbSet,alpha=0.8){
	ph <- pheno(rnbSet)
	age <- rnb.getOption("inference.age.column")
	if(age %in% colnames(ph)){
		ages <- ph[,age]
		if(is.factor(ages)){
			ages <- as.character(ages)
		}
		if(is.character(ages)){
			ages[!is.na(ages)] <- unlist(lapply(ages[!is.na(ages)],convert.string.ages))
		}
		ages <- as.numeric(ages)
		if(length(ages)==0 || is.null(ages) || sum(is.na(ages)) == length(ages)){
			logger.warning("Could not read age column")
			return(NULL)
		}
		missingAges <- is.na(ages)
		ages <- ages[!missingAges]
		rnbSet <- remove.samples(rnbSet,missingAges)
		ages <- age.transformation(ages)
		set.seed(Sys.time())
		samples <- samples(rnbSet)
		sampled <- sample(samples)
		methData <- meth(rnbSet)
		match <- match(sampled,samples)
		methData <- methData[,match]
		ages <- ages[match]
		anno <- annotation(rnbSet)
		start <- anno$Start
		end <- anno$End
		chromosome <- anno$Chromosome
		names <- paste0(chromosome,"_",start,"_",end)
		row.names(methData) <- names
		Xchrom <- is.element(anno$Chromosome,"chrX")
		methData <- methData[!Xchrom,]
		anno <- anno[!Xchrom,]
		Ychrom <- is.element(anno$Chromosome,"chrY")
		methData <- methData[!Ychrom,]
		anno <- anno[!Ychrom,]
		methData <- imputeBiseq(methData)
		result <- general.cv(simpleGlmnetEvaluate,ages,methData,alpha=alpha)
		result <- as.data.frame(result)
		return(result)
	}
	return(NULL)
}

###################################################################################
#' simpleGlmnetEvaluate
#'
#' This function is needed to perform cross-validation. In contrast to simpleGlmnet
#' it does not write the predictor to a csv-file but returns a prediction
#' function that is to be applied in each fold.
#'
#' @param	methData input methylation data as a data.frame for age prediction
#' @param	ages	  training ages
#' @param alpha   alpha parameter used in the elastic net regression
#'
#' @return	the age prediction function to be applied in each fold
#'
#' @author	Michael Scherer
#' @noRd

simpleGlmnetEvaluate <- function(methData,ages,alpha=0.8){
	methData <- t(methData)
	missingAges <- is.na(ages)
	methData <- methData[!missingAges,]
	ages <- ages[!missingAges]
	cv <- cv.glmnet(methData,ages,parallel=TRUE)
	model <- glmnet(methData,ages,alpha=alpha,lambda=cv$lambda.min)
	finalModel <- createPredictor(model,lambda=cv$lambda.min)
	return(finalModel)
}

###################################################################################
#' createPredictor
#'
#' This function is needed to perform cross-validation. It creates a prediction
#' function in contrast to writing the predictor to a csv-file.
#'
#' @param	linearModel	an output of glmnet from which the predictor should
#'			     be created
#'
#' @return	the age prediction function
#'
#' @author	Michael Scherer
#' @noRd

createPredictor <- function(linearModel,lambda=NULL){
	ret <- function(rnbSet){
		if(inherits(rnbSet,"RnBSet")){
			methData <- meth(rnbSet)
			coeffs <- as.matrix(coef(linearModel))
			intercept <- coeffs[1]
			anno <- annotation(rnbSet)
			coeffs <- coeffs[-1]
			usedCpGs <- row.names(coeffs)
			match <- match(sort(usedCpGs),usedCpGs)
			usedCpGs <- sort(usedCpGs)
			usedCoeffs <- coeffs
			usedCoeffs <- usedCoeffs[match]
			selectCpGs <- row.names(methData) %in% usedCpGs
			existingCpGs <- usedCpGs  %in% row.names(methData)
			selected <- methData[selectCpGs,]
			selected <- selected[sort(row.names(selected)),]
			options(warn=-1)
			dummy <- capture.output(selected <- (impute.knn(selected,colmax=1))$data)
			options(warn=1)
			selected <- t(selected)
			if(length(usedCoeffs) > dim(selected)[2]){
				usedCoeffs <- usedCoeffs[existingCpGs]
			}
			predictedAges <- intercept + selected%*%usedCoeffs
			predictedAges <- age.anti.transformation(predictedAges)
			return(predictedAges)
		}else{
			methData <- rnbSet
			methData <- t(methData)
			methData <- as.matrix(methData)
			predictedAges <- predict(linearModel,methData,s=lambda)
			predictedAges <- as.vector(predictedAges)
			predictedAges <- age.anti.transformation(predictedAges)
			return(predictedAges)
		}
	}
	return(ret)
}


#######################################################################################
#' writePredictorToCsv
#'
#' This function writes the coefficients of a linear model into a csv file.
#'
#' @param linearModel	linear Model which should be stored as a csv file
#' @param path		path in which the new predictor should be written
#'
#' @return	TRUE if the writing was successful
#'
#' @author	Michael Scherer
#' @noRd

writePredictorToCsv <- function(linearModel,path){
	coeffs <- as.matrix(linearModel)
	write.csv(coeffs,path)
	return(TRUE)
}

#######################################################################################
#' imputeBiseq
#'
#' Since the k-nearest-neighboring approach does not converge for all CpGs present in
#' sequencing data, a simple mean imputation is performed for this kind of data.
#'
#' @param methData	methylation data as a data frame or matix to be imputed
#'
#' @return	the changed methylation data frame
#'
#' @author	Michael Scherer
#' @noRd

imputeBiseq <- function(methData){
	samples <- dim(methData)[2]
	nas <- is.na(methData)
	newMean <- function(x){
		mean(x,na.rm=TRUE)
	}
	means <- apply(methData,1,newMean)
	for(i in (1:(dim(nas)[1]))){
		methData[i,nas[i,]] <- means[i]
	}
	return(methData)
}

#######################################################################################
#' convert.string.ages
#'
#' This function converts different formats of age information into numeric format to
#' be used in the following functions.
#'
#' @param age_information	input ages
#'
#' @return	the input ages in numeric format
#'
#' @author	Michael Scherer
#' @noRd

convert.string.ages <- function(age_information){
	temp <- strsplit(age_information,"[.]")
	temp <- temp[[1]][1]
	temp <- unique(na.omit(as.numeric(unlist(strsplit(temp,"[^0-9]+")))))
	if(length(temp)>1){
		temp <- temp[1]
	}
	temp
}

