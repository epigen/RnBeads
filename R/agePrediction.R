#######################################################################################
##		agePrediction.R							     ##
##	This file contains all the relevant function that are needed for the 	     ##
##	age prediction section of RnBeads.					     ##
#######################################################################################

#######################################################################################################################

#' rnb.section.ageprediction
#'
#' Adds section for the age prediction part to the report. In this function we check the available opions 
#' \bold{\code{inference.age.predidction}},\bold{\code{inference.age.prediction.training}} and \bold{\code{inference.age.prediction.predictor}}
#' and accordingly either the predefined predictor (in extdata), a newly trained predictor (in extdata) or a
#' formerly trained predictor that is stored in the path given by \bold{\code{inference.age.prediction.predictor}}.
#'
#' @param report           Report in which the age prediction section should be included. This must be an object of type \code{\linkS4class{Report}}.
#' @param object           Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @return The (possibly modified) report.
#'
#' @author Michael Scherer

rnb.section.ageprediction <- function(object,report){
	title <- "Age Prediction"
	descr <- "Plots for the visualization of predicted ages based on the DNA methylation pattern"
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
			infoText <- "Until now there is no predefined predictor available for sequencing based data. You may want to train your own age predictor based on your dataset by setting the option <em>inference.age.prediction.training</em>."
			rnb.add.paragraph(report,infoText)
			return(report)
		}
	}
	if(!is.null(predictedAges) && sum(is.na(predictedAges))==length(samples(object))){
		txt <- "The Age prediction was not successful, therefore no age prediction Section was added to the report"
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
					rnb.add.paragraph(report,txt)
			}
			header <- "The predefined predictor was used"
			prediction_path <- system.file(file.path("extdata", "predefined_predictor.csv"), package="RnBeads")
		}else{
			if(!rnb.getOption("inference.age.prediction.training")){
				header <- "A default predictor that was used"
			}
		}
	}else if(inherits(object,"RnBiseqSet") && rnb.getOption("inference.age.prediction.biseq")){
		header <- "A new predictor was trained for the specified data set"
		if(is.null(prediction_path) || prediction_path==""){
			if(rnb.getOption("inference.age.prediction.training")){
					txt <- "No age information was available for the column age in the annotation, therefore the predefined predictor is used in the further calculation."
					rnb.add.paragraph(report,txt)
			}
			txt <- "Until now there is no predefined predictor for sequencing based data available."
			rnb.add.paragraph(report,txt)
			prediction_path <- ""
		}else{
			if(!rnb.getOption("inference.age.prediction.training")){
				header <- "A default predictor that was used"
			}
		}
	}
	txt <- c(header, " and can be found as a <a href=\"",prediction_path,"\">","comma-separated file</a>.")
	rnb.add.paragraph(report, txt)

	add.info <- function(stitle, ffunction, txt, actualAges, predictedAges) {
		result <- rnb.add.section(report, stitle, txt, level = 2)
		result <- ffunction(report, object,actualAges,predictedAges)
		rnb.status(c("Added", stitle))
		result
	}

	ph <- pheno(object)
	actualAges <- ph$age
	if(!is.null(actualAges)){
		options(warn=-1)
		actualAgesNumeric <- as.numeric(actualAges)
		if(sum(is.na(actualAgesNumeric)) != length(actualAges)){
			predictedAges <- ph$predicted_ages
			txt <- "Plotting annotated ages versus predicted ages and indicating different traits with different colors."
			report <- add.info("Comparison Plot",add.agecomparison.plot,txt,actualAgesNumeric,predictedAges)
			txt <- "Plotting differences between predicted ages for each sample."
			report <- add.info("Error Plot",add.combination.plot,txt,actualAgesNumeric,predictedAges)
			#txt <- "Plotting quantiles for the difference between predicted and annotated ages."
			#report <- add.info("Quantile Plot",add.quantile.plot,txt,actualAgesNumeric,predictedAges)
		}else{
			if(is.character(actualAges)){
				fun <- function(s){
					temp <- strsplit(s,"[.]")
					temp <- temp[[1]][1]
 					unique(na.omit(as.numeric(unlist(strsplit(temp,"[^0-9]+")))))
				}
				actualAges <- unlist(lapply(actualAges,fun))
				predictedAges <- ph$predicted_ages
				txt <- "Plotting annotated ages versus predicted ages and indicating different tissues with different colors."
				report <- add.info("Comparison Plot",add.agecomparison.plot,txt,actualAges,predictedAges)
				txt <- "Plotting differences between predicted ages for each sample."
				report <- add.info("Error Plot",add.error.plot,txt,actualAges,predictedAges)
			}else{
				txt <- "The age annotation column of the dataset has a format that can not be read."
				report <- rnb.add.paragraph(report,txt)
			}
		}
	}else{
		txt <- "There was no age annotation available for the data set, therefore no inference based on age prediction could be done."
		report <- rnb.add.paragraph(report,txt)
	}
	options(warn=1)
	return(report)
}

#######################################################################################################################

#' rnb.step.ageprediction
#'
#' Creates section for the age prediction part in the report
#'
#' @param object a \code{\linkS4class{RnBeadSet}} object
#' @param report Report in which the age prediction section should be created. This must be an object of
#'               type \code{\linkS4class{Report}}.
#'
#' @return the modified report object
#'
#' @author Michael Scherer

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

	return(report)
}

#######################################################################################################################

#' rnb.execute.age.prediction
#'
#' Adds a column called 'predicted_ages' to the phenotypic table, in which the predicted ages of
#' the corresponding sample reside
#'
#' @param object a \code{\linkS4class{RnBeadSet}} object for which age prediction should
#' 	be performed
#'
#' @return the modified  \code{\linkS4class{RnBeadSet}} object
#'
#' @author Michael Scherer

rnb.execute.age.prediction <- function(object){
	if(!inherits(object,"RnBSet")){
		stop("Supplied object is not of the class inheriting from RnBSet")
	}
	prediction_path <- rnb.getOption("inference.age.prediction.predictor")
	if(inherits(object,"RnBeadSet")){
		if(!is.null(prediction_path) && prediction_path != ""){
			logger.start("Performing Age Prediction on specified predictor")
			rnbSet <- agePredictor(object,prediction_path)
			logger.completed()
		}else{
			logger.start("Performing Age Prediction")
			rnbSet <- agePredictor(object)
			logger.completed()
		}
	}else if(inherits(object,"RnBiseqSet") && rnb.getOption("inference.age.prediction.biseq")){
		if(!is.null(prediction_path) && prediction_path != ""){
			logger.start("Performing Age Prediction")
			rnbSet <- agePredictor(object,prediction_path)
			logger.completed()
		}else{
			logger.error("Currently no predefined RRBS predictor available.")
		}
	}
	return(rnbSet)
}

#######################################################################################################################

#' rnb.execute.age.training
#'
#' Trains a new predictor on the specified data set and writes it to the given path
#'
#' @param object a \code{\linkS4class{RnBeadSet}} object for which the predictor should
#' 	be trained
#' @param path path to which the predictor should be written
#'
#' @return the modified  \code{\linkS4class{RnBeadSet}} object
#'
#' @author Michael Scherer

rnb.execute.training <- function(object,path=""){
	if(path==""){
		stop("A path to store the predictor has to be specified.")
	}else{
		logger.start("Training new age predictor")
		prediction_path <- trainPredictor(object,path)
		logger.completed()
	}
}

#######################################################################################################################

#' add.agecomparison.plot
#'
#' Adds a plot that compares predicted ages by the age predictor with the annotated ages of the samples. One can compare
#' different traits with each other. Until now the following traits are (if annotated in the sample annotation sheet)
#' supported:	tissue, gender, hivstatus, sex, disease state, smoking status, ethnicity
#'
#' @param report	Report in which the corresponding plots should be integrated. 
#' @param object	Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param actualAges	Annotated ages as found in the sample annotation sheet
#' @param predictedAges	Ages as predicted by the age prediction algorithm
#' @return	The modified report with the age plots integrated
#'
#' @author	Michael Scherer

add.agecomparison.plot <- function(report, object, actualAges, predictedAges){
	
	descr <- "Comparing predicted ages by age predictor with annotated ages. Points that are labeled with their identifiers have a larger difference than 15 years between predicted and annotated age."

	ph <- pheno(object)
	Default <- rep(NA,dim(ph)[1])
	ph <- cbind(ph,Default=Default)
	copyActual <- actualAges

	traits <- names(rnb.sample.groups(object))
	traits <- c("Default",traits)

	report.plots <- list()
	if(length(traits)<=1){
		print("only one or less trait")
		actualAges <- as.numeric(actualAges)
		naAges <- is.na(actualAges)
		actualAges <- actualAges[!naAges]
		predictedAges <- predictedAges[!naAges]
		notPredicted <- is.na(predictedAges)
		actualAges <- actualAges[!notPredicted]
		predictedAges <- predictedAges[!notPredicted]
		report.plot <- createReportPlot("age_comparison_",report)
		data <- data.frame(actualAges,predictedAges,Sample=row.names(ph))
		diff <- abs(predictedAges-actualAges)
		plot <- ggplot(data,aes(x=actualAges,y=predictedAges),environment=environment())+geom_point()+geom_text(aes(label=ifelse(diff > 15, as.character(Sample),"")),hjust=0,vjust=0,size=3)+xlim(1,100)+ylim(1,100)+geom_abline(intercept=0,slope=1)+xlab("Annotated Ages")+ylab("Predicted Age")+theme(legend.title=element_blank())
		print(plot)
		report.plot <- off(report.plot)
		report.plots <- c(report.plots,report.plot)
	}else{
	for(trait in traits){
		if(trait %in% colnames(ph)){
			placeholder <- ph[,trait]
			actualAges <- copyActual
			if(!is.null(placeholder)){
				for(secondTrait in traits){
					if(secondTrait %in% colnames(ph)){
					placeholder2 <- ph[,secondTrait]
					actualAges <- copyActual
					if(!is.null(placeholder2)){
						trait <- gsub(" ","",trait)
						trait <- gsub("[[:punct:]]","",trait)
						secondTrait <- gsub(" ","",secondTrait)
						secondTrait <- gsub("[[:punct:]]","",secondTrait)
						actualAges <- as.numeric(actualAges)
						naAges <- is.na(actualAges)
						actualAges <- actualAges[!naAges]
						predictedAges <- predictedAges[!naAges]
						notPredicted <- is.na(predictedAges)
						actualAges <- actualAges[!notPredicted]
						predictedAges <- predictedAges[!notPredicted]
						placeholder <- placeholder[!naAges]
						placeholder <- placeholder[!notPredicted]
						placeholder2 <- placeholder2[!naAges]
						placeholder2 <- placeholder2[!notPredicted]
						names <- row.names(ph)
						names <- names[!naAges]
						names <- names[!notPredicted]
						report.plot <- createReportPlot(paste0("age_comparison_",trait,"_",secondTrait), report)
						if(length(placeholder)==0){
							placeholder <- rep(NA,length(actualAges))
						}
						if(length(placeholder2)==0){
							placeholder2 <- rep(NA,length(actualAges))
						}
						data <- data.frame(actualAges,predictedAges,placeholder,Sample=names)
						diff <- abs(predictedAges-actualAges)
						cvalues <- rep(rnb.getOption("colors.category"),length.out=length(placeholder))
						ptvalues <- rep(rnb.getOption("points.category"),length.out=length(placeholder2))
						plot <- ggplot(data,aes(x=actualAges,y=predictedAges,group=placeholder),environment=environment())+geom_point(aes(colour=placeholder,shape=placeholder2))+geom_text(aes(label=ifelse(diff > 15, as.character(Sample),"")),hjust=0,vjust=0,size=3)+xlim(1,100)+ylim(1,100)+geom_abline(intercept=0,slope=1)+xlab("Annotated Ages")+ylab("Predicted Age")+scale_colour_manual(name=trait,na.value = ifelse(trait=="Default","black","#C0C0C0"), values = cvalues)+scale_shape_manual(name=secondTrait,na.value=ifelse(trait=="Default",20,1L),values=ptvalues)
						print(plot)
						report.plot <- off(report.plot)
						report.plots <- c(report.plots,report.plot)
						}else{
							txt <- "There was no supported trait annotation avaible, therefore no tissue comparison for age prediction was performed. Supported traits (column names in the sample annotation sheet are: tissue, gender, hivstatus, sex, disease state, smoking status, ethnicity)"
							report <- rnb.add.paragraph(report, txt)
						}
					}
					}
				}	
			}
		}
	}
	traits <- gsub(" ","",traits)
	traits <- gsub("[[:punct:]]","",traits)
	setting_names <- list("Compare first trait"=traits,"Compare second trait"=traits[])
	names(setting_names[[1]]) <- traits
	names(setting_names[[2]]) <- traits
	report <- rnb.add.figure(report, descr, report.plots,setting_names)

	return(report)
}


#######################################################################################################################

#' add.error.plot
#'
#' Adds a plot to the given report, that describes the differences between the predicted ages by the age prediction
#' algorithm and the ages as annotated in the sample annotation sheet(if available). Possible outliers are marked by
#' their sample names.
#'
#' @param report	Report in which the corresponding plot should be integrated. 
#' @param object	Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param actualAges	Annotated ages as found in the sample annotation sheet
#' @param predictedAges	Ages as predicted by the age prediction algorithm
#' @return	The modified report with the age plot integrated
#'
#' @author	Michael Scherer


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
#' Adds a plot to the fiven report object, that shows the distribution of the differences between
#' predicted and annotated ages. Point inside the 1- and 99-quantile, respectively are marked 
#' by their sample identifiers.
#'
#' @param report	Report in which the corresponding plot should be integrated. 
#' @param object	Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param actualAges	Annotated ages as found in the sample annotation sheet
#' @param predictedAges	Ages as predicted by the age prediction algorithm
#' @return	The modified report with the quantile plot integrated
#'
#' @author	Michael Scherer

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

	plot <- ggplot(data,aes(x=Difference,y=0.0))+geom_density(aes(y=..density..))+geom_point(aes(label=Sample))+geom_text(aes(label=Sample),hjust=.5,vjust=1.5,size=3,colour='#BE1E2D')+geom_area(data=subset(density.frame,Difference>=q1 & Difference<=q99),aes(x=Difference,y=Density),fill='#75BC1C')+geom_area(data=subset(density.frame,Difference<=q1),aes(x=Difference,y=Density),fill='#BE1E2D')+geom_area(data=subset(density.frame,Difference>=q99),aes(x=Difference,y=Density),fill='#BE1E2D')+geom_vline(xintercept=median)+ylab("Density")
	report.plot <- createReportPlot("age_prediction_quantile", report)
	print(plot)
	report.plot <- off(report.plot)

	report <- rnb.add.figure(report,descr,report.plot)

	return(report)	
}

#######################################################################################################################

#' add.combination.plot
#'
#' This function adds a plot that is the combination of the upper two plots in one plot.The upper 
#' panel shows the distribution and the quantiles for the differences between predicted and annotated
#' ages and the lowe panel shows the corresponding errors
#'
#' @param report	Report in which the corresponding plot should be integrated. 
#' @param object	Methylation dataset as an object of type \code{\linkS4class{RnBeadSet}}.
#' @param actualAges	Annotated ages as found in the sample annotation sheet
#' @param predictedAges	Ages as predicted by the age prediction algorithm
#' @return	The modified report with the quantile plot integrated
#'
#' @author	Michael Scherer

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
	
	diff <- predictedAges - actualAges
	q1 <- quantile(diff,0.01,na.rm=TRUE)
	q99 <- quantile(diff,0.99,na.rm=TRUE)
	density <- density(diff)
	Sample <- rep(NA,length(density$x))
	Set <- rep("Density",length(density$y))
	density.frame <- data.frame(Sample,Difference=as.numeric(density$x),Density=as.numeric(density$y),Set)
	median <- median(diff)

	mean <- mean(diff,na.rm=TRUE)
	data <- data.frame(row.names(ph),as.numeric(diff))
	data <- cbind(data,as.numeric(row.names(data)))
	Set <- rep("Difference",length(diff))
	data <- cbind(data,Set)	
	colnames(data) <- c("Sample","Difference","Density","Set")
	hline_mean <- data.frame(yint=mean,Measure="Mean")
	hline_quantiles <- data.frame(yint=c(q1,q99),Measure="1- and 99-quantiles")

	complete_data <- rbind(density.frame,data)
	cvalues <- rep(rnb.getOption("colors.category"))

	plot <- ggplot(complete_data,aes(x=Difference,y=Density,colour=Set))+geom_point()+facet_grid(Set~.,scale="free",drop=TRUE)+geom_text(aes(label=ifelse(Difference <= q1, as.character(Sample),""),hjust=.5),vjust=0,size=3,colour="black")+geom_text(aes(x=Difference,y=Density,label=ifelse(Difference >= q99, as.character(Sample),"")),hjust=.5,vjust=0,size=3,colour="black")+geom_vline(data=hline_mean,aes(xintercept=yint,linetype=Measure),show_guide=TRUE)+geom_vline(data=hline_quantiles,aes(xintercept=yint,linetype=Measure),show_guide=TRUE)+ylab("Sample Number / Density")+xlab("Difference between predicted age and annotated age")+scale_colour_manual(na.value = "#C0C0C0", values = cvalues, guide=FALSE)
	report.plot <- createReportPlot("age_prediction_quantile", report)
	print(plot)
	report.plot <- off(report.plot)

	report <- rnb.add.figure(report,descr,report.plot)

	return(report)	
}

#######################################################################################
#' age.transformation
#'
#' This function is used to transform the inpute ages, to account for the fact, that
#' the human DNA methylation pattern goes through vast changes in the beginning of human
#' life and milder changes later on.
#'
#' The function is directly taken from Steve Horvath's prediction algorithm
#' DNA methylation age of human tissues and cell types, Steve Horvath, Genome Biology, 2013
#'
#' @param x		Input age about to be transformed
#' @param adult.age	Threshold for deciding between 'fast' and 'slow' changes in DNA methylation
#' 
#' @return	Transformed age
age.transformation <- function(x,adult.age=20)
{
	x=(x+1)/(1+adult.age);
	y=ifelse(x<=1, log( x),x-1);
	y
}

#######################################################################################
#' age.anti.transformation
#'
#' This function is used to transform the values back to actual ages, that are corresponding
#' to the predicted ages by the prediction algorithm
#'
#' The function is directly taken from Steve Horvath's prediction algorithm
#' DNA methylation age of human tissues and cell types, Steve Horvath, Genome Biology, 2013
#'
#' @param x		Input age about to be transformed back
#' @param adult.age	Threshold for deciding between 'fast' and 'slow' changes in DNA methylation
#' 
#' @return	Backtransformed age
age.anti.transformation <- function(x,adult.age=20)
{
	ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age)
}

#######################################################################################
#' agePredictor
#'
#' This function is the core of the age prediction algorithm. It takes a rnb.set object
#' and loads a given predictor from a csv file to directly predict the age of the samples
#' by this predictor.
#'
#' @param rnbSet	An \code{RnBSet} object containing the relevant methylation infos.
#' @param path		Path to a csv file in which a trained predictor should lay.
#'			DEFAULT:	The predefined predictor is loaded.
#' 
#' @return	RnBeadSet with a new column in the phenotypic table in which the predicted ages
#'		are annotated

agePredictor <- function(rnbSet, path=system.file(file.path("extdata", "predefined_predictor.csv"), package="RnBeads")){
	if(!file.exists(path)){
		logger.info(c("The specified predictor does not exists at ",path))
		return(rnbSet)
	}
	if(inherits(rnbSet,"RnBeadSet")){
		rnbSet <- agePredictor450(rnbSet,path)
	}else if(inherits(rnbSet,"RnBiseqSet")){
		rnbSet <- agePredictorRRBS(rnbSet,path)
	}
	return(rnbSet)
}

#######################################################################################
#' agePredictor450
#'
#' This function is the corresponding instance of the general age prediction algorithm for data
#' created by Illumina Infinium Bead Chips 
#'
#' @param rnbSet	An \code{RnBeadSet} object containing the relevant methylation infos.
#' @param path		Path to a csv file in which a trained predictor should lay.
#'			DEFAULT:	The predefined predictor is loaded.
#' 
#' @return	RnBeadSet with a new column in the phenotypic table in which the predicted ages
#'		are annotated

agePredictor450 <- function(rnbSet, path){
	require(impute)
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
	dummy <- capture.output(selected <- (impute.knn(selected))$data)
	selected <- t(selected)
	if(length(usedCoeffs) > dim(selected)[2]){
		usedCoeffs <- usedCoeffs[existingCpGs]
	}
	predictedAges <- intercept + selected%*%usedCoeffs
	predictedAges <- age.anti.transformation(predictedAges)
	predictedAges <- as.numeric(predictedAges)
	ages <- ph$predicted_ages
	if(is.null(ages)){
		rnbSet <- addPheno(rnbSet,predictedAges,"predicted_ages")
		actualAges <- ph$age
		if(!is.null(actualAges)){
			if(is.character(actualAges)){
				fun <- function(s){
					temp <- strsplit(s,"[.]")
					temp <- temp[[1]][1]
 					unique(na.omit(as.numeric(unlist(strsplit(temp,"[^0-9]+")))))
				}
				actualAges <- unlist(lapply(actualAges,fun))
			}
			actualAges <- as.numeric(actualAges)
			difference <- predictedAges-actualAges
			rnbSet <- addPheno(rnbSet,difference,"age_increase")
		}
	}
	options(warn=1)
	return(rnbSet)
}

#######################################################################################
#' agePredictorRRBS
#'
#' This function is the corresponding instance of the general age prediction algorithm for data
#' created by Reduced Representation Bisulphite Sequencing. WGBS is not supported until now. 
#'
#' @param rnbSet	An \code{RnBeadSet} object containing the relevant methylation infos.
#' @param path		Path to a csv file in which a trained predictor should lay.
#'			DEFAULT:	The predefined predictor is loaded.
#' 
#' @return	RnBeadSet with a new column in the phenotypic table in which the predicted ages
#'		are annotated

agePredictorRRBS <- function(rnbSet, path){
	require(impute)
	options(warn=-1)
	coeffs <- read.csv(path)
	methData <- meth(rnbSet)
	anno <- annotation(rnbSet)
	ph <- pheno(rnbSet)
	ages <- ph$age
	ages <- age.transformation(ages)
	missingAges <- is.na(ages)
	methData <- methData[,!missingAges]
	ages <- ages[!missingAges]
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
	dummy <- capture.output(selected <- (impute.knn(selected))$data)
	selected <- t(selected)
	if(dim(selected)[1]==0){
		logger.info("Age prediction could not be performed; NAs introduced")
		predictedAges <- rep(NA,length(samples(rnbSet)))
	}else{
		if(length(usedCoeffs) > dim(selected)[2]){
			usedCoeffs <- usedCoeffs[existingCpGs]
		}
		predictedAges <- intercept + selected%*%usedCoeffs
		predictedAges <- age.anti.transformation(predictedAges)
		predictedAges <- as.numeric(predictedAges)
	}
	ages <- ph$predicted_ages
	if(is.null(ages)){
		rnbSet <- addPheno(rnbSet,predictedAges,"predicted_ages")
		actualAges <- ph$age
		if(!is.null(actualAges)){
			if(is.character(actualAges)){
				fun <- function(s){
					temp <- strsplit(s,"[.]")
					temp <- temp[[1]][1]
 					unique(na.omit(as.numeric(unlist(strsplit(temp,"[^0-9]+")))))
				}
				actualAges <- unlist(lapply(actualAges,fun))
			}
			actualAges <- as.numeric(actualAges)
			difference <- predictedAges-actualAges
			rnbSet <- addPheno(rnbSet,difference,"age_increase")
		}
	}
	options(warn=1)
	return(rnbSet)
}

#######################################################################################
#' trainPredictor
#'
#' This function is the starting point for training a new predictor based on the dataset
#'
#' @param rnbSet	An \code{RnBeadSet} object containing the methylation info on which 
#'			the new predictor should be trained
#' 
#' @return	Absolute path to the corresponding predictor. The path is also set as the 
#'		option \code{inference.age.prediction.predictor}
trainPredictor <- function(rnbSet,data.dir){
	if(!file.exists(data.dir)){
		stop("The specified directory does not exist!")
	}
	if(inherits(rnbSet,"RnBeadSet")){
		path <- simpleGlmnet(rnbSet,file.path(data.dir,"trained_predictor.csv"))
	}else if(inherits(rnbSet,"RnBiseqSet")){
		path <- simpleGlmnetRRBS(rnbSet,file.path(data.dir,"trained_predictor_RRBS.csv"))
		rnb.options(inference.age.prediction.biseq=TRUE)	
	}
	rnb.options(inference.age.prediction.predictor=path)
	return(path)
}

#######################################################################################
#' simpleGlmnet
#'
#' This function actually fits the metylation data to the available age data by fitting
#' a general regularized linear model (alpha parameter = 0.8) and then building a linear
#' regression model upon the coefficients that had a non-zero coefficient in the glm. Specifically
#' dedicated to BeadChip data.
#'
#' @param trainRnBSet	An \code{RnBeadSet} object containing the methylation info on which 
#'			the new predictor should be trained
#' @param path		Path in which the new predictor should be written
#' 
#' @return	Absolute path to the corresponding predictor. Null if the function was 
#'		unable to create such an predictor.

simpleGlmnet <- function(trainRnBSet,filePath=""){
	require(glmnet)
	require(impute)
	methData <- meth(trainRnBSet)
	anno <- annotation(trainRnBSet)
	ph <- pheno(trainRnBSet)
	ages <- ph$age
	if(!is.null(ages)){
		if(is.character(ages)){
			fun <- function(s){
					temp <- strsplit(s,"[.]")
					temp <- temp[[1]][1]
 					unique(na.omit(as.numeric(unlist(strsplit(temp,"[^0-9]+")))))
				}
			ages <- unlist(lapply(ages,fun))
		}
		ages <- as.numeric(ages)
		ages <- age.transformation(ages)
		row.names(methData) <- row.names(anno)
		Xchrom <- is.element(anno$Chromosome,"chrX")
		methData <- methData[!Xchrom,]
		anno <- anno[!Xchrom,]
		Ychrom <- is.element(anno$Chromosome,"chrY")
		methData <- methData[!Ychrom,]
		anno <- anno[!Ychrom,]
		#Imputing of missing values is performed by a k-nearest-neighbors approach
		options(warn=-1)
		dummy <- capture.output(methData <- (impute.knn(methData))$data)
		options(warn=1)
		methData <- t(methData)
		missingAges <- is.na(ages)
		methData <- methData[!missingAges,]
		ages <- ages[!missingAges]
		cv <- cv.glmnet(methData,ages,parallel=TRUE)
		model <- glmnet(methData,ages,alpha=0.8,lambda=cv$lambda.min)
		coeffs <- as.matrix(coef(model))
		coeffs <- coeffs[-1,]
		non_zero <- which(coeffs!=0)
		coeffs <- coeffs[non_zero]
		max <- max(abs(coeffs))
		scaled <- abs(coeffs)/max
		sorted <- sort(scaled,decreasing=TRUE)
		cpGs <- names(sorted)
		selected <- methData[,cpGs]
		selected <- as.data.frame(selected)
		selectedPredictors <- methData[,cpGs]
		selectedPredictors <- as.data.frame(selectedPredictors)
		linearModel <- lm(ages~.,selectedPredictors)
		if(writePredictorToCsv(linearModel,filePath)){
			return(filePath)
		}
	}
	return(NULL)
}

#######################################################################################
#' simpleGlmnet
#'
#' This function actually fits the metylation data to the available age data by fitting
#' a general regularized linear model (alpha parameter = 0.8) and then building a linear
#' regression model upon the coefficients that had a non-zero coefficient in the glm.
#' Dedicated to sequencing data.
#'
#' @param trainRnBSet	An \code{RnBeadSet} object containing the methylation info on which 
#'			the new predictor should be trained
#' @param path		Path in which the new predictor should be written
#' 
#' @return	Absolute path to the corresponding predictor. Null if the function was 
#'		unable to create such an predictor.

simpleGlmnetRRBS <- function(trainRnBSet,filePath=""){
	require(glmnet)
	methData <- meth(trainRnBSet)
	anno <- annotation(trainRnBSet)
	ph <- pheno(trainRnBSet)
	ages <- ph$age
	if(!is.null(ages)){
		if(is.character(ages)){
			fun <- function(s){
					temp <- strsplit(s,"[.]")
					temp <- temp[[1]][1]
 					unique(na.omit(as.numeric(unlist(strsplit(temp,"[^0-9]+")))))
				}
			ages <- unlist(lapply(ages,fun))
		}
		ages <- as.numeric(ages)
		ages <- age.transformation(ages)
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
		methData <- imputeRRBS(methData)
		options(warn=1)
		methData <- t(methData)
		missingAges <- is.na(ages)
		methData <- methData[!missingAges,]
		ages <- ages[!missingAges]
		cv <- cv.glmnet(methData,ages)
		model <- glmnet(methData,ages,alpha=0.8,lambda=cv$lambda.min)
		coeffs <- as.matrix(coef(model))
		coeffs <- coeffs[-1,]
		non_zero <- which(coeffs!=0)
		coeffs <- coeffs[non_zero]
		max <- max(abs(coeffs))
		scaled <- abs(coeffs)/max
		sorted <- sort(scaled,decreasing=TRUE)
		cpGs <- names(sorted)
		selected <- methData[,cpGs]
		selected <- as.data.frame(selected)
		selectedPredictors <- methData[,cpGs]
		selectedPredictors <- as.data.frame(selectedPredictors)
		linearModel <- lm(ages~.,selectedPredictors)
		if(writePredictorToCsv(linearModel,filePath)){
			return(filePath)
		}
	}
	return(NULL)
}

#######################################################################################
#' writePredictorToCsv
#'
#' This function writes the coefficients of a linear model into a csv file.
#'
#' @param linearModel	Linear Model which should be stored as a csv file
#' @param path		Path in which the new predictor should be written
#' 
#' @return	TRUE if the writing is successful
writePredictorToCsv <- function(linearModel,path){
	coeffs <- as.matrix(coef(linearModel))
	write.csv(coeffs,path)
	return(TRUE)
}

mapToCpG <- function(methData,anno){
	anno <- GRanges(Rle(anno$Chromosome),IRanges(start=anno$Start,end=anno$End),anno$Strand)
	anno_chip <- rnb.get.annotation("probes450")
	anno_chip <- unlist(anno_chip)
	names <- names(anno_chip)
	names <- strsplit(names,"[[:punct:]]")
	names <- lapply(names,function(x){x[2]})
	names <- unlist(names)
	overlap <- findOverlaps(anno_chip,anno)
	hits_anno <- subjectHits(overlap)
	hits_chip <- queryHits(overlap)
	methData <- methData[hits_anno,]
	names <- names[hits_chip]
	row.names(methData) <- names
	return(methData)
}

#######################################################################################
#' imputeRRBS
#'
#' Since the k-nearest-neighbroing approach does not converge for all CpGs present in
#' sequencing data, a simple mean imputation is performed for this kind of data
#'
#' @param methData	Methylation data as a data frame or matix to be imputed
#' 
#' @return	the changed methyaltion data frame

imputeRRBS <- function(methData){
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
