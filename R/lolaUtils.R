########################################################################################################################
## lolaUtils.R
## created: 2017-07-06
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for LOLA enrichment analysis
########################################################################################################################

#' getCellTypesFromLolaDb
#'
#' retrieve or guess cell types from a LOLA DB object
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @return character vector with cell types
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getCellTypesFromLolaDb(lolaDb)
#' }
getCellTypesFromLolaDb <- function(lolaDb){
	tt <- data.frame(lolaDb$regionAnno)
	res <- tt$cellType

	# special treatment for "sheffield_dnase" from LOLACore
	isInCollection <- tt$collection=="sheffield_dnase"
	if (sum(isInCollection) > 0) res[isInCollection] <- tt$description[isInCollection]

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub("^(.+)-(.+)$", "\\1", tt$filename[isInCollection])
	}

	return(res)
}

#' getTargetFromLolaDb
#'
#' retrieve or guess the target from a LOLA DB object. Here, target typically
#' refers to antibodies for ChIP-seq experiments, but could also refer to other annotations
#' (e.g. motifs in TF motif databases, annotation according to UCSC features etc.)
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @return character vector with targets
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getTargetFromLolaDb(lolaDb)
#' }
getTargetFromLolaDb <- function(lolaDb){
	tt <- data.frame(lolaDb$regionAnno)
	res <- tt$antibody

	# special treatment for "ucsc_features" from LOLACore
	isInCollection <- tt$collection=="ucsc_features"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("^UCSC ", "", tt$description[isInCollection])

	# special treatment for "sheffield_dnase" from LOLACore
	isInCollection <- tt$collection=="sheffield_dnase"
	if (sum(isInCollection) > 0) res[isInCollection] <- "DNase"

	# special treatment for "encode_segmentation" from LOLACore
	isInCollection <- tt$collection=="encode_segmentation"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub(" Segments$", "", tt$description[isInCollection])
	}

	# special treatment for "jaspar_motifs" from LOLAExt
	isInCollection <- tt$collection=="jaspar_motifs"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("\\.bed$", "", tt$filename[isInCollection])

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub("(\\.(hotspot\\.(fdr0\\.01|all)|macs2))?\\.(peaks\\.bed|narrowPeak)$", "", 
			gsub("^(.+)-(.+)$", "\\2", tt$filename[isInCollection])
		)
	}

	return(res)
}

#' getNamesFromLolaDb
#'
#' get human readable names from a LOLA DB object
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param addCollectionNames  attach the name of the collection to the name
#' @param addDbId  attach the index of the item in the LOLA DB object to the name
#' @return character vector with human readable names
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getNamesFromLolaDb(lolaDb)
#' }
getNamesFromLolaDb <- function(lolaDb, addCollectionNames=FALSE, addDbId=TRUE){
	tt <- data.frame(lolaDb$regionAnno)
	res <- paste0(tt$description)

	cellTypes <- getCellTypesFromLolaDb(lolaDb)
	targets <- getTargetFromLolaDb(lolaDb)

	# if there is a valid "target" (see getTargetFromLolaDb) inferred from the annotation, construct a name from it
	hasDesc <- !is.na(targets)
	if (sum(hasDesc) > 0) res[hasDesc] <- targets[hasDesc]

	# if both cellType and target annotation were inferred, construct a name from them
	hasDesc <- !is.na(cellTypes) & !is.na(targets)
	if (sum(hasDesc) > 0) res[hasDesc] <- paste0(targets[hasDesc], " - ", cellTypes[hasDesc])

	# add suffixes if requested
	if (addCollectionNames){
		res <- paste0(res, " [", tt$collection, "]")
	}
	if (addDbId){
		res <- paste0(res, " [db", 1:nrow(tt), "]")
	}
	return(res)
}

#' lolaPrepareDataFrameForPlot
#'
#' Helper function for preparing a data.frame for multiple LOLA plots
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}.
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut p-value cutoff to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @noRd
lolaPrepareDataFrameForPlot <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, perUserSet=FALSE, groupByCollection=TRUE, orderDecreasing=NULL){
	#dedect by column name whether decreasing order needs to be used
	if (is.null(orderDecreasing)){
		oset.dec <- c("pValueLog", "qValueLog", "logOddsRatio", "oddsRatio")
		oset.inc <- c("maxRnk", "meanRnk", "qValue")
		if (!is.element(orderCol, c(oset.dec, oset.inc))){
			logger.error(c("Could not determine whether to use increasing or decreasing order for column:", orderCol))
		}
		orderDecreasing <- is.element(orderCol, oset.dec)
	}

	if (!is.element(signifCol, c("pValueLog", "qValue"))){
		logger.error(c("Invalid significance column name:", signifCol))
	}

	lolaRes <- as.data.frame(lolaRes)

	if (signifCol == "qValue"){
		lolaRes[["qValueLog"]] <- -log10(lolaRes[["qValue"]])
		signifCol <- "qValueLog"
	}

	#check if column name is in the LOLA results
	# take older versions of LOLA into account (oddsRatio was named logOddsRatio before)
	if (is.element("logOddsRatio", colnames(lolaRes))){
		logger.warning("Detected old version of LOLA. Renaming 'logOddsRatio' to 'oddsRatio'")
		colnames(lolaRes)[colnames(lolaRes)=="logOddsRatio"] <- "oddsRatio"
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	if (orderCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		orderCol <- "oddsRatio"
	}

	lolaRes$name <- getNamesFromLolaDb(lolaDb, addCollectionNames=!groupByCollection)[lolaRes$dbSet]
	lolaRes$name <- trimChar(lolaRes$name, len.out=50, trim.str="...", len.pref=30)
	lolaRes$target <- getTargetFromLolaDb(lolaDb)[lolaRes$dbSet]

	if (!is.element(scoreCol, colnames(lolaRes))){
		logger.error(c("Missing score column from LOLA result:", scoreCol))
	}
	if (!is.element(orderCol, colnames(lolaRes))){
		logger.error(c("Missing order column from LOLA result:", orderCol))
	}

	if (length(includedCollections) > 0){
		lolaRes <- lolaRes[lolaRes[["collection"]] %in% includedCollections,]
	}

	calledSignif <- lolaRes[[signifCol]] > -log10(pvalCut)

	if (sum(calledSignif) < 1){
		# lolaRes.signif <- lolaRes
		logger.warning("No significant association found")
		return(NULL)
	} else {
		lolaRes.signif <- lolaRes[calledSignif,]
	}
	lolaRes.signif <- lolaRes.signif[order(lolaRes.signif[,orderCol], decreasing=orderDecreasing),]

	# data frame for plotting sorted according to significance/odds ratio
	df2p <- lolaRes.signif
	dbSets2keep <- unique(df2p[["dbSet"]])
	if (length(dbSets2keep) > maxTerms){
		if (perUserSet){
			logger.info(c("Reduced the number of region sets in the plot to", maxTerms, "per userSet"))
			uSetF <- lolaRes[["userSet"]]
			if (!is.factor(uSetF)) uSetF <- factor(uSetF)
			dbSets2keep <- c()
			for (ss in levels(uSetF)){
				dbSets.cur <- unique(df2p[uSetF==ss,"dbSet"])
				dbSets2keep <- union(dbSets2keep, dbSets.cur[1:min(maxTerms, length(dbSets.cur))])
			}
		} else {
			logger.info(c("Reduced the number of region sets in the plot to", maxTerms))
			dbSets2keep <- dbSets2keep[1:maxTerms]
		}
		df2p <- df2p[df2p[["dbSet"]] %in% dbSets2keep, ]
	}

	df2p[["term"]] <- factor(df2p[["name"]], levels=unique(df2p[["name"]])) #sort factor by occurence (data.frame is ordered)
	df2p[["target"]] <- factor(df2p[["target"]])

	#bound score at maximum (excluding +/-Inf)
	scoresFinite <- df2p[[scoreCol]][is.finite(df2p[[scoreCol]])]
	cutVal <- max(scoresFinite, na.rm=TRUE)
	df2p[[scoreCol]][df2p[[scoreCol]]>cutVal] <- cutVal
	cutVal <- min(scoresFinite, na.rm=TRUE)
	df2p[[scoreCol]][df2p[[scoreCol]]<cutVal] <- cutVal

	return(df2p)
}


#' lolaVolcanoPlot
#'
#' plot a volcano plot showing LOLA enrichment results: LOLA p-value against the log-odds score. Colored by rank
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}.
#' @param colorBy  annotation/column in the the LOLA DB that should be used for point coloring
#' @param colorpanel colors to be used for coloring the points
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaVolcanoPlot(res$lolaDb, lolaRes, signifCol="qValue")
#' }
lolaVolcanoPlot <- function(lolaDb, lolaRes, includedCollections=c(), signifCol="qValue", colorBy="maxRnk", colorpanel=c()){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}

	#prepare data.frame for plotting
	# adjust maxTerms, pvalCut to not filter anything
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol="pValueLog", signifCol=signifCol, orderCol="maxRnk", includedCollections=includedCollections, pvalCut=1.1, maxTerms=Inf, perUserSet=FALSE, groupByCollection=TRUE, orderDecreasing=NULL)
	#reverse the order of points s.t. the plot looks nice
	df2p <- df2p[nrow(df2p):1,]

	if (signifCol == "qValue"){
		signifCol <- "qValueLog"
	}

	is.color.gradient <- FALSE
	is.color.discrete <- FALSE
	if (!is.null(colorBy) && nchar(colorBy) > 0 && !is.element(colorBy, colnames(df2p))){
		logger.warning(c("Invalid colorBy argument:", colorBy, "--> using black as color"))
	} else {
		is.color.gradient <- is.numeric(df2p[[colorBy]])
		is.color.discrete <- !is.color.gradient
	}

	oddsRatioCol <- "oddsRatio"
	if (!is.element(oddsRatioCol, colnames(df2p))){
		logger.error("Invalid LOLA result. Could not find a valid column containing odds ratios.")
	}
	pp <- ggplot(df2p) + aes_string(oddsRatioCol, signifCol, color=colorBy) + geom_point()
	cpanel <- colorpanel
	if (length(cpanel) < 1){
		if (is.color.gradient){
			cpanel <- rev(rnb.getOption("colors.gradient"))
		}
		if (is.color.discrete){
			if (is.factor(df2p[[colorBy]])){
				cpanel <- rep(rnb.getOption("colors.category"), length.out=nlevels(df2p[[colorBy]]))
				names(cpanel) <- levels(df2p[[colorBy]])
			} else if (is.character(df2p[[colorBy]])){
				cpanel.names <- unique(df2p[[colorBy]])
				cpanel <- rep(rnb.getOption("colors.category"), length.out=length(cpanel.names))
				names(cpanel) <- cpanel.names
			} else {
				logger.error("invalid discrete coloring column type")
			}
		}
	}
	if (is.color.gradient){
		pp <- pp + scale_color_gradientn(colors=cpanel)
	}
	if (is.color.discrete){
		pp <- pp + scale_color_manual(na.value="#C0C0C0", values=cpanel)
	}
		
	return(pp)
}

#' lolaBarPlot
#'
#' plot a barplot of LOLA enrichment results
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut p-value cutoff to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaBarPlot(res$lolaDb, lolaRes, scoreCol="oddsRatio", orderCol="maxRnk", pvalCut=0.05)
#' }
lolaBarPlot <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, colorpanel=sample(rainbow(maxTerms,v=0.5)), groupByCollection=TRUE, orderDecreasing=NULL){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	#prepare data.frame for plotting
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, signifCol=signifCol, includedCollections=includedCollections, pvalCut=pvalCut, maxTerms=maxTerms, perUserSet=FALSE, groupByCollection=groupByCollection, orderDecreasing=orderDecreasing)

	if (is.null(df2p)){
		return(rnb.message.plot("No significant association found"))
	}

	aesObj <- aes_string("term", scoreCol)
	if (length(colorpanel) > 0) aesObj <- aes_string("term", scoreCol, fill="target")

	pp <- ggplot(df2p) + aesObj + geom_col() + 
		  scale_x_discrete(name="") + 
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

	if (length(colorpanel) > 0){
		# check if the defined color panel already contains a mapping of targets to colors
		# otherwise create an arbitrary mapping
		if (!is.null(names(colorpanel)) && all(levels(df2p[["target"]]) %in% names(colorpanel))){
			cpanel <- colorpanel
		} else {
			cpanel <- rep(colorpanel, length.out=nlevels(df2p[["target"]]))
			names(cpanel) <- levels(df2p[["target"]])
		}
		pp <- pp + scale_fill_manual(na.value="#C0C0C0", values=cpanel, guide=FALSE)
	}

	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}

#' lolaBarPlot.hyp
#'
#' Adaptation of  \code{\link{lolaBarPlot}} to plot hypo and hpermethylated regions side by side
#' Hyper-/hypomethylation is assumed to be annotated in the userSet column of lolaRes
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut p-value cutoff (negative log10) to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @noRd
lolaBarPlot.hyp <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, colorpanel=c("#a6611a", "#018571"), groupByCollection=TRUE, orderDecreasing=NULL){
	isHypo  <- grepl("hypo", lolaRes[["userSet"]], ignore.case=TRUE)
	isHyper <- grepl("hyper", lolaRes[["userSet"]], ignore.case=TRUE)
	if (sum(isHypo) < 1) logger.error("LOLA Result does not contain hypOmethylated userSet")
	if (sum(isHyper) < 1) logger.error("LOLA Result does not contain hypERmethylated userSet")
	lolaRes$userSet[isHypo]  <- "hypomethylated"
	lolaRes$userSet[isHyper] <- "hypermethylated"
	if (sum(isHypo) + sum(isHyper) != nrow(lolaRes)){
		lolaRes <- lolaRes[isHypo | isHyper,]
		logger.warning("LOLA Result contains userSets not annotated as hyper or hypo")
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	#prepare data.frame for plotting
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, signifCol=signifCol, includedCollections=includedCollections, pvalCut=pvalCut, maxTerms=maxTerms, perUserSet=TRUE, groupByCollection=groupByCollection, orderDecreasing=orderDecreasing)

	if (is.null(df2p)){
		return(rnb.message.plot("No significant association found"))
	}

	cpanel <- c(hypomethylated=colorpanel[1], hypermethylated=colorpanel[2])
	aesObj <- aes_string("term", scoreCol, fill="userSet")

	pp <- ggplot(df2p) + aesObj + geom_col(position="dodge") + 
		  scale_x_discrete(name="") + 
		  scale_fill_manual(values=cpanel) +
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}

#' lolaBoxPlotPerTarget
#'
#' plot a boxplot showing LOLA enrichment results per "target" group (see \code{\link{getTargetFromLolaDb}} for an explanation of
#' "target").
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{\link{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param signifCol column name of the significance score in \code{lolaRes}. Should be one of \code{c("pValueLog", "qValue")}
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut  p-value cutoff to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \code{\link{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @param scoreDecreasing flag indicating whether the value in \code{scoreCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaBoxPlotPerTarget(res$lolaDb, lolaRes, scoreCol="oddsRatio", orderCol="maxRnk", pvalCut=0.05)
#' }
lolaBoxPlotPerTarget <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, signifCol="qValue", includedCollections=c(), pvalCut=0.01, maxTerms=50, colorpanel=c(), groupByCollection=TRUE, orderDecreasing=NULL, scoreDecreasing=NULL){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}
	# detect by column name whether decreasing order needs to be used for the score column
	if (is.null(scoreDecreasing)){
		oset.dec <- c("pValueLog", "logOddsRatio", "oddsRatio")
		oset.inc <- c("maxRnk", "meanRnk")
		if (!is.element(scoreCol, c(oset.dec, oset.inc))){
			logger.error(c("Could not determine whether to use increasing or decreasing order for score column:", scoreCol))
		}
		scoreDecreasing <- is.element(scoreCol, oset.dec)
	}
	if (scoreCol=="logOddsRatio") {
		logger.warning("In newer versions of LOLA the odds ratio column is called 'oddsRatio' (no longer 'logOddsRatio')")
		scoreCol <- "oddsRatio"
	}
	#prepare data.frame for plotting
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, signifCol=signifCol, includedCollections=includedCollections, pvalCut=pvalCut, maxTerms=Inf, perUserSet=FALSE, groupByCollection=groupByCollection, orderDecreasing=orderDecreasing)

	if (is.null(df2p)){
		return(rnb.message.plot("No significant association found"))
	}

	# maxTerms now applies to the targets, not to LOLA DB items
	if (nlevels(df2p$target) > maxTerms){
		# if there are too many targets, select the targets from the top of the ordered table
		targets2keep <- unique(as.character(df2p$target))[1:maxTerms]
		df2p <- df2p[df2p$target %in% targets2keep,]
		df2p[["target"]] <- factor(as.character(df2p[["target"]]))
		logger.info(c("Reduced the number of targets in the plot to", maxTerms))
	}

	# sort targets by median score
	targetScore <- tapply(df2p[[scoreCol]], df2p[["target"]], FUN=function(x){median(x, na.rm=TRUE)})
	df2p[["target"]] <- factor(as.character(df2p[["target"]]),
		levels=names(targetScore)[order(targetScore, decreasing=scoreDecreasing)]
	)

	aesObj <- aes_string("target", scoreCol)
	if (length(colorpanel) > 0) aesObj <- aes_string("target", scoreCol, fill="target")

	pp <- ggplot(df2p) + aesObj + geom_boxplot(outlier.shape=NA) +
		  geom_dotplot(aes(fill=NULL), binaxis="y", stackdir="center") +
		  scale_x_discrete(name="") + 
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

	if (length(colorpanel) > 0){
		# check if the defined color panel already contains a mapping of targets to colors
		# otherwise create an arbitrary mapping
		if (!is.null(names(colorpanel)) && all(levels(df2p[["target"]]) %in% names(colorpanel))){
			cpanel <- colorpanel
		} else {
			cpanel <- rep(colorpanel, length.out=nlevels(df2p[["target"]]))
			names(cpanel) <- levels(df2p[["target"]])
		}
		pp <- pp + scale_fill_manual(na.value="#C0C0C0", values=cpanel, guide=FALSE)
	} 
	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}
