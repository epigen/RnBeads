########################################################################################################################
## filteringSummary.R
## created: 2013-12-12
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the summary section after all filtering steps have been performed.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

## rnb.get.filtered.sites.samples
##
## Validates the given datasets before and after filtering and extracts the lists of filtered sites and samples.
##
## @param old.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
## @param new.set Methylation dataset after filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
## @return \code{list} of three elements: \code{"mm"}, \code{"samples"} and \code{"sites"}.
##
## @author Yassen Assenov
rnb.get.filtered.sites.samples <- function(old.set, new.set) {

	## Validate that new.set is a subset of old.set
	mm.old <- meth(old.set, row.names = TRUE)
	mm.new <- meth(new.set, row.names = TRUE)
	sites.old <- rownames(mm.old)
	sites.new <- rownames(mm.new)
	samples.old <- colnames(mm.old)
	samples.new <- colnames(mm.new)
	rm(mm.new)
	if (length(setdiff(sites.new, sites.old)) != 0) {
		stop("inconsistent sites in old.set and new.set")
	}
	if (length(setdiff(samples.new, samples.old)) != 0) {
		stop("inconsistent samples in old.set and new.set")
	}

	removed.samples <- which(!(samples.old %in% samples.new))
	removed.sites <- which(!(sites.old %in% sites.new))
	return(list(mm = mm.old, samples = removed.samples, sites = removed.sites))
}

########################################################################################################################

#' rnb.filter.dataset
#' 
#' Modifies the given methylation dataset after a sequence of filtering steps has been performed and summarized.
#' 
#' @param rnb.set   Methylation dataset as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param r.samples Indices of all samples to be removed as a result of performing filtering steps.
#' @param r.sites   Indices of all sites/probes to be removed as a result of performing filtering steps.
#' @param mask      Optionally, masking of low quality measures to be applied. If specified, this parameter must be a
#'                  \code{matrix} of \code{logical} values with dimensions corresponding to the size of \code{rnb.set}.
#' @return The possibly modified dataset.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.filter.dataset <- function(rnb.set, r.samples, r.sites, mask = NULL) {
	logger.start("Manipulating the object")
	needs.summary <- (!is.null(mask))
	if (needs.summary) {
		rnb.set@meth.sites[,][mask] <- NA
	}
	if (length(r.samples) != 0) {
		if(rnb.getOption("enforce.destroy.disk.dumps")){
			rnb.set@status$discard.ff.matrices<-TRUE
		}
		rnb.set <- remove.samples(rnb.set, r.samples)
		if(isTRUE(rnb.set@status$discard.ff.matrices)){
			rnb.set@status$discard.ff.matrices<-NULL
		}
		# needs.summary <- FALSE
		logger.status(sprintf("Removed %d samples", length(r.samples)))
	}
	if (length(r.sites) != 0) {
		if(rnb.getOption("enforce.destroy.disk.dumps")){
			rnb.set@status$discard.ff.matrices<-TRUE
		}
		rnb.set <- remove.sites(rnb.set, r.sites)
		if(isTRUE(rnb.set@status$discard.ff.matrices)){
			rnb.set@status$discard.ff.matrices<-NULL
		}
		needs.summary <- FALSE
		logger.status(sprintf("Removed %d sites (probes)", length(r.sites)))
	}
	if (needs.summary) {
		rnb.set <- updateRegionSummaries(rnb.set)
		logger.status("Updated region-level data")
	}
	logger.completed()
	rnb.set
}

########################################################################################################################

#' rnb.execute.filter.summary
#'
#' Calculates a table summarizing the effect of the applied filtering procedures.
#'
#' @param old.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param new.set Methylation dataset after filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @return \code{matrix} summarizing the number of removed and retained sites, samples, and (optionally) reliable and
#'         unreliable measurements.
#'
#' @details
#' This function expects that the sites and samples in \code{new.set} are subsets of the sites and samples in
#' \code{old.set}, respectively. If this is not the case, it exists with an error.
#'
#' @seealso \code{\link{rnb.run.preprocessing}} for running the whole preprocessing module
#'
#' @author Yassen Assenov
#' @export
rnb.execute.filter.summary <- function(old.set, new.set) {
	if (!inherits(old.set, "RnBSet")) {
		stop("invalid value for old.set")
	}
	if (!inherits(new.set, "RnBSet")) {
		stop("invalid value for new.set")
	}

	removed <- rnb.get.filtered.sites.samples(old.set, new.set)
	rnb.execute.filter.summary.internal(old.set, removed$samples, removed$sites)
}

rnb.execute.filter.summary.internal <- function(rnb.set, removed.samples, removed.sites) {
	cont.matrix <- rbind(
		"Sites" = as.double(c(nsites(rnb.set) - length(removed.sites), length(removed.sites))),
		"Samples" = as.double(c(length(samples(rnb.set)) - length(removed.samples), length(removed.samples))))
	rownames(cont.matrix)[1] <- capitalize(rnb.get.row.token(class(rnb.set), plural = TRUE))
	colnames(cont.matrix) <- c("Retained", "Removed")
	if (rnb.has.reliability.info(rnb.set)) {
		cmatrix <- matrix(0, nrow = 2, ncol = 2)
		rownames(cmatrix) <- paste(c("Reliable", "Unreliable"), "measurements")

		hasRemovedSites <- length(removed.sites) > 0
		hasRemovedSamples <- length(removed.samples) > 0
		if (hasRemovedSites){
			relCounts.remSites <- rnb.get.reliability.counts.per.sample(rnb.set, siteIndices=removed.sites)
			notRelCounts.remSites <- length(removed.sites) - relCounts.remSites
			relCounts.retSites <- rnb.get.reliability.counts.per.sample(rnb.set, siteIndices=-removed.sites)
			notRelCounts.retSites <- nsites(rnb.set) - length(removed.sites) - relCounts.retSites
			if (hasRemovedSamples){
				cmatrix[1,1] <- sum(as.numeric(relCounts.retSites[-removed.samples]))
				cmatrix[1,2] <- sum(as.numeric(relCounts.remSites[-removed.samples])) + sum(as.numeric(relCounts.retSites[removed.samples]) + as.numeric(relCounts.remSites[removed.samples]))
				cmatrix[2,1] <- sum(as.numeric(notRelCounts.retSites[-removed.samples]))
				cmatrix[2,2] <- sum(as.numeric(notRelCounts.remSites[-removed.samples])) + sum(as.numeric(notRelCounts.retSites[removed.samples]) + as.numeric(notRelCounts.remSites[removed.samples]))
			} else {
				cmatrix[1,1] <- sum(as.numeric(relCounts.retSites))
				cmatrix[1,2] <- sum(as.numeric(relCounts.remSites))
				cmatrix[2,1] <- sum(as.numeric(notRelCounts.retSites))
				cmatrix[2,2] <- sum(as.numeric(notRelCounts.remSites))
			}
		} else {
			relCounts <- rnb.get.reliability.counts.per.sample(rnb.set, siteIndices=NULL)
			notRelCounts <- nsites(rnb.set) - relCounts
			if (hasRemovedSamples){
				cmatrix[1,1] <- sum(as.numeric(relCounts[-removed.samples]))
				cmatrix[1,2] <- sum(as.numeric(relCounts[ removed.samples]))
				cmatrix[2,1] <- sum(as.numeric(notRelCounts[-removed.samples]))
				cmatrix[2,2] <- sum(as.numeric(notRelCounts[ removed.samples]))
			} else {
				cmatrix[1,1] <- sum(as.numeric(relCounts))
				cmatrix[1,2] <- as.numeric(0)
				cmatrix[2,1] <- sum(as.numeric(notRelCounts))
				cmatrix[2,2] <- as.numeric(0)
			}
		}
		cont.matrix <- rbind(cont.matrix, cmatrix)
	}
	cbind(cont.matrix, "Total" = rowSums(cont.matrix))
}

########################################################################################################################

#' rnb.step.filter.summary
#'
#' Calculates a table summarizing the effect of the applied filtering procedures and adds a corresponding section to the
#' given report.
#'
#' @param old.set Methylation dataset before filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param new.set Methylation dataset after filtering as an object of type inheriting \code{\linkS4class{RnBSet}}.
#' @param report  Report to summarize the outcome of this procedure. This must be an object of type
#'                \code{\linkS4class{Report}}.
#' @return The report, possibly modified.
#'
#' @details
#' This function expects that the sites and samples in \code{new.set} are subsets of the sites and samples in
#' \code{old.set}, respectively. If this is not the case, it exists with an error. If \code{old.set} and \code{new.set}
#' are identical, no information is added to the given report.
#'
#' @seealso \code{\link{rnb.execute.filter.summary}}, \code{\link{rnb.run.filtering}}
#'
#' @author Yassen Assenov
#' @noRd
rnb.step.filter.summary <- function(old.set, new.set, report) {
	if (!inherits(old.set, "RnBSet")) {
		stop("invalid value for old.set")
	}
	if (!inherits(new.set, "RnBSet")) {
		stop("invalid value for new.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	if (rnb.getOption("logging") && logger.isinitialized() == FALSE) {
		logger.start(fname = NA) # initialize console logger
	}

	removed <- rnb.get.filtered.sites.samples(old.set, new.set)
	rnb.step.filter.summary.internal(old.set, removed$samples, removed$sites, report,
		TRUE)
}

rnb.step.filter.summary.internal <- function(rnb.set, removed.samples, removed.sites, report,
	log.section = FALSE, section.name="Filtering Summary", section.order=1) {

	if (log.section) {
		logger.start("Summary of Filtering Procedures")
	}

	## Create a summary table of removed sites, samples and unreliable measurements
	logger.status("Creating summary table of removed sites, samples and unreliable measurements...")
	table.summary <- rnb.execute.filter.summary.internal(rnb.set, removed.samples, removed.sites)
	if (all(table.summary[, "Removed"] == 0)) {
		if (log.section) logger.completed()
		return(report)
	}

	## Save table and create figure
	logger.status("Saving table and figures...")
	fname <- sprintf("summary%d.csv", section.order)
	utils::write.csv(table.summary, file = file.path(rnb.get.directory(report, "data", TRUE), fname))
	dframe <- table.summary[, -3]
	colnames(dframe) <- tolower(colnames(dframe))
	dframe <- data.frame(
		x = as.vector(dframe),
		item = factor(rep(rownames(dframe), ncol(dframe)), levels = rev(rownames(dframe))),
		group = factor(rep(colnames(dframe), each = nrow(dframe)), levels = colnames(dframe)))
	pp <- ggplot(dframe, aes_string("item", weight = "x", fill = "factor(group, levels = rev(levels(group)))")) +
		scale_y_continuous(limits = c(0, 1), expand = c(0, 0), labels = percent_format()) +
		labs(x = NULL, y = "Fraction", fill = "Group") + theme(axis.ticks.y = element_blank()) +
		ggplot2::geom_bar(position = "fill", width = 0.7) + coord_flip() +
		theme(plot.margin = unit(0.1 + c(0, 0, 0, 0), "in"))
	pname <- sprintf("summary%d_barchart", section.order)
	rplot <- createReportPlot(pname, report, width = 7, height = 0.55 + 0.4 * nrow(table.summary))
	suppressWarnings(print(pp))
	rplot <- off(rplot)
	fname <- paste(rnb.get.directory(report, "data"), fname, sep = "/")
	dataset.class <- class(rnb.set)
	txt.site <- rnb.get.row.token(dataset.class)
	txt.sites <- rnb.get.row.token(dataset.class, plural = TRUE)
	rems <- table.summary[c(capitalize(txt.sites), "Samples"), "Removed"]
	txt <- c("As a final outcome of the filtering procedures, ", rems[1], " ",
		ifelse(rems[1] == 1, txt.site, txt.sites), " and ", rems[2], " sample", ifelse(rems[2] != 1, "s", ""), " were ",
		"removed. These statistics are presented in <a href=\"", fname, "\">a dedicated table</a> that accompanies ",
		"this report and visualized in the figure below.")
	report <- rnb.add.section(report, section.name, txt)
	txt <- "Fractions of removed values in the dataset after applying filtering procedures."
	report <- rnb.add.figure(report, txt, rplot)
	rm(table.summary, fname, dframe, pp, pname, rplot, rems, txt)
	logger.status("Added summary table of removed and retained items")

	## Construct vectors of removed and retained betas
	siteInds <- NULL
	# subsample sites from removed and retained sites in order to reduce memory footprint 
	if (rnb.getOption("distribution.subsample") > 1){
		numSites <- rnb.getOption("distribution.subsample")
		if (length(removed.sites) != 0) {
			if (length(removed.sites) > numSites){
				siteInds <- sort(sample(removed.sites, size=numSites))
			} else {
				siteInds <- removed.sites
			}
			retained.sites <- setdiff(1:nsites(rnb.set), removed.sites)
			if (length(retained.sites) > numSites){
				siteInds <- sort(c(siteInds, sample(retained.sites, size=numSites)))
			} else {
				siteInds <- sort(c(siteInds, retained.sites))
			}
		} else {
			siteInds <- sort(sample.int(nsites(rnb.set), size=numSites))
		}
	}
	if (!is.null(siteInds)){
		logger.info(c("Subsampling", length(siteInds), "sites for plotting density distributions"))
		# new site indices in subsampled set
		indexMap <- rep(NA, nsites(rnb.set))
		indexMap[siteInds] <- 1:length(siteInds)
		removed.sites <- na.omit(indexMap[removed.sites])
	}
	
	mm <- meth(rnb.set, i=siteInds)
	if (length(removed.samples) != 0) {
		betas.removed <- as.vector(mm[, removed.samples])
		mm <- mm[, -removed.samples]
	} else {
		betas.removed <- double()
	}
	if (length(removed.sites) != 0) {
		betas.removed <- c(betas.removed, mm[removed.sites, ])
		mm <- mm[-removed.sites, ]
	}
	beta.values <- list(Removed = betas.removed, Retained = as.vector(mm))
	for (i in 1:length(beta.values)) {
		beta.values[[i]] <- beta.values[[i]][!is.na(beta.values[[i]])]
	}
	logger.status("Constructed sequences of removed and retained methylation values")
	rm(betas.removed, mm, i)
	rnb.cleanMem()

	## Compare removed vs. retained methylation beta values
	if ((!is.null(beta.values)) && min(sapply(beta.values, length)) >= 501) {
		report.plots <- rnb.plot.beta.comparison(beta.values, sprintf("summary%d_betas", section.order), report)
		setting.names <- list("Plot type" =
			c("density" = "density estimation", "histogram" = "histograms", "qq" = "quantile-quantile plot"))
		txt <- c("The figure below compares the distributions of the removed methylation &beta; values and of the ",
			"retained ones.")
		rnb.add.paragraph(report, txt)
		txt <- "Comparison of removed and retained &beta; values."
		txt <- c(txt, add.text.subsampling(attr(report.plots, "subsampled"), paste(names(beta.values), "betas")))
		report <- rnb.add.figure(report, txt, report.plots, setting.names)
		logger.status("Added comparison between removed and retained beta values")
	}

	if (log.section) {
		logger.completed()
	}
	return(report)
}

## E N D ###############################################################################################################
