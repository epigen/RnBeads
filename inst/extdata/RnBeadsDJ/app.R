################################################################################
# RnBeadsDJ - RnBeads Data Juggler
# A Shiny app for configuring and running RnBeads analyses
#-------------------------------------------------------------------------------
# Main app file
#-------------------------------------------------------------------------------
# created: 2017-10-04
# creator: Fabian Mueller
#-------------------------------------------------------------------------------
# Run using:
# shiny::runApp(file.path('RnBeadsDJ'))
################################################################################
library(shiny)
library(shinyjs)
# library(shinyFiles)
library(RnBeads)

################################################################################
# Options
################################################################################
options(shiny.maxRequestSize=50*1024^2)

################################################################################
# globals
################################################################################
REFRESH.TIME <- 5000
RNB.MODULES <- c(
	"Data Import"="data_import",
	"Quality Control"="quality_control",
	"Preprocessing"="preprocessing",
	"Tracks and Tables"="tracks_and_tables",
	"Covariate Inference"="covariate_inference",
	"Exploratory Analysis"="exploratory_analysis",
	"Differential Methylation"="differential_methylation"
)
RNB.MODULES.LOG.MSG <- c(
	"data_import"="Loading Data",
	"quality_control"="Quality Control",
	"preprocessing"="Preprocessing",
	"tracks_and_tables"="Tracks and Tables",
	"covariate_inference"="Covariate Inference",
	"exploratory_analysis"="Exploratory Analysis",
	"differential_methylation"="Differential Methylation"
)
RNB.PLATFORMS <- c("Bisulfite Sequencing"="biseq", "Illumina EPIC"="illEpic", "Illumina EPIC v2"="illEpicv2", "Illumina 450k"="ill450k", "Illumina27k"="ill27k")
RNB.ASSEMBLIES <- rnb.get.assemblies()
RNB.TABLE.SEPS <- c("comma" = ",", "tab"="\t")
RNB.BED.STYLES <- c("BisSNP"="BisSNP", "ENCODE"="Encode", "EPP"="EPP", "Bismark cytosine"="bismarkCytosine", "Bismark coverage"="bismarkCov")
RNB.FILTERING.SNP <- c("No filtering"="no", "3 SNPs"="3", "5 SNPs"="5", "Any SNPs"="any")
RNB.NORMALIZATION.METHODS=c("none", "bmiq", "illumina", "swan", "minfi.funnorm", "wm.dasen", "wm.nasen", "wm.betaqn", "wm.naten", "wm.nanet", "wm.nanes", "wm.danes", "wm.danet", "wm.danen", "wm.daten1", "wm.daten2", "wm.tost", "wm.fuks", "wm.swan")
RNB.NORMALIZATION.BG.METHODS <- c("none", "methylumi.noob", "methylumi.goob", "enmix.oob")
RNB.IMPUTATION.METHODS <- c("none", "mean.cpgs", "mean.samples", "random", "knn")
RNB.TRACKHUB.FORMATS <- c("bigBed", "bigWig")
RNB.SVA.NUM.METHODS <- c("leek", "be")
RNB.DIFFMETH.TEST.METHODS <- c("limma", "refFreeEWAS")
RNB.DIFFVAR.METHODS <- c("diffVar", "iEVORA")
RNB.COLSCHEMES.CATEGORY <- list(
	default=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
	extended=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#2166AC","#B2182B","#00441B","#40004B","#053061","#003D7C","#D50911")
)
RNB.COLSCHEMES.METH <- list(
	default=c("#AD0021","#909090","#39278C"),
	YlBl=c("#EDF8B1","#41B6C4","#081D58")
)

RNB.OPTION.DESC.TAB <- RnBeads:::rnb.options.description.table()
RNB.OPTION.DESC <- sapply(names(rnb.options()), FUN=function(x){
	if (is.element(x, rownames(RNB.OPTION.DESC.TAB))){
		return(RNB.OPTION.DESC.TAB[x, "desc"])
	} else {
		return("See the help pages: '?rnb.options'")
	}
})

RNB.GROUP.SIZE.RANGE  <- c(1, 20)
RNB.GROUP.COUNT.RANGE <- c(2, 20)

RNB.OPTION.PROFILES.PATH <- system.file(file.path("extdata", "optionProfiles"), package="RnBeads")
RNB.OPTION.PROFILES <- gsub("\\.xml$", "", list.files(path=RNB.OPTION.PROFILES.PATH, pattern="\\.xml$"))
################################################################################
# Choose local file or directory
# adapted from https://github.com/wleepang/shiny-directory-input
################################################################################
# Interactive Choosers for MacOS and Linux (For Windos this is already implemented in the utils package)
#' Choose a Folder Interactively (Mac OS X)
#'
#' Display a folder selection dialog under Mac OS X
#'
#' @param default which folder to show initially
#' @param caption the caption on the selection dialog
#'
#' @details
#' Uses an Apple Script to display a folder selection dialog.  With \code{default = NA},
#' the initial folder selection is determined by default behavior of the
#' "choose folder" Apple Script command.  Otherwise, paths are expanded with
#' \link{path.expand}.
#'
#' @return
#' A length one character vector, character NA if 'Cancel' was selected.
#'
if (Sys.info()['sysname'] == 'Darwin') {
	choose.dir <- function(default = NA, caption = NA, ...) {
		command = 'osascript'
		args = '-e "POSIX path of (choose folder{{prompt}}{{default}})"'

		if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
			prompt = sprintf(' with prompt \\"%s\\"', caption)
		} else {
			prompt = ''
		}
		args = sub('{{prompt}}', prompt, args, fixed = T)

		if (!is.null(default) && !is.na(default) && nzchar(default)) {
			default = sprintf(' default location \\"%s\\"', path.expand(default))
		} else {
			default = ''
		}
		args = sub('{{default}}', default, args, fixed = T)

		suppressWarnings({
			path = system2(command, args = args, stderr = TRUE)
			path <- path[length(path)]
		})
		if (!file.exists(path)) {
			# user canceled
		 	path = NA
		}

		return(path)
	}

	choose.files <- function(default = NA, caption = NA, ...) {
		command = 'osascript'
		args = '-e "POSIX path of (choose file{{prompt}}{{default}})"'

		if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
			prompt = sprintf(' with prompt \\"%s\\"', caption)
		} else {
			prompt = ''
		}
		args = sub('{{prompt}}', prompt, args, fixed = T)

		if (!is.null(default) && !is.na(default) && nzchar(default)) {
			default = sprintf(' default location \\"%s\\"', path.expand(default))
		} else {
			default = ''
		}
		args = sub('{{default}}', default, args, fixed = T)

		suppressWarnings({
			path = system2(command, args = args, stderr = TRUE)
			path <- path[length(path)]
		})
		if (!file.exists(path)) {
			# user canceled
			path = NA
		}

		return(path)
	}
} else if (Sys.info()['sysname'] == 'Linux') {
	library(tcltk)
	choose.dir <- tk_choose.dir
	choose.files <- tk_choose.files
}

#' Directory Selection Control
#'
#' Create a directory selection control to select a directory on the server
#'
#' @param inputId The \code{input} slot that will be used to access the value
#' @param label Display label for the control, or NULL for no label
#' @param value Initial value.  Paths are exapnded via \code{\link{path.expand}}.
#'
#' @details
#' This widget relies on \link{\code{choose.dir}} to present an interactive
#' dialog to users for selecting a directory on the local filesystem.  Therefore,
#' this widget is intended for shiny apps that are run locally - i.e. on the
#' same system that files/directories are to be accessed - and not from hosted
#' applications (e.g. from shinyapps.io).
#'
#' @return
#' A directory input control that can be added to a UI definition.
#'
#' @seealso
#' \link{updateDirectoryInput}, \link{readDirectoryInput}, \link[utils]{choose.dir}
directoryInput = function(inputId, label, value = NULL) {
	if (!is.null(value) && !is.na(value)) {
		value = path.expand(value)
	}

	tagList(
		singleton(
			tags$head(
				tags$script(src = 'js/directory_input_binding.js')
			)
		),

		div(
			class = 'form-group directory-input-container',
			tags$label(label),
			div(
				span(
					class = 'col-xs-9 col-md-11',
					style = 'padding-left: 0; padding-right: 5px;',
					div(
						class = 'input-group shiny-input-container',
						style = 'width:100%;',
						div(class = 'input-group-addon', icon('folder-o')),
						tags$input(
							id = sprintf('%s__chosen_dir', inputId),
							value = value,
							type = 'text',
							class = 'form-control directory-input-chosen-dir',
							readonly = 'readonly'
						)
					)
				),
				span(
					class = 'shiny-input-container',
					tags$button(
						id = inputId,
						class = 'btn btn-default directory-input',
						'...'
					)
				)
			)
		)
	)
}

#' Local File Selection Control
#'
#' Create a local file selection control to select a local file on the server
#'
#' @param inputId The \code{input} slot that will be used to access the value
#' @param label Display label for the control, or NULL for no label
#' @param value Initial value.  Paths are exapnded via \code{\link{path.expand}}.
#'
#' @details
#' This widget relies on \link{\code{choose.files}} to present an interactive
#' dialog to users for selecting files on the local filesystem.  Therefore,
#' this widget is intended for shiny apps that are run locally - i.e. on the
#' same system that files/directories are to be accessed - and not from hosted
#' applications (e.g. from shinyapps.io).
#'
#' @return
#' A local file input control that can be added to a UI definition.
#'
#' @seealso
#' \link{updateLocalFileInput}, \link{readLocalFileInput}, \link[utils]{choose.files}
localFileInput = function(inputId, label, value = NULL) {
	if (!is.null(value) && !is.na(value)) {
		value = path.expand(value)
	}

	tagList(
		singleton(
			tags$head(
				tags$script(src = 'js/localfile_input_binding.js')
			)
		),

		div(
			class = 'form-group localfile-input-container',
			tags$label(label),
			div(
				span(
					class = 'col-xs-9 col-md-11',
					style = 'padding-left: 0; padding-right: 5px;',
					div(
						class = 'input-group shiny-input-container',
						style = 'width:100%;',
						div(class = 'input-group-addon', icon('file-o')),
						tags$input(
							id = sprintf('%s__chosen_file', inputId),
							value = value,
							type = 'text',
							class = 'form-control localfile-input-chosen-file',
							readonly = 'readonly'
						)
					)
				),
				span(
					class = 'shiny-input-container',
					tags$button(
						id = inputId,
						class = 'btn btn-default localfile-input',
						'...'
					)
				)
			)
		)
	)
}

#' Change the value of a directoryInput on the client
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object.
#' @param value A directory path to set
#' @param ... Additional arguments passed to \link{\code{choose.dir}}.  Only used
#'    if \code{value} is \code{NULL}.
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object.  For \code{directoryInput} objects, this changes the value displayed
#' in the text-field and triggers a client-side change event.  A directory
#' selection dialog is not displayed.
#'
updateDirectoryInput = function(session, inputId, value = NULL, ...) {
  if (is.null(value)) {
    value = choose.dir(...)
  }
  session$sendInputMessage(inputId, list(chosen_dir = value))
}

#' Change the value of a localFileInput on the client
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object.
#' @param value A file path to set
#' @param ... Additional arguments passed to \link{\code{choose.files}}.  Only used
#'    if \code{value} is \code{NULL}.
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object.  For \code{localFileInput} objects, this changes the value displayed
#' in the text-field and triggers a client-side change event.  A file
#' selection dialog is not displayed.
#'
updateLocalFileInput = function(session, inputId, value = NULL, ...) {
  if (is.null(value)) {
    value = choose.files(...)
  }
  session$sendInputMessage(inputId, list(chosen_file = value))
}

#' Read the value of a directoryInput
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object
#'
#' @details
#' Reads the value of the text field associated with a \code{directoryInput}
#' object that stores the user selected directory path.
#'
readDirectoryInput = function(session, inputId) {
  session$input[[sprintf('%s__chosen_dir', inputId)]]
}

#' Read the value of a localFileInput
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object
#'
#' @details
#' Reads the value of the text field associated with a \code{localFileInput}
#' object that stores the user selected file path.
#'
readLocalFileInput = function(session, inputId) {
  session$input[[sprintf('%s__chosen_file', inputId)]]
}

observeDirectoryInput <- function(input, session, inputId){
	observeEvent(ignoreNULL=TRUE, eventExpr={input[[inputId]]},
		handlerExpr={
			if (input[[inputId]] > 0) {	      
				# launch the directory selection dialog with initial path read from the widget
				defaultPath <- readDirectoryInput(session, inputId)
				if (defaultPath == "[NULL]") defaultPath <- NA
				path = choose.dir(default = defaultPath)
				# update the widget value
				updateDirectoryInput(session, inputId, value=path)
			}
		}
	)
}

observeLocalFileInput <- function(input, session, inputId){
	observeEvent(ignoreNULL=TRUE, eventExpr={input[[inputId]]},
		handlerExpr={
			if (input[[inputId]] > 0) {
				# launch the file selection dialog with initial path read from the widget 
				filt <- matrix(c("comma-separated", ".csv", "tab-separated", ".tsv", "Text", ".txt", "All files", "*"), 4, 2, byrow = TRUE)
				defaultPath <- readLocalFileInput(session, inputId)
				if (defaultPath == "[NULL]") defaultPath <- NA
				path = choose.files(default = defaultPath, filters=filt)
				# update the widget value
				updateLocalFileInput(session, inputId, value=path)
			}
		}
	)
}

################################################################################
# Little Helpers
################################################################################
plotColPal <- function(colpal){
	require(plotrix)
	par(mar=c(0,0,0,0))
	plot.new()
	gradient.rect(0,0,1,1,col=colpal,nslices=length(colpal),gradient="x",border=NA)
}

getRnbStatusFromLog <- function(logFile){
	if (!file.exists(logFile)) return(NULL)
	ll <- readLines(logFile)
	res <- data.frame(
		module=RNB.MODULES,
		status=NA,
		scheduled=NA
	)
	rownames(res) <- RNB.MODULES
	for (mm in RNB.MODULES){
		if(sum(grepl(paste0("COMPLETED ", RNB.MODULES.LOG.MSG[mm]), ll))>0){
			res[mm, "status"] <- "completed"
		} else if(sum(grepl(paste0("STARTED ", RNB.MODULES.LOG.MSG[mm]), ll))>0){
			res[mm, "status"] <- "started"
		}
	}
	#check the options if the anlysis is scheduled
	res["data_import", "scheduled"]              <- rnb.getOption("import")
	res["quality_control", "scheduled"]          <- rnb.getOption("qc")
	res["preprocessing", "scheduled"]            <- rnb.getOption("preprocessing")
	res["tracks_and_tables", "scheduled"]        <- rnb.getOption("export.to.csv") || rnb.getOption("export.to.bed") || length(rnb.getOption("export.to.trackhub")) > 0
	res["covariate_inference", "scheduled"]      <- rnb.getOption("inference")
	res["exploratory_analysis", "scheduled"]     <- rnb.getOption("exploratory")
	res["differential_methylation", "scheduled"] <- rnb.getOption("import")
	return(res)
}

checkReportDir <- function(repDir){
	res <- list(valid=FALSE)
	res$reportHtml <- rep(NA, length(RNB.MODULES))
	names(res$reportHtml) <- RNB.MODULES
	res$moduleStatus <- list(NULL)
	res$createdByDj <- FALSE # check if this report directory was created by RnBeadsDJ

	if (!dir.exists(repDir)) return(res)
	contents <- list.files(repDir, include.dirs=TRUE, all.files=TRUE)
	res$valid <- all(c("configuration", "analysis.log") %in% contents)
	if (!res$valid) return(res)
	res$createdByDj <- is.element(".RnBeadsDJ", contents) 
	res$valid <- res$createdByDj || all(c("index.html", "analysis_options.xml") %in% contents)
	if (!res$valid) return(res)
	htmlModules <- paste0(RNB.MODULES, ".html")
	htmlModules.exist <- htmlModules %in% contents
	res$reportHtml[htmlModules.exist] <- htmlModules[htmlModules.exist]
	res$logFile <- file.path(repDir, "analysis.log")
	res$moduleStatus <- getRnbStatusFromLog(res$logFile)
	return(res)
}
# mark a report directory as created by RnBeadsDJ
markDirDJ <- function(repDir){
	repStatus <- checkReportDir(repDir)
	if (!repStatus$createdByDj){
		file.create(file.path(repDir, ".RnBeadsDJ"))
	}
	invisible(NULL)
}
################################################################################
# UI configuration
################################################################################
ui <- tagList(useShinyjs(), navbarPage(
	windowTitle="RnBeadsDJ",
	tags$p(tags$a(href="https://rnbeads.org", tags$img(width=145, height=50, src="img/rnbeads_logo.png")), "DJ"),
	# tabPanel("Sandbox", icon=icon("dropbox"),
	# 	verbatimTextOutput("sandboxOut"),
	# 	localFileInput('sandboxIn2', label='select a local file'),
	# 	verbatimTextOutput("sandboxOut2")
	# ),
	tabPanel("About", icon=icon("book"),
		tags$h1("Welcome to the RnBeads Data Juggler"),
		tags$p("Here, you can configure and run your RnBeads analyses. You can ..."),
		tags$ul(
			tags$li("... run new analyses by specifying a non-existing data directory in the", "'Analysis'", "tab."),
			tags$li("... view the status of an RnBeads analysis by specifying an existing report directory in the", "'Analysis'", "tab."),
			tags$li("... configure, load and save option settings for your analyses in the", "'Analysis Options'", "tab."),
			tags$li("... run the complete RnBeads pipeline using the", "'Input'",  " and 'Run'", "tabs."),
			tags$li("... run individual RnBeads modules for new or existing RnBeads analysis runs via the", "'Modules'", "tab.")
		)
	),
	tabPanel(tags$div(title="Specify general parameters, such as report (output) directories, number of cores to use, etc. You can configure new analyses or inspect the status of previous RnBeads runs.", icon("bar-chart"), "Analysis"),
		sidebarPanel(
			tags$div(title="Specify parent directory where the analysis reports will be stored.", directoryInput('outDir', label = 'Select analysis directory', value = '~')),
			tags$div(title="Specify the report directory (subdirectory of the analysis directory). Specify a non-existing directory for a new RnBeads analysis or an existing one to inspect the status and/or continue a previously started analysis.", textInput("reportSubDir", "Choose the name of the report directory", "rnbeads_report")),
			tags$div(title="Select the number of cores to use for the analysis. Note that selecting more cores will speed up analysis, but might lead to significantly increased memory usage.", sliderInput('numCores', "Select the number of cores to use", min=1, max=detectCores(), value=1, step=1)),
			tags$div(title="Select a color scheme for the plots generated by RnBeads.", selectInput('ggplotTheme', "Select a theme for the plots", c("Black & White"="bw", "Grey"="grey")))
		),
		mainPanel(
			uiOutput("anaStatus")
		)
	),
	tabPanel(tags$div(title="Specify input parameters, such as where to find data and which format and assay platform to expect.", icon("sign-in"), "Input"),
		sidebarPanel(
			tags$div(title="Select a file which contains a valid annotation table with rows corresponding to each sample in the analysis.", localFileInput("sampleAnnotFile", "Select sample annotation file")),
			tags$div(title="Select the character that is used to separate the columns in the sample annotation table.", selectInput("rnbOptsI.import.table.separator", "Separator:", RNB.TABLE.SEPS)),
			tags$div(title="Select the experimental platform that is used for DNA methylation quantification.", selectInput("platform", "Platform", RNB.PLATFORMS)),
			tags$div(title="Specify the directory where the input methylation data (IDAT or methylation call files) can be found.", directoryInput('dataDir', label='Select input data directory'))
		),
		mainPanel(
			uiOutput("inputStatus")
		)
	),
	tabPanel(tags$div(title="Specify analysis parameters in order to finetune your methylation analysis. Parameter settings can be configured, saved and loaded.", icon("sliders"), "Analysis Options"),
		fluidRow(
			column(9,
				wellPanel(
					tags$h3("Load Option Profile"),
					fluidRow(
						column(3, tags$div(title="Load option profile store from a preexisting analysis. Existing analysis directory must be specified.", wellPanel(
							actionButton("loadOptsAnaDirDo", "Load from Analysis Directory")
						))),
						column(5, tags$div(title="Load option profile from XML file.", wellPanel(
							localFileInput("loadOptsXmlFile", "XML file"),
							actionButton("loadOptsXmlDo", "Load from XML")
						))),
						column(4, tags$div(title="Load option profile from preset. Several options for more or less comprehensive and resource-demanding analyses are available.", wellPanel(
							selectInput("loadOptsProfileSel", "Predefined Option Profile", RNB.OPTION.PROFILES),
							actionButton("loadOptsProfileDo", "Load Option Profile")
						)))
					)
				)
			),
			column(3,
				tags$div(title="Save the current option profile to an XML file. This file can be loaded later for running a different analysis with the same parameters.", wellPanel(
					tags$h3("Saving Option Profile"),
					downloadButton("saveOptsXml", "Save to XML")
				))
			)
		),
		tabsetPanel(
			tabPanel("General",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["analysis.name"], tags$code("analysis.name"))
							),
							tags$td(
								textInput("rnbOptsI.analysis.name", NULL, "RnBeads Analysis")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.analysis.name")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["assembly"], tags$code("assembly"))
							),
							tags$td(
								selectInput("rnbOptsI.assembly", NULL, RNB.ASSEMBLIES, selected="hg38")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.assembly", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["region.types"], tags$code("region.types"))
							),
							tags$td(
								uiOutput("selRegionTypes")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.region.types", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["identifiers.column"], tags$code("identifiers.column"))
							),
							tags$td(
								uiOutput('selColumn.id')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.identifiers.column", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["min.group.size"], tags$code("min.group.size"))
							),
							tags$td(
								sliderInput("rnbOptsI.min.group.size", NULL, min=RNB.GROUP.SIZE.RANGE[1], max=RNB.GROUP.SIZE.RANGE[2], value=2)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.min.group.size")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["max.group.count"], tags$code("max.group.count"))
							),
							tags$td(
								sliderInput("rnbOptsI.max.group.count", NULL, min=RNB.GROUP.COUNT.RANGE[1], max=RNB.GROUP.COUNT.RANGE[2], value=RNB.GROUP.COUNT.RANGE[2])
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.max.group.count")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["colors.category"], tags$code("colors.category"))
							),
							tags$td(
								selectInput("rnbOptsI.colors.category", NULL, names(RNB.COLSCHEMES.CATEGORY))
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.colors.category"),
								plotOutput("rnbOptsOP.colors.category", height="30px")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["colors.meth"], tags$code("colors.meth"))
							),
							tags$td(
								selectInput("rnbOptsI.colors.meth", NULL, names(RNB.COLSCHEMES.METH))
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.colors.meth"),
								plotOutput("rnbOptsOP.colors.meth", height="30px")
							)
						)
					)
				)
			),
			tabPanel("Import",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["import.default.data.type"], tags$code("import.default.data.type"))
							),
							tags$td(
								"Defined per 'Platform' in the 'Input' section"
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.import.default.data.type")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["import.table.separator"], tags$code("import.table.separator"))
							),
							tags$td(
								"Defined per 'Separator' in the 'Input' section"
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.import.table.separator")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["import.bed.style"], tags$code("import.bed.style"))
							),
							tags$td(
								selectInput("rnbOptsI.import.bed.style", NULL, RNB.BED.STYLES)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.import.bed.style")
							)
						)
					)
				)
			),
			tabPanel("Quality Control",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["qc"], tags$code("qc"))
							),
							tags$td(
								checkboxInput("rnbOptsI.qc", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.qc")
							)
						)
					)
				)
			),
			tabPanel("Preprocessing",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["preprocessing"], tags$code("preprocessing"))
							),
							tags$td(
								checkboxInput("rnbOptsI.preprocessing", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.preprocessing")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.coverage.threshold"], tags$code("filtering.coverage.threshold"))
							),
							tags$td(
								sliderInput("rnbOptsI.filtering.coverage.threshold", NULL, min=1, max=100, value=5)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.coverage.threshold")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.low.coverage.masking"], tags$code("filtering.low.coverage.masking"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.low.coverage.masking", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.low.coverage.masking")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.high.coverage.outliers"], tags$code("filtering.high.coverage.outliers"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.high.coverage.outliers", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.high.coverage.outliers")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.missing.value.quantile"], tags$code("filtering.missing.value.quantile"))
							),
							tags$td(
								sliderInput("rnbOptsI.filtering.missing.value.quantile", NULL, min=0, max=1, value=1)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.missing.value.quantile")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.greedycut"], tags$code("filtering.greedycut"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.greedycut", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.greedycut")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.sex.chromosomes.removal"], tags$code("filtering.sex.chromosomes.removal"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.sex.chromosomes.removal", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.sex.chromosomes.removal")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.snp"], tags$code("filtering.snp"))
							),
							tags$td(
								selectInput("rnbOptsI.filtering.snp", NULL, RNB.FILTERING.SNP, selected="3")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.snp")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.cross.reactive"], tags$code("filtering.cross.reactive"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.cross.reactive", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.cross.reactive")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.whitelist"], tags$code("filtering.whitelist"))
							),
							tags$td(
									localFileInput("rnbOptsI.filtering.whitelist", NULL),
									actionButton("rnbOptsResetWhitelist", "Reset whitelist")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.whitelist")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.blacklist"], tags$code("filtering.blacklist"))
							),
							tags$td(
									localFileInput("rnbOptsI.filtering.blacklist", NULL),
									actionButton("rnbOptsResetblacklist", "Reset blacklist")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.blacklist")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["normalization.method"], tags$code("normalization.method"))
							),
							tags$td(
								selectInput("rnbOptsI.normalization.method", NULL, RNB.NORMALIZATION.METHODS, selected="wm.dasen")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.normalization.method")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["normalization.background.method"], tags$code("normalization.background.method"))
							),
							tags$td(
								selectInput("rnbOptsI.normalization.background.method", NULL, RNB.NORMALIZATION.BG.METHODS, selected="none")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.normalization.background.method")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["imputation.method"], tags$code("imputation.method"))
							),
							tags$td(
								selectInput("rnbOptsI.imputation.method", NULL, RNB.IMPUTATION.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.imputation.method")
							)
						)
					)
				)
			),
			tabPanel("Export",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["export.to.bed"], tags$code("export.to.bed"))
							),
							tags$td(
								checkboxInput("rnbOptsI.export.to.bed", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.export.to.bed")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["export.to.csv"], tags$code("export.to.csv"))
							),
							tags$td(
								checkboxInput("rnbOptsI.export.to.csv", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.export.to.csv")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["export.to.trackhub"], tags$code("export.to.trackhub"))
							),
							tags$td(
								selectInput("rnbOptsI.export.to.trackhub", NULL, RNB.TRACKHUB.FORMATS, multiple=TRUE, selected=RNB.TRACKHUB.FORMATS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.export.to.trackhub", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["export.types"], tags$code("export.types"))
							),
							tags$td(
								uiOutput('selRegions.export')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.export.types", placeholder=TRUE)
							)
						)
					)
				)
			),
			tabPanel("Inference",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference"], tags$code("inference"))
							),
							tags$td(
								checkboxInput("rnbOptsI.inference", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.age.prediction"], tags$code("inference.age.prediction"))
							),
							tags$td(
								checkboxInput("rnbOptsI.inference.age.prediction", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.age.prediction")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.age.column"], tags$code("inference.age.column"))
							),
							tags$td(
								uiOutput('selColumn.agepred')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.age.column", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.targets.sva"], tags$code("inference.targets.sva"))
							),
							tags$td(
								uiOutput('selColumn.sva')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.targets.sva", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.sva.num.method"], tags$code("inference.sva.num.method"))
							),
							tags$td(
								selectInput("rnbOptsI.inference.sva.num.method", NULL, RNB.SVA.NUM.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.sva.num.method", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.reference.methylome.column"], tags$code("inference.reference.methylome.column"))
							),
							tags$td(
								uiOutput('selColumn.cellTypeRef')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.reference.methylome.column", placeholder=TRUE)
							)
						)
					)
				)
			),
			tabPanel("Exploratory",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory"], tags$code("exploratory"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.columns"], tags$code("exploratory.columns"))
							),
							tags$td(
								uiOutput('selColumn.ex')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.intersample"], tags$code("exploratory.intersample"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.intersample", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.intersample")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.beta.distribution"], tags$code("exploratory.beta.distribution"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.beta.distribution", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.beta.distribution")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.correlation.qc"], tags$code("exploratory.correlation.qc"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.correlation.qc", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.correlation.qc")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.region.profiles"], tags$code("exploratory.region.profiles"))
							),
							tags$td(
								uiOutput('selRegionProfiles.ex')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.region.profiles", placeholder=TRUE)
							)
						)
					)
				)
			),
			tabPanel("Differential",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential"], tags$code("differential"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.comparison.columns"], tags$code("differential.comparison.columns"))
							),
							tags$td(
								uiOutput('selColumn.diff')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.comparison.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["covariate.adjustment.columns"], tags$code("covariate.adjustment.columns"))
							),
							tags$td(
								uiOutput('selAdjColumns')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.covariate.adjustment.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.site.test.method"], tags$code("differential.site.test.method"))
							),
							tags$td(
								selectInput("rnbOptsI.differential.site.test.method", NULL, RNB.DIFFMETH.TEST.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.site.test.method", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.report.sites"], tags$code("differential.report.sites"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.report.sites", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.report.sites")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.variability"], tags$code("differential.variability"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.variability", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.variability")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.variability.method"], tags$code("differential.variability.method"))
							),
							tags$td(
								selectInput("rnbOptsI.differential.variability.method", NULL, RNB.DIFFVAR.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.variability.method")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.enrichment.go"], tags$code("differential.enrichment.go"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.enrichment.go", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.enrichment.go")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.enrichment.lola"], tags$code("differential.enrichment.lola"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.enrichment.lola", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.enrichment.lola")
							)
						)
					)
				)
			),
			tabPanel("R Object",
				verbatimTextOutput("rnbOpts")
			)
		)
	),
	tabPanel(tags$div(title="Run a full RnBeads analysis from scratch once you specified your analysis parameters.", icon("play"), "Run"),
		uiOutput("runHint"),
		actionButton("runRnb", "Run Analysis", class="btn-primary"),
		textOutput("runStatusMsg")
	),
	tabPanel(tags$div(title="Run individual steps of an RnBeads analysis using preexisting datasets as input.", icon("cubes"), "Modules"),
		tabsetPanel(
			tabPanel("Data Import",
				sidebarPanel(
					tags$h1("Import from ..."),
					wellPanel(
						tags$h4("... new Dataset"),
						tags$p("Import using the inputs and options specified in other tabs."),
						uiOutput("modImportNew.about"),
						actionButton("modImportNew", "Import New Data", class="btn-primary")
					),
					wellPanel(
						tags$h4("... existing Analysis Directory"),
						uiOutput("modImportAnaDir.about"),
						actionButton("modImportAnaDir", "Import from Analysis Directory", class="btn-primary")
					),
					wellPanel(
						tags$h4("... RnBeads Objects"),
						tags$div(title="Specify the directory of an existing RnBSet object. Be sure to select a directory (not an RData file).", directoryInput('modImportRnBSetDir', label = 'Select directory containing an RnBSet object', value = '~')),
						# localFileInput('modImportOptionsFile', label = 'Select Options XML file (optional)'),
						actionButton("modImportRObjects", "Import", class="btn-primary")
					),
					wellPanel(
						tags$h4("Reset"),
						actionButton("modImportReset", "Unload dataset")
					)
				),
				mainPanel(
					tags$h1("Loading Status"),
					uiOutput("modImport.status"),
					uiOutput("rnbSetInfo")
				)
			),
			tabPanel("Quality Control",
				sidebarPanel(
					checkboxInput("modQC.overwrite", "Overwrite existing report", value=FALSE),
					actionButton("modQC.run", "Run Quality Control", class="btn-primary")
				),
				mainPanel(
					uiOutput("modQC.status")
				)
			),
			tabPanel("Preprocessing",
				sidebarPanel(
					checkboxInput("modPreprocessing.overwrite", "Overwrite existing report", value=FALSE),
					checkboxInput("modPreprocessing.save", "Save RnBSet object", value=TRUE),
					actionButton("modPreprocessing.run", "Run Preprocessing", class="btn-primary")
				),
				mainPanel(
					uiOutput("modPreprocessing.status")
				)
			),
			tabPanel("Tracks and Tables",
				sidebarPanel(
					checkboxInput("modTNT.overwrite", "Overwrite existing report", value=FALSE),
					actionButton("modTNT.run", "Run Tracks and Tables", class="btn-primary")
				),
				mainPanel(
					uiOutput("modTNT.status")
				)
			),
			tabPanel("Covariate Inference",
				sidebarPanel(
					checkboxInput("modInference.overwrite", "Overwrite existing report", value=FALSE),
					checkboxInput("modInference.save", "Save RnBSet object", value=TRUE),
					actionButton("modInference.run", "Run Covariate Inference", class="btn-primary")
				),
				mainPanel(
					uiOutput("modInference.status")
				)
			),
			tabPanel("Exploratory Analysis",
				sidebarPanel(
					checkboxInput("modExploratory.overwrite", "Overwrite existing report", value=FALSE),
					actionButton("modExploratory.run", "Run Exploratory Analysis", class="btn-primary")
				),
				mainPanel(
					uiOutput("modExploratory.status")
				)
			),
			tabPanel("Differential Methylation",
				sidebarPanel(
					checkboxInput("modDifferential.overwrite", "Overwrite existing report", value=FALSE),
					actionButton("modDifferential.run", "Run Differential Methylation", class="btn-primary")
				),
				mainPanel(
					uiOutput("modDifferential.status")
				)
			)
		)
	)
))
################################################################################
# Server configuration
################################################################################
server <- function(input, output, session) {
	############################################################################
	# SANDBOX
	# file and directory handling
	# shinyDirChoose(input, "outDirSel", roots=c(disk="/", home="~", wd="."), session=session)

	# dynamic outputs
	observeLocalFileInput(input, session, 'sandboxIn2')
	output$sandboxOut2 <- renderPrint({readLocalFileInput(session, 'sandboxIn2')})
	############################################################################

	# output$outDirSel <- renderPrint({input$outDirSel})
	# output$outDirSel <- renderPrint({readDirectoryInput(session, 'outDirSel')})


	############################################################################
	# Analysis status
	############################################################################
	refreshTimer <- reactiveTimer(REFRESH.TIME)
	observeDirectoryInput(input, session, 'outDir')
	reportDir <- reactive({
		file.path(readDirectoryInput(session, 'outDir'), input$reportSubDir)
	})
	setNumCoresFromSlider <- reactive({
		res <- input$numCores
		if (res > 1){
			parallel.setup(res)
		} else {
			parallel.teardown()
		}
		logger.close()
		res
	})
	anaStatus <- reactive({
		refreshTimer()
		res <- list(status="invalid", statusTab=NULL, rnbSet.paths=c(), logFile=NA)
		if (dir.exists(reportDir())){
			reportStatus <- checkReportDir(reportDir())
			if (reportStatus$valid){
				#Valid existing report directory
				statusTab <- data.frame(
					Module=names(RNB.MODULES),
					Report=NA,
					Status=NA
				)
				rownames(statusTab) <- RNB.MODULES
				for (mm in RNB.MODULES){
					if (!is.na(reportStatus$reportHtml[mm])){
						statusTab[mm, "Report"] <- as.character(tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml[mm])), tags$code("html")))
						if (reportStatus$moduleStatus[mm, "status"]=="completed"){
							statusTab[mm, "Status"] <- as.character(tags$span(style="color:green", icon("check")))
						} else if (reportStatus$moduleStatus[mm, "status"]=="started"){
							statusTab[mm, "Status"] <- as.character(tags$span(style="color:orange", icon("play")))
						}
						
					} else {
						statusTab[mm, "Status"] <- as.character(tags$span(style="color:red", icon("times")))
					}
				}
				res$rnbSet.paths <- file.path(reportDir(), grep("rnbSet_", list.dirs(reportDir(), full.names=FALSE, recursive=FALSE), value=TRUE))
				if (!modImportStatus$dataset.loaded) shinyjs::enable("modImportAnaDir")
				res$logFile <- reportStatus$logFile
				res$status <- "reportDir"
				res$statusTab <- statusTab
				shinyjs::enable("loadOptsAnaDirDo")
			} else {
				# Existing directory, but no valid RnBeads report
				shinyjs::disable("modImportAnaDir")
				shinyjs::disable("loadOptsAnaDirDo")
			}	
		} else {
			# New RnBeads analysis
			shinyjs::disable("modImportAnaDir")
			shinyjs::disable("loadOptsAnaDirDo")
			res$status <- "new"
		}
		res
	})
	output$anaStatus <- renderUI({
		curStatus <- anaStatus()
		if (curStatus$status=="new"){
			tagList(
				tags$h1("Configure new analysis"),
				tags$p(
					"Configure the input and analysis options by clicking on the tabs above",
					"and selecting corresponding values. Finally go to 'Run' to start your analysis."
				),
				tags$p(
					"A new report will be created in", tags$code(reportDir())
				)
			)
		} else if (curStatus$status=="reportDir"){
			tagList(
				tags$h1("Existing RnBeads analysis"),
				renderTable(curStatus$statusTab, sanitize.text.function=function(x){x}, striped=TRUE, hover=TRUE, bordered=TRUE),
				tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("file-text-o"), tags$a(href=paste0("file://", curStatus$logFile), "Log File"))))
			)
		} else if (curStatus$status=="invalid"){
			tagList(
				tags$h1("Invalid analysis directory"),
				tags$p("Please select a non-existing analysis directory for new analyses or an existing RnBeads report directory.")
			)
		}
	})

	############################################################################
	# Input
	############################################################################
	observeLocalFileInput(input, session, 'sampleAnnotFile')
	observeDirectoryInput(input, session, 'dataDir')

	tabSep <- reactive({input$rnbOptsI.import.table.separator})
	sampleAnnotFile <- reactive({
		readLocalFileInput(session, 'sampleAnnotFile')
	})
	inputDataDir <- reactive({
		readDirectoryInput(session, 'dataDir')
	})

	sannot.fromFile <- reactive({
		sampleAnnotFn <- sampleAnnotFile()
		if (is.null(sampleAnnotFn) || sampleAnnotFn=="") return(NULL)
		tryCatch(
			{
				RnBeads::read.sample.annotation(sampleAnnotFn, sep=tabSep())
			},
			error = function(err) {
				data.frame(Error="Unable to load table", Why=err$message)
				NULL
			}
		)
	})
	sannot <- reactive({
		# first look if an RnBSet object has been loaded. If so, take the sample annotation from there
		# otherwise require the user to upload an annotation table via file input
		if (!is.null(rnbData$rnbSet)){
			pheno(rnbData$rnbSet)
		} else {
			sannot.fromFile()
		}
	})
	sannot.nsamples <- reactive({
		if (is.null(sannot())){
			0
		} else {
			nrow(sannot())
		}
	})
	sannot.cols <- reactive({
		cnames <- colnames(sannot())
		if (length(cnames)==2 && cnames==c("Error", "Why")) cnames <- NULL #error --> no column names
		cnames
	})
	sannot.cols.plusNone <- reactive(c("[None]", sannot.cols()))
	sannot.cols.plusDefault <- reactive(c("[default]", sannot.cols()))
	sannot.cols.grps <- reactive({
		# print("DEBUG: updated sannot.cols.grps")
		depDummy <- input$rnbOptsI.min.group.size
		depDummy <- input$rnbOptsI.max.group.count
		if(length(sannot.cols())>0) {
			names(rnb.sample.groups(sannot()))
		} else {
			NULL
		}
	})
	sannot.cols.plusAutomatic <- reactive(c("[automatic]", sannot.cols()))
	output$sampleAnnotContent <- renderTable({sannot.fromFile()}, striped=TRUE, hover=TRUE, bordered=TRUE)

	isBiseq <- reactive({
		selPlatform <- input$platform
		if (!is.null(rnbData$rnbSet) && (inherits(rnbData$rnbSet, "RnBeadSet") || inherits(rnbData$rnbSet, "RnBiseqSet"))){
			res <- inherits(rnbData$rnbSet, "RnBiseqSet")
		} else {
			res <- selPlatform == "biseq"
		}
		if (res) {
			shinyjs::enable("rnbOptsI.assembly")
			shinyjs::enable("rnbOptsI.import.bed.style")
			shinyjs::enable("rnbOptsI.filtering.coverage.threshold")
			shinyjs::enable("rnbOptsI.filtering.low.coverage.masking")
			shinyjs::enable("rnbOptsI.filtering.high.coverage.outliers")
			updateCheckboxInput(session, "rnbOptsI.filtering.greedycut", value=FALSE)
			shinyjs::disable("rnbOptsI.filtering.greedycut")
			# shinyjs::disable("rnbOptsI.filtering.whitelist")
			# shinyjs::disable("rnbOptsI.filtering.blacklist")
			shinyjs::disable("rnbOptsI.filtering.cross.reactive")
			shinyjs::disable("rnbOptsI.normalization.method")
			shinyjs::disable("rnbOptsI.normalization.background.method")
			shinyjs::disable("rnbOptsI.exploratory.correlation.qc")
		} else {
			shinyjs::disable("rnbOptsI.assembly")
			shinyjs::disable("rnbOptsI.import.bed.style")
			shinyjs::disable("rnbOptsI.filtering.coverage.threshold")
			shinyjs::disable("rnbOptsI.filtering.low.coverage.masking")
			shinyjs::disable("rnbOptsI.filtering.high.coverage.outliers")
			updateCheckboxInput(session, "rnbOptsI.filtering.greedycut", value=TRUE)
			shinyjs::enable("rnbOptsI.filtering.greedycut")
			shinyjs::enable("rnbOptsI.filtering.cross.reactive")
			# shinyjs::enable("rnbOptsI.filtering.whitelist")
			# shinyjs::enable("rnbOptsI.filtering.blacklist")
			shinyjs::enable("rnbOptsI.normalization.method")
			shinyjs::enable("rnbOptsI.normalization.background.method")
			shinyjs::enable("rnbOptsI.exploratory.correlation.qc")
		}
		res
	})

	output$selColumn.id <- renderUI({
		selectInput('rnbOptsI.identifiers.column', NULL, sannot.cols.plusNone(), selected=optionSettingObserver$selColumn.id)
	})
	output$selColumn.agepred <- renderUI({
		selectInput('rnbOptsI.inference.age.column', NULL, sannot.cols.plusDefault(), selected=optionSettingObserver$selColumn.agepred)
	})
	output$selColumn.sva<- renderUI({
		selectInput('rnbOptsI.inference.targets.sva', NULL, sannot.cols.grps(), multiple=TRUE, selected=optionSettingObserver$selColumn.sva)
	})
	output$selColumn.cellTypeRef <- renderUI({
		selectInput('rnbOptsI.inference.reference.methylome.column', NULL, sannot.cols.plusNone(), selected=optionSettingObserver$selColumn.cellTypeRef)
	})
	output$selColumn.ex <- renderUI({
		selectInput('rnbOptsI.exploratory.columns', NULL, sannot.cols.plusAutomatic(), multiple=TRUE, selected=optionSettingObserver$selColumn.ex)
	})
	output$selColumn.diff <- renderUI({
		selectInput('rnbOptsI.differential.comparison.columns', NULL, sannot.cols.plusAutomatic(), multiple=TRUE, selected=optionSettingObserver$selColumn.diff)
	})
	output$selAdjColumns <- renderUI({
		selectInput('rnbOptsI.covariate.adjustment.columns', NULL, sannot.cols(), multiple=TRUE, selected=optionSettingObserver$selColumn.adj)
	})

	output$inputStatus <- renderUI({
		res <- list()
		showAnnotTab <- TRUE
		sampleAnnotFn <- sampleAnnotFile()
		if (is.null(sampleAnnotFn) || sampleAnnotFn=="") showAnnotTab <- FALSE
		if (showAnnotTab){
			res <- c(res,
				list(tags$h1("Preview of the sample annotation table")),
				list(tableOutput("sampleAnnotContent"))
			)
		} else {
			res <- c(res,
				list(tags$p("No sample annotation loaded"))
			)
		}
		tagList(res)
	})

	############################################################################
	# Analysis Options
	############################################################################
	# autoOpts <- c("analysis.name", "import.table.separator")
	# for (oo in autoOpts){
	# 	output[[paste0("rnbOptsO.", oo)]] <- renderText({
	# 		args <- list()
	# 		args[[oo]] <- input[[paste0("rnbOptsI.", oo)]]
	# 		do.call("rnb.options", args)
	# 		rnb.getOption(oo)
	# 	})
	# }

	# a helper container for storing variable, option-related information
	optionSettingObserver <- reactiveValues(
		group.count.range.max=RNB.GROUP.COUNT.RANGE[2],
		colors.category.list=RNB.COLSCHEMES.CATEGORY,
		selRegs=NULL,
		selRegs.export=NULL,
		selRegs.exploratory.profiles=NULL,
		selRegs.differential=NULL,
		selColumn.id="[None]",
		selColumn.agepred="[default]",
		selColumn.sva=character(0),
		selColumn.cellTypeRef="[None]",
		selColumn.ex="[automatic]",
		selColumn.diff="[automatic]",
		selColumn.adj=character(0)
	)

	optList <- reactive({
		depDummy <- assemblySel()
		depDummy <- isBiseq() #dummy to create dependency upon update of platform selection
		depDummy <- optsImportDataType()
		# auto update on any option change to rnbOptsI
		for (oo in grep("^rnbOptsI.", names(input), value=TRUE)){
			input[[oo]]
		}
		rnb.options()
	})
	output$rnbOpts <- renderPrint({
		optList()
	})
	output$rnbOptsO.analysis.name <- renderText({
		rnb.options(analysis.name=input$rnbOptsI.analysis.name)
		rnb.getOption("analysis.name")
	})
	#default setting for assembly if platform is not bisulfite and thus no assembly input has been given
	assemblySel <- reactive({
		interfaceSetting <- input$rnbOptsI.assembly
		res <- "hg19"
		if (isBiseq() && !is.null(interfaceSetting)){
			res <- interfaceSetting
		}
		rnb.options(assembly=res)
		res
	})
	output$rnbOptsO.assembly <- renderText({
		depDummy <- assemblySel()
		rnb.getOption("assembly")
	})
	regTypes.all <- reactive({
		depDummy <- assemblySel() # dummy to update on dependency: assembly
		rnb.region.types(rnb.getOption("assembly"))
	})
	output$selRegionTypes <- renderUI({
		rrs <- optionSettingObserver$selRegs
		if (is.null(rrs)){
			rrs <- regTypes.all()
			optionSettingObserver$selRegs <- rrs
		}
		selectInput('rnbOptsI.region.types', NULL, regTypes.all(), multiple=TRUE, selected=rrs)
	})
	regTypes <- reactive({
		input$rnbOptsI.region.types
	})
	regTypes.plus.sites <- reactive({c("sites", regTypes())})
	output$rnbOptsO.region.types <- renderText({
		rnb.options(region.types=regTypes())
		rnb.getOption("region.types")
	})
	output$selRegionProfiles.ex <- renderUI({
		selRegs <- optionSettingObserver$selRegs
		rrs <- optionSettingObserver$selRegs.exploratory.profiles
		if (is.null(rrs)) {
			rrs <- intersect(c("genes", "promoters", "cpgislands"), selRegs)
			optionSettingObserver$selRegs.exploratory.profiles <- rrs
		}
		selectInput('rnbOptsI.exploratory.region.profiles', NULL, selRegs, selected=rrs, multiple=TRUE)
	})
	output$selRegions.export <- renderUI({
		selRegs <- regTypes.plus.sites()
		rrs <- optionSettingObserver$selRegs.export
		if (is.null(rrs)) {
			rrs <- "sites"
			optionSettingObserver$selRegs.export <- rrs
		}
		selectInput('rnbOptsI.export.types', NULL, selRegs, selected=optionSettingObserver$selRegs.export, multiple=TRUE)
	})
	output$rnbOptsO.identifiers.column <- renderText({
		cname <- input$rnbOptsI.identifiers.column
		if (is.null(cname) || cname=="[None]") cname <- NULL
		rnb.options(identifiers.column=cname)
		rnb.getOption("identifiers.column")
	})
	output$rnbOptsO.min.group.size <- renderText({
		rnb.options(min.group.size=input$rnbOptsI.min.group.size)
		rnb.getOption("min.group.size")
	})
	observe({
		val.min.group.size <- input$rnbOptsI.min.group.size
		val.max.group.count <- input$rnbOptsI.max.group.count
		if (sannot.nsamples() > 0) optionSettingObserver$group.count.range.max <- sannot.nsamples()
		updateSliderInput(session, "rnbOptsI.min.group.size", value=val.min.group.size, min=RNB.GROUP.SIZE.RANGE[1], max=RNB.GROUP.SIZE.RANGE[2], step=1)
		updateSliderInput(session, "rnbOptsI.max.group.count", value=val.max.group.count, min=RNB.GROUP.COUNT.RANGE[1], max=optionSettingObserver$group.count.range.max, step=1)
	})
	output$rnbOptsO.max.group.count <- renderText({
		rnb.options(max.group.count=input$rnbOptsI.max.group.count)
		rnb.getOption("max.group.count")
	})
	output$rnbOptsO.colors.category <- renderText({
		cols <- optionSettingObserver$colors.category.list[[input$rnbOptsI.colors.category]]
		rnb.options(colors.category=cols)
		rnb.getOption("colors.category")
	})
	output$rnbOptsOP.colors.category <- renderPlot({
		cols <- optionSettingObserver$colors.category.list[[input$rnbOptsI.colors.category]]
		plotColPal(cols)
	})
	output$rnbOptsO.colors.meth <- renderText({
		cols <- RNB.COLSCHEMES.METH[[input$rnbOptsI.colors.meth]]
		rnb.options(colors.meth=cols)
		rnb.getOption("colors.meth")
	})
	output$rnbOptsOP.colors.meth <- renderPlot({
		cols <- RNB.COLSCHEMES.METH[[input$rnbOptsI.colors.meth]]
		# rnb.options(colors.meth=cols)
		# rnb.getOption("colors.meth")
		plotColPal(cols)
	})
	optsImportDataType <- reactive({
		res <- ""
		if (isBiseq()){
			res <- "bed.dir"
		} else {
			res <- "idat.dir"
		}
		rnb.options(import.default.data.type=res)
		res
	})
	output$rnbOptsO.import.default.data.type <- renderText({
		rnb.getOption("import.default.data.type")
	})
	output$rnbOptsO.import.table.separator <- renderText({
		rnb.options(import.table.separator=tabSep())
		rnb.getOption("import.table.separator")
	})
	output$rnbOptsO.import.bed.style <- renderText({
		interfaceSetting <- input$rnbOptsI.import.bed.style
		res <- rnb.getOption("import.bed.style")
		if (isBiseq()){
			res <- interfaceSetting
			rnb.options(import.bed.style=res)
		}
		res
	})
	doQc <- reactive({
		res <- input$rnbOptsI.qc
		oNames <- grep("^rnbOptsI.qc.", names(input), value=TRUE)
		if (res){
			for (oo in oNames){
				shinyjs::enable(oo)
			}
		} else {
			for (oo in oNames){
				shinyjs::disable(oo)
			}
		}
		res
	})
	output$rnbOptsO.qc <- renderText({
		rnb.options(qc=doQc())
		rnb.getOption("qc")
	})
	doPreprocessing <- reactive({
		res <- input$rnbOptsI.preprocessing
		oNames <- grep("^rnbOptsI.preprocessing.", names(input), value=TRUE)
		oNames <- union(oNames, grep("^rnbOptsI.normalization.", names(input), value=TRUE))
		oNames <- union(oNames, grep("^rnbOptsI.filtering.", names(input), value=TRUE))
		if (res){
			for (oo in oNames){
				shinyjs::enable(oo)
			}
		} else {
			for (oo in oNames){
				shinyjs::disable(oo)
			}
		}
		res
	})
	output$rnbOptsO.preprocessing <- renderText({
		rnb.options(preprocessing=doPreprocessing())
		rnb.getOption("preprocessing")
	})
	output$rnbOptsO.filtering.coverage.threshold <- renderText({
		interfaceSetting <- input$rnbOptsI.filtering.coverage.threshold
		res <- rnb.getOption("filtering.coverage.threshold")
		if (isBiseq()){
			res <- interfaceSetting
			rnb.options(filtering.coverage.threshold=res)
		}
		res
	})
	output$rnbOptsO.filtering.low.coverage.masking <- renderText({
		interfaceSetting <- input$rnbOptsI.filtering.low.coverage.masking
		res <- rnb.getOption("filtering.low.coverage.masking")
		if (isBiseq()){
			res <- interfaceSetting
			rnb.options(filtering.low.coverage.masking=res)
		}
		res
	})
	output$rnbOptsO.filtering.high.coverage.outliers <- renderText({
		interfaceSetting <- input$rnbOptsI.filtering.high.coverage.outliers
		res <- rnb.getOption("filtering.high.coverage.outliers")
		if (isBiseq()){
			res <- interfaceSetting
			rnb.options(filtering.high.coverage.outliers=res)
		}
		res
	})
	output$rnbOptsO.filtering.missing.value.quantile <- renderText({
		rnb.options(filtering.missing.value.quantile=input$rnbOptsI.filtering.missing.value.quantile)
		rnb.getOption("filtering.missing.value.quantile")
	})
	output$rnbOptsO.filtering.greedycut <- renderText({
		rnb.options(filtering.greedycut=input$rnbOptsI.filtering.greedycut)
		rnb.getOption("filtering.greedycut")
	})
	output$rnbOptsO.filtering.sex.chromosomes.removal <- renderText({
		rnb.options(filtering.sex.chromosomes.removal=input$rnbOptsI.filtering.sex.chromosomes.removal)
		rnb.getOption("filtering.sex.chromosomes.removal")
	})
	output$rnbOptsO.filtering.snp <- renderText({
		rnb.options(filtering.snp=input$rnbOptsI.filtering.snp)
		rnb.getOption("filtering.snp")
	})
	output$rnbOptsO.filtering.cross.reactive <- renderText({
		rnb.options(filtering.cross.reactive=input$rnbOptsI.filtering.cross.reactive)
		rnb.getOption("filtering.cross.reactive")
	})
	output$rnbOptsO.filtering.whitelist <- renderText({
		wl <- rnbOpts.filtering.whitelist.fn()
		if (is.null(wl) || nchar(wl) < 1 || wl=="[NULL]") wl <- NULL
		rnb.options(filtering.whitelist=wl)
		res <- rnb.getOption("filtering.whitelist")
		if (is.null(res)) res <- "NULL"
		res
	})
	output$rnbOptsO.filtering.blacklist <- renderText({
		bl <- rnbOpts.filtering.blacklist.fn()
		if (is.null(bl) || nchar(bl) < 1 || bl=="[NULL]") bl <- NULL
		rnb.options(filtering.blacklist=bl)
		res <- rnb.getOption("filtering.blacklist")
		if (is.null(res)) res <- "NULL"
		res
	})
	output$rnbOptsO.normalization.method <- renderText({
		interfaceSetting <- input$rnbOptsI.normalization.method
		res <- rnb.getOption("normalization.method")
		if (!isBiseq()){
			res <- interfaceSetting
			rnb.options(normalization.method=res)
		}
		res
	})
	output$rnbOptsO.normalization.background.method <- renderText({
		interfaceSetting <- input$rnbOptsI.normalization.background.method
		res <- rnb.getOption("normalization.background.method")
		if (!isBiseq()){
			res <- interfaceSetting
			rnb.options(normalization.background.method=res)
		}
		res
	})
	output$rnbOptsO.imputation.method <- renderText({
		rnb.options(imputation.method=input$rnbOptsI.imputation.method)
		rnb.getOption("imputation.method")
	})
	output$rnbOptsO.export.to.bed<- renderText({
		rnb.options(export.to.bed=input$rnbOptsI.export.to.bed)
		rnb.getOption("export.to.bed")
	})
	output$rnbOptsO.export.to.csv <- renderText({
		rnb.options(export.to.csv=input$rnbOptsI.export.to.csv)
		rnb.getOption("export.to.csv")
	})
	output$rnbOptsO.export.to.trackhub <- renderText({
		rnb.options(export.to.trackhub=input$rnbOptsI.export.to.trackhub)
		rnb.getOption("export.to.trackhub")
	})
	output$rnbOptsO.export.types <- renderText({
		rts <- input$rnbOptsI.export.types
		if (length(rts)<1) rts <- NULL
		rnb.options(export.types=rts)
		rnb.getOption("export.types")
	})
	doInference <- reactive({
		res <- input$rnbOptsI.inference
		res.agepred <- input$rnbOptsI.inference.age.prediction
		inferenceOptNames <- grep("^rnbOptsI\\.inference\\.", names(input), value=TRUE)
		inferenceOptNames.agepred <- setdiff(grep("^rnbOptsI\\.inference\\.age\\.", names(input), value=TRUE), "rnbOptsI.inference.age.prediction")
		if (res){
			for (oo in inferenceOptNames){
				if (is.element(oo, inferenceOptNames.agepred)){
					if (res.agepred){
						shinyjs::enable(oo)
					} else {
						shinyjs::disable(oo)
					}
				} else {
					shinyjs::enable(oo)
				}
			}
		} else {
			for (oo in inferenceOptNames){
				shinyjs::disable(oo)
			}
		}
		res
	})
	output$rnbOptsO.inference <- renderText({
		rnb.options(inference=doInference())
		rnb.getOption("inference")
	})
	output$rnbOptsO.inference.age.prediction <- renderText({
		rnb.options(inference.age.prediction=input$rnbOptsI.inference.age.prediction)
		rnb.getOption("inference.age.prediction")
	})
	output$rnbOptsO.inference.age.column <- renderText({
		cname <- input$rnbOptsI.inference.age.column
		if (is.null(cname) || cname=="[default]") cname <- "age"
		rnb.options(inference.age.column=cname)
		rnb.getOption("inference.age.column")
	})
	output$rnbOptsO.inference.targets.sva <- renderText({
		cnames <- input$rnbOptsI.inference.targets.sva
		if (length(cnames)<1) cnames <- character(0)
		rnb.options(inference.targets.sva=cnames)
		rnb.getOption("inference.targets.sva")
	})
	output$rnbOptsO.inference.sva.num.method <- renderText({
		rnb.options(inference.sva.num.method=input$rnbOptsI.inference.sva.num.method)
		rnb.getOption("inference.sva.num.method")
	})
	output$rnbOptsO.inference.reference.methylome.column <- renderText({
		cname <- input$rnbOptsI.inference.reference.methylome.column
		if (is.null(cname) || cname=="[None]") cname <- NULL
		rnb.options(inference.reference.methylome.column=cname)
		rnb.getOption("inference.reference.methylome.column")
	})
	doExploratory <- reactive({
		res <- input$rnbOptsI.exploratory
		oNames <- grep("^rnbOptsI.exploratory.", names(input), value=TRUE)
		if (res){
			for (oo in oNames){
				shinyjs::enable(oo)
			}
		} else {
			for (oo in oNames){
				shinyjs::disable(oo)
			}
		}
		res
	})
	output$rnbOptsO.exploratory <- renderText({
		rnb.options(exploratory=doExploratory())
		rnb.getOption("exploratory")
	})
	output$rnbOptsO.exploratory.columns <- renderText({
		cnames <- input$rnbOptsI.exploratory.columns
		if (length(cnames)<1) cnames <- NULL
		if (length(cnames)==1 && cnames=="[automatic]") cnames <- NULL
		cnames <- setdiff(cnames, "[automatic]")
		rnb.options(exploratory.columns=cnames)
		rnb.getOption("exploratory.columns")
	})
	output$rnbOptsO.exploratory.intersample <- renderText({
		rnb.options(exploratory.intersample=input$rnbOptsI.exploratory.intersample)
		rnb.getOption("exploratory.intersample")
	})
	output$rnbOptsO.exploratory.beta.distribution <- renderText({
		rnb.options(exploratory.beta.distribution=input$rnbOptsI.exploratory.beta.distribution)
		rnb.getOption("exploratory.beta.distribution")
	})
	output$rnbOptsO.exploratory.correlation.qc <- renderText({
		interfaceSetting <- input$rnbOptsI.exploratory.correlation.qc
		res <- rnb.getOption("exploratory.correlation.qc")
		if (!isBiseq()){
			res <- interfaceSetting
			rnb.options(exploratory.correlation.qc=res)
		}
		res
	})
	output$rnbOptsO.exploratory.region.profiles <- renderText({
		rts <- input$rnbOptsI.exploratory.region.profiles
		if (length(rts)<1) rts <- character(0)
		rnb.options(exploratory.region.profiles=rts)
		rnb.getOption("exploratory.region.profiles")
	})
	doDifferential <- reactive({
		res <- input$rnbOptsI.differential
		oNames <- grep("^rnbOptsI.differential.", names(input), value=TRUE)
		if (res){
			for (oo in oNames){
				shinyjs::enable(oo)
			}
		} else {
			for (oo in oNames){
				shinyjs::disable(oo)
			}
		}
		res
	})
	output$rnbOptsO.differential <- renderText({
		rnb.options(differential=doDifferential())
		rnb.getOption("differential")
	})
	output$rnbOptsO.differential.comparison.columns <- renderText({
		cnames <- input$rnbOptsI.differential.comparison.columns
		if (length(cnames)<1) cnames <- NULL
		if (length(cnames)==1 && cnames=="[automatic]") cnames <- NULL
		cnames <- setdiff(cnames, "[automatic]")
		rnb.options(differential.comparison.columns=cnames)
		rnb.getOption("differential.comparison.columns")
	})
	output$rnbOptsO.covariate.adjustment.columns <- renderText({
		cnames <- input$rnbOptsI.covariate.adjustment.columns
		if (length(cnames)<1) cnames <- NULL
		rnb.options(covariate.adjustment.columns=cnames)
		rnb.getOption("covariate.adjustment.columns")
	})
	output$rnbOptsO.differential.site.test.method <- renderText({
		rnb.options(differential.site.test.method=input$rnbOptsI.differential.site.test.method)
		rnb.getOption("differential.site.test.method")
	})
	output$rnbOptsO.differential.report.sites <- renderText({
		rnb.options(differential.report.sites=input$rnbOptsI.differential.report.sites)
		rnb.getOption("differential.report.sites")
	})
	doDiffVar <- reactive({
		res <- input$rnbOptsI.differential.variability
		if (res){
			shinyjs::enable("rnbOptsI.differential.variability.method")
		} else {
			shinyjs::disable("rnbOptsI.differential.variability.method")
		}
		res
	})
	output$rnbOptsO.differential.variability <- renderText({
		rnb.options(differential.variability=doDiffVar())
		rnb.getOption("differential.variability")
	})
	output$rnbOptsO.differential.variability.method <- renderText({
		rnb.options(differential.variability.method=input$rnbOptsI.differential.variability.method)
		rnb.getOption("differential.variability.method")
	})
	output$rnbOptsO.differential.enrichment.go <- renderText({
		rnb.options(differential.enrichment.go=input$rnbOptsI.differential.enrichment.go)
		rnb.getOption("differential.enrichment.go")
	})
	output$rnbOptsO.differential.enrichment.lola <- renderText({
		rnb.options(differential.enrichment.lola=input$rnbOptsI.differential.enrichment.lola)
		rnb.getOption("differential.enrichment.lola")
	})

	#apply the option setting 'ovalue' for option with name 'oname'
	applyOptValue <- function(oname, ovalue, fallback=FALSE){
		if (oname=="analysis.name") {
			updateTextInput(session, "rnbOptsI.analysis.name", value=ovalue)
		} else if (oname=="assembly") {
			if (is.element(ovalue, RNB.ASSEMBLIES)){
				updateSelectInput(session, "rnbOptsI.assembly", selected=ovalue)
			} else {
				stop(paste0("Invalid assembly: ", ovalue))
			}
		} else if (oname=="region.types") {
			if (is.null(ovalue) || all(ovalue %in% regTypes.all())){
				if (is.null(ovalue)) ovalue <- regTypes.all()
				optionSettingObserver$selRegs <- ovalue
				updateSelectInput(session, "rnbOptsI.region.types", selected=ovalue)
			} else {
				stop(paste0("Region type(s) not supported by current assembly (", rnb.getOption("assembly"), "): ", paste(setdiff(ovalue, regTypes.all()), collapse=", ")))
			}
		} else if (oname=="identifiers.column") {
			if (is.null(ovalue) || is.element(ovalue, sannot.cols.plusNone())){
				if (is.null(ovalue)) ovalue <- "[None]"
				optionSettingObserver$selColumn.id <- ovalue
				updateSelectInput(session, "rnbOptsI.identifiers.column", selected=ovalue)
			} else {
				stop(paste0("Sample annotation column not supported"))
			}
		} else if (oname=="min.group.size") {
			if (ovalue >= RNB.GROUP.SIZE.RANGE[1] && ovalue <= RNB.GROUP.SIZE.RANGE[2]){
				# print("DEBUG: updating min.group.size option")
				updateSliderInput(session, "rnbOptsI.min.group.size", value=ovalue)
			} else {
				stop(paste0("Not within expected range [", RNB.GROUP.SIZE.RANGE[1], "-", RNB.GROUP.SIZE.RANGE[2],"]: ", ovalue))
			}
		} else if (oname=="max.group.count") {
			if (is.null(ovalue) || (ovalue >= RNB.GROUP.COUNT.RANGE[1] && ovalue <= optionSettingObserver$group.count.range.max)){
				updateSliderInput(session, "rnbOptsI.max.group.count", value=ovalue)
			} else {
				stop(paste0("Not within expected range [", RNB.GROUP.COUNT.RANGE[1], "-", optionSettingObserver$group.count.range.max,"]: ", ovalue))
			}
		} else if (oname=="colors.category") {
			selVal <- ovalue
			if (length(ovalue)>1){
				selVal <- "[custom]"
				optionSettingObserver$colors.category.list[[selVal]] <- ovalue
			}
			updateSelectInput(session, "rnbOptsI.colors.category", choices=names(optionSettingObserver$colors.category.list), selected=selVal)
		} else if (oname=="import.default.data.type") {
			if (ovalue != optsImportDataType()){
				stop(paste0("Incompatible option with currently selected platform"))
			}
		} else if (oname=="import.table.separator") {
			if (is.element(ovalue, RNB.TABLE.SEPS)){
				updateSelectInput(session, "rnbOptsI.import.table.separator", selected=ovalue)
			} else {
				stop(paste0("Invalid table separator: ", ovalue))
			}
		} else if (oname=="import.bed.style") {
			if (is.element(ovalue, RNB.BED.STYLES)){
				updateSelectInput(session, "rnbOptsI.import.bed.style", selected=ovalue)
			} else {
				stop(paste0("Invalid import.bed.style: ", ovalue))
			}
		} else if (oname=="filtering.coverage.threshold") {
			if (ovalue >= 1 && ovalue <= 100){
				updateSliderInput(session, "rnbOptsI.filtering.coverage.threshold", value=ovalue)
			} else {
				stop(paste0("Not within expected range [", 1, "-", 100, "]: ", ovalue))
			}
		} else if (oname=="filtering.low.coverage.masking") {
			updateCheckboxInput(session, "rnbOptsI.filtering.low.coverage.masking", value=ovalue)
		} else if (oname=="filtering.high.coverage.outliers") {
			updateCheckboxInput(session, "rnbOptsI.filtering.high.coverage.outliers", value=ovalue)
		} else if (oname=="filtering.missing.value.quantile") {
			if (ovalue >= 0 && ovalue <= 1){
				updateSliderInput(session, "rnbOptsI.filtering.missing.value.quantile", value=ovalue)
			} else {
				stop(paste0("Not within expected range [", 0, "-", 1, "]: ", ovalue))
			}
		} else if (oname=="filtering.greedycut") {
			updateCheckboxInput(session, "rnbOptsI.filtering.greedycut", value=ovalue)
		} else if (oname=="filtering.sex.chromosomes.removal") {
			updateCheckboxInput(session, "rnbOptsI.filtering.sex.chromosomes.removal", value=ovalue)
		} else if (oname=="filtering.snp") {
			if (is.element(ovalue, RNB.FILTERING.SNP)){
				updateSelectInput(session, "rnbOptsI.filtering.snp", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="filtering.cross.reactive") {
			updateCheckboxInput(session, "rnbOptsI.filtering.cross.reactive", value=ovalue)
		} else if (oname=="normalization.method") {
			if (is.element(ovalue, RNB.NORMALIZATION.METHODS)){
				updateSelectInput(session, "rnbOptsI.normalization.method", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="normalization.background.method") {
			if (is.element(ovalue, RNB.NORMALIZATION.BG.METHODS)){
				updateSelectInput(session, "rnbOptsI.normalization.background.method", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="imputation.method") {
			if (is.element(ovalue, RNB.IMPUTATION.METHODS)){
				updateSelectInput(session, "rnbOptsI.imputation.method", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="export.to.trackhub") {
			if (is.null(ovalue) || all(ovalue %in% RNB.TRACKHUB.FORMATS)){
				if (is.null(ovalue)) ovalue <- character(0)
				updateSelectInput(session, "rnbOptsI.export.to.trackhub", selected=ovalue)
			} else {
				stop(paste0("Trackhub formats not supported: ", paste(setdiff(ovalue, RNB.TRACKHUB.FORMATS), collapse=", ")))
			}
		} else if (oname=="export.types") {
			# allowedVals <- regTypes.plus.sites()
			allowedVals <- c("sites", optionSettingObserver$selRegs)
			if (fallback) ovalue <- intersect(ovalue, allowedVals) #if resetting to old options, make sure that the old regions are contained in the allowed regions
			if (is.null(ovalue) || all(ovalue %in% allowedVals)){
				if (is.null(ovalue)) ovalue <- character(0)
				optionSettingObserver$selRegs.export <- ovalue
				updateSelectInput(session, "rnbOptsI.export.types", selected=ovalue)
			} else {
				stop(paste0("Export types not supported: ", paste(setdiff(ovalue, allowedVals), collapse=", ")))
			}
		} else if (oname=="inference") {
			updateCheckboxInput(session, "rnbOptsI.inference", value=ovalue)
		} else if (oname=="inference.age.prediction") {
			updateCheckboxInput(session, "rnbOptsI.inference.age.prediction", value=ovalue)
		} else if (oname=="inference.targets.sva") {
			print(ovalue)
			print(sannot.cols.grps())
			if (all(ovalue %in% sannot.cols.grps())){
				# print("DEBUG: updating inference.targets.sva option")
				optionSettingObserver$selColumn.sva <- ovalue
				updateSelectInput(session, "rnbOptsI.inference.targets.sva", selected=ovalue)
			} else {
				stop(paste0("Sample annotation column(s) not supported"))
			}
		} else if (oname=="inference.sva.num.method") {
			if (is.element(ovalue, RNB.SVA.NUM.METHODS)){
				updateSelectInput(session, "rnbOptsI.inference.sva.num.method", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="inference.reference.methylome.column") {
			if (is.null(ovalue) || is.element(ovalue, sannot.cols.plusNone())){
				if (is.null(ovalue)) ovalue <- "[None]"
				optionSettingObserver$selColumn.cellTypeRef <- ovalue
				updateSelectInput(session, "rnbOptsI.inference.reference.methylome.column", selected=ovalue)
			} else {
				stop(paste0("Sample annotation column not supported"))
			}
		} else if (oname=="exploratory.columns") {
			if (is.null(ovalue) || all(ovalue %in% sannot.cols.plusAutomatic())){
				if (is.null(ovalue)) ovalue <- "[automatic]"
				optionSettingObserver$selColumn.ex <- ovalue
				updateSelectInput(session, "rnbOptsI.exploratory.columns", selected=ovalue)
			} else {
				stop(paste0("Sample annotation column(s) not supported"))
			}
		} else if (oname=="exploratory.intersample") {
			updateCheckboxInput(session, "rnbOptsI.exploratory.intersample", value=ovalue)
		} else if (oname=="exploratory.beta.distribution") {
			updateCheckboxInput(session, "rnbOptsI.exploratory.beta.distribution", value=ovalue)
		} else if (oname=="exploratory.correlation.qc") {
			updateCheckboxInput(session, "rnbOptsI.exploratory.correlation.qc", value=ovalue)
		} else if (oname=="exploratory.region.profiles") {
			allowedVals <- optionSettingObserver$selRegs
			if (fallback) ovalue <- intersect(ovalue, allowedVals) #if resetting to old options, make sure that the old regions are contained in the allowed regions
			if (is.null(ovalue) || all(ovalue %in% allowedVals)){
				if (is.null(ovalue)) ovalue <- intersect(c("genes", "promoters", "cpgislands"), allowedVals)
				optionSettingObserver$selRegs.exploratory.profiles <- ovalue
				updateSelectInput(session, "rnbOptsI.exploratory.region.profiles", selected=ovalue)
			} else {
				stop(paste0("Region profiles not supported for region types: ", paste(setdiff(ovalue, allowedVals), collapse=", ")))
			}
		} else if (oname=="differential.comparison.columns") {
			if (is.null(ovalue) || all(ovalue %in% sannot.cols.plusAutomatic())){
				if (is.null(ovalue)) ovalue <- "[automatic]"
				optionSettingObserver$selColumn.diff <- ovalue
				updateSelectInput(session, "rnbOptsI.differential.comparison.columns", selected=ovalue)
			} else {
				stop(paste0("Sample annotation column(s) not supported"))
			}
		} else if (oname=="covariate.adjustment.columns") {
			if (is.null(ovalue) || all(ovalue %in% sannot.cols())){
				if (is.null(ovalue)) ovalue <- character(0)
				optionSettingObserver$selColumn.adj <- ovalue
				updateSelectInput(session, "rnbOptsI.covariate.adjustment.columns", selected=ovalue)
			} else {
				stop(paste0("Sample annotation column(s) not supported"))
			}
		} else if (oname=="differential.site.test.method") {
			if (is.element(ovalue, RNB.DIFFMETH.TEST.METHODS)){
				updateSelectInput(session, "rnbOptsI.differential.site.test.method", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="differential.report.sites") {
			updateCheckboxInput(session, "rnbOptsI.differential.report.sites", value=ovalue)
		} else if (oname=="differential.variability") {
			updateCheckboxInput(session, "rnbOptsI.differential.variability", value=ovalue)
		} else if (oname=="differential.variability.method") {
			if (is.element(ovalue, RNB.DIFFVAR.METHODS)){
				updateSelectInput(session, "rnbOptsI.differential.variability.method", selected=ovalue)
			} else {
				stop(paste0("Invalid selection: ", ovalue))
			}
		} else if (oname=="differential.enrichment.go") {
			updateCheckboxInput(session, "rnbOptsI.differential.enrichment.go", value=ovalue)
		} else if (oname=="differential.enrichment.lola") {
			updateCheckboxInput(session, "rnbOptsI.differential.enrichment.lola", value=ovalue)
		} 
	}

	#apply the options stored in list 'ol'
	applyOptList <- function(ol, ol.old=list()){
		for (oname in names(ol)){
			# print(paste("DEBUG: Reading XML option:", oname))
			rr <- tryCatch(
				applyOptValue(oname, ol[[oname]]),
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Could not update option '", oname, " (", err$message, ")")))
					if (is.element(oname, names(ol.old))){
						# print(paste("DEBUG: FAILED: Resetting to old option"))
						applyOptValue(oname, ol.old[[oname]], fallback=TRUE)
						optSettingList <- list(ol.old[[oname]])
						names(optSettingList) <- oname
						do.call("rnb.options", optSettingList)
					}
				}
			)
		}
		showNotification(tags$span(style="color:green", icon("check"), paste0("Option settings applied")))
	}

	observeLocalFileInput(input, session, 'loadOptsXmlFile')
	loadOptsXml.fName <- reactive({
		readLocalFileInput(session, 'loadOptsXmlFile')
	})
	observeEvent(input$loadOptsAnaDirDo, {
		curStatus <- anaStatus()
		isValid <- curStatus$status=="reportDir"
		xmlFile <- file.path(reportDir(), "analysis_options.xml")
		if (isValid && file.exists(xmlFile)){
			rnbOpts.old <- rnb.options()
			optList <- tryCatch({
					dummy <- rnb.xml2options(xmlFile)
					rnb.options()
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Could not load option file from analysis directory:", err$message)))
					NULL
				}
			)
			if (length(optList) > 0){
				applyOptList(optList, rnbOpts.old)
			}
		} else {
			showNotification(tags$span(style="color:red", icon("warning"), paste0("Option file does not exist in analysis directory")))
		}
	})
	observeEvent(input$loadOptsXmlDo, {
		# print("DEBUG: Reading XML file")
		xmlFile <- loadOptsXml.fName()
		if (file.exists(xmlFile)){
			rnbOpts.old <- rnb.options()
			# print("DEBUG: Old option settings")
			# print(rnbOpts.old)
			optList <- tryCatch({
					dummy <- rnb.xml2options(xmlFile)
					rnb.options()
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Could not load option file: ", err$message)))
					NULL
				}
			)
			if (length(optList) > 0){
				applyOptList(optList, rnbOpts.old)
			}
		} else {
			showNotification(tags$span(style="color:red", icon("warning"), paste0("Option file does not exist")))
		}
	})
	observeEvent(input$loadOptsProfileDo, {
		xmlFile <- file.path(RNB.OPTION.PROFILES.PATH, paste0(input$loadOptsProfileSel, ".xml"))
		if (file.exists(xmlFile)){
			rnbOpts.old <- rnb.options()
			optList <- tryCatch({
					dummy <- rnb.xml2options(xmlFile)
					rnb.options()
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Could not load option profile:", err$message)))
					NULL
				}
			)
			if (length(optList) > 0){
				applyOptList(optList, rnbOpts.old)
			}
		} else {
			showNotification(tags$span(style="color:red", icon("warning"), paste0("Option profile does not exist")))
		}
	})

	observeLocalFileInput(input, session, 'rnbOptsI.filtering.whitelist')
	rnbOpts.filtering.whitelist.fn <- reactive({
		readLocalFileInput(session, 'rnbOptsI.filtering.whitelist')
	})
	observeEvent(input$rnbOptsResetWhitelist, {
		updateLocalFileInput(session, 'rnbOptsI.filtering.whitelist', value='[NULL]')
	})
	observeLocalFileInput(input, session, 'rnbOptsI.filtering.blacklist')
	rnbOpts.filtering.blacklist.fn <- reactive({
		readLocalFileInput(session, 'rnbOptsI.filtering.blacklist')
	})
	observeEvent(input$rnbOptsResetblacklist, {
		updateLocalFileInput(session, 'rnbOptsI.filtering.blacklist', value='[NULL]')
	})

	output$saveOptsXml <- downloadHandler(
		filename="rnb_options.xml",
		content=function(fname){
			cat(rnb.options2xml(), file=fname)
		}
	)

	############################################################################
	# Run
	############################################################################
	enableRun <- reactive({
		res <- TRUE
		if (dir.exists(reportDir())) res <- FALSE #directory exists
		if (!file.exists(sampleAnnotFile())) res <- FALSE #sample annotation file does not exist
		if (!dir.exists(inputDataDir())) res <- FALSE # data directory does not exist
		return(res)
	})
	output$runHint <- renderUI({
		if (enableRun()){
			shinyjs::enable("runRnb")
			tags$p(tags$span(style="color:green", icon("play"), "Ready to run"))
		} else {
			shinyjs::disable("runRnb")
			tags$p(tags$span(style="color:red", icon("warning"), "Unable to run the analysis. Please make sure that you selected a non-existing report directory ('Analysis' tab) and that you specified the sample annotation file and data directory correctly ('Input' tab)"))
		}
	})
	observeEvent(input$ggplotTheme, {
		eval(parse(text=paste0("theme_set(theme_", input$ggplotTheme, "())")))
	})
	observeEvent(input$runRnb, {
		if(input$runRnb == 0) return()
        shinyjs::disable("runRnb")
    	setNumCoresFromSlider()

    	logger.close() # close the server's current logger s.t. a new log file is created in the report directory
		withProgress({
			rnb.run.analysis(
				dir.reports=reportDir(),
				sample.sheet=sampleAnnotFile(),
				data.dir=inputDataDir(),
				data.type=rnb.getOption("import.default.data.type"),
				save.rdata=TRUE
			)
		}, message="Runnning RnBeads analysis")
	})
	# output$runStatusMsg <- renderText({
	# 	refreshTimer()
	# 	print(input$runRnb)
	# 	if(input$runRnb == 0) return("Analysis not started")
	# 	res <- "Analysis started"
	# 	# logFn <- file.path(reportDir(), "analysis.log")
	# 	# status <- getRnbStatusFromLog(logFn)
	# 	# status <- status[status$scheduled,]
	# 	# if (all(status$status=="completed")) return("Analysis completed")
	# 	# startedModules <- which(status$status=="started")
	# 	# if (length(startedModules) > 0){
	# 	# 	mm <- status$module[startedModules[length(startedModules)]] # pick last started module
	# 	# 	res <- paste0("Running ", mm)
	# 	# }
	# 	return(res)
	# })

	############################################################################
	# Modules
	############################################################################
	# IMPORT
	output$modImportNew.about <- renderUI({
		if (enableRun()){
			if (!modImportStatus$dataset.loaded){
				shinyjs::enable("modImportNew")
				shinyjs::enable("modImportNew.save")
			}
			tagList(
				tags$p(tags$span(icon("play"), "Ready to run")),
				checkboxInput("modImportNew.save", "Save RnBSet object", value=TRUE)
			)
		} else {
			shinyjs::disable("modImportNew")
			shinyjs::disable("modImportNew.save")
			tags$p(tags$span(icon("warning"), "Unable to run. Please make sure that you selected a non-existing report directory ('Analysis' tab) and that you specified the sample annotation file and data directory correctly ('Input' tab)"))
		}
	})
	output$modImportAnaDir.about <- renderUI({
		curStatus <- anaStatus()
		if (curStatus$status=="reportDir"){
			if (length(curStatus$rnbSet.paths) > 0){
				pp <- curStatus$rnbSet.paths
				names(pp) <- basename(pp)
				tagList(
					tags$p("An existing RnBeads report has been specified. Ready to import."),
					selectInput("modImportAnaDir.rnbSet", "Select RnBSet to load", pp)
				)
			} else {
				tagList(tags$p(icon("warning"), "The analysis directory does not contain RnBSet objects."))
			}
		} else {
			tagList(tags$p(icon("warning"), "Please specify an existing RnBeads report in the 'Analysis' tab to continue."))
		}
	})
	output$rnbSetInfo <- renderUI({
		if (is.null(rnbData$rnbSet)) return(NULL)
		tagList(
			tags$h2("RnBSet Object:"),
			renderPrint({methods::show(rnbData$rnbSet)}),
			tags$h2("Sample Annotation:"),
			renderTable({sannot()}, striped=TRUE, hover=TRUE, bordered=TRUE)
		)
	})
	modImportStatus <- reactiveValues(
		dataset.loaded=FALSE,
		dataset.loaded.nsamples=-1
	)
	rnbData <- reactiveValues(
		rnbSet=NULL
	)
	# watch if dataset has been loaded
	reactDataset <- observeEvent(modImportStatus$dataset.loaded, {
		if (modImportStatus$dataset.loaded){
			shinyjs::disable("modImportNew")
			shinyjs::disable("modImportNew.save")
			shinyjs::disable("modImportAnaDir")
			shinyjs::disable("modImportRObjects")
			shinyjs::enable("modImportReset")
		} else {
			if (enableRun()) {
				shinyjs::enable("modImportNew")
				shinyjs::enable("modImportNew.save")
			}
			if (anaStatus()$status=="reportDir" && length(anaStatus()$rnbSet.paths) > 0) shinyjs::enable("modImportAnaDir")
			shinyjs::enable("modImportRObjects")
			shinyjs::disable("modImportReset")
		}
	})
	output$modImport.status <- renderUI({
		# curStatus <- modImportStatus()
		res <- list()
		if (modImportStatus$dataset.loaded){
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Loaded dataset containing", modImportStatus$dataset.loaded.nsamples, "samples"))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded"))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		tagList(res)
	})
	observeEvent(input$modImportNew, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.import(c(inputDataDir(), sampleAnnotFile()), data.type=rnb.getOption("import.default.data.type"), dir.reports=reportDir())
					logger.close()
					rnbData$rnbSet <- res$rnb.set
					markDirDJ(reportDir())
					if (input$modImportNew.save) save.rnb.set(rnbData$rnbSet, file.path(reportDir(), "rnbSet_import"), archive=FALSE)
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbData$rnbSet))
				},
				error = function(err) {
					rnbData$rnbSet <- NULL
					modImportStatus$dataset.loaded <- FALSE
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Failed to import data: ", err$message)))
				}
			)
		}, message="Importing dataset")
	})
	observeEvent(input$modImportAnaDir, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					rnbData$rnbSet <- load.rnb.set(input$modImportAnaDir.rnbSet)
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbData$rnbSet))
				},
				error = function(err) {
					rnbData$rnbSet <- NULL
					modImportStatus$dataset.loaded <- FALSE
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Failed to load RnBSet object: ", err$message)))
				}
			)
		}, message="Loading dataset")
		#TODO: import options
	})
	observeDirectoryInput(input, session, 'modImportRnBSetDir')
	# observeLocalFileInput(input, session, 'modImportOptionsFile')
	observeEvent(input$modImportRObjects, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					rnbData$rnbSet <- load.rnb.set(readDirectoryInput(session, 'modImportRnBSetDir'))
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbData$rnbSet))
				},
				error = function(err) {
					rnbData$rnbSet <- NULL
					modImportStatus$dataset.loaded <- FALSE
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Failed to load RnBSet object: ", err$message)))
				}
			)
		}, message="Loading dataset")
		#TODO: import options
	})
	observeEvent(input$modImportReset, {
		modImportStatus$dataset.loaded <- FALSE
		modImportStatus$dataset.loaded.nsamples <- -1
		rnbData$rnbSet <- NULL
	})

	# QUALITY CONTROL
	reportExists.quality_control <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["quality_control", "Report"])
	})
	output$modQC.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.quality_control()){
			if (rdy){
				shinyjs::enable("modQC.overwrite")
			}
			if (!input$modQC.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["quality_control"])), "View report"))))))
		} else {
			shinyjs::disable("modQC.overwrite")
		}
		if (rdy){
			shinyjs::enable("modQC.run")
		} else {
			shinyjs::disable("modQC.run")
		}
		tagList(res)
	})
	observeEvent(input$modQC.run, {
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.quality_control() && input$modQC.overwrite){
						unlink(file.path(reportDir(), paste0("quality_control","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modQC.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.qc(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Quality Control) completed")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Quality Control) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Quality Control")
	})

	# PREPROCESSING
	reportExists.preprocessing <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["preprocessing", "Report"])
	})
	output$modPreprocessing.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.preprocessing()){
			if (rdy){
				shinyjs::enable("modPreprocessing.overwrite")
			}
			if (!input$modPreprocessing.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["preprocessing"])), "View report"))))))
		} else {
			shinyjs::disable("modPreprocessing.overwrite")
		}
		if (rdy){
			shinyjs::enable("modPreprocessing.run")
		} else {
			shinyjs::disable("modPreprocessing.run")
		}
		tagList(res)
	})
	observeEvent(input$modPreprocessing.run, {
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.preprocessing() && input$modPreprocessing.overwrite){
						unlink(file.path(reportDir(), paste0("preprocessing","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modPreprocessing.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.preprocessing(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					rnbData$rnbSet <- res$rnb.set
					if (input$modPreprocessing.save) {
						rnbsPath <- file.path(reportDir(), "rnbSet_preprocessed")
						if (dir.exists(rnbsPath)){
							unlink(rnbsPath, recursive=TRUE)
							showNotification(tags$span(icon("trash"), paste0("Deleted previous RnBSet object from report directory")))
						}
						save.rnb.set(rnbData$rnbSet, rnbsPath, archive=FALSE)
					}
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Preprocessing) completed")))
					showNotification(tags$span(style="color:green", icon("warning"), paste0("Replaced loaded RnBSet object with preprocessed object")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Preprocessing) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Preprocessing")
	})

	# TRACKS AND TABLES
	reportExists.tracks_and_tables <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["tracks_and_tables", "Report"])
	})
	output$modTNT.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.tracks_and_tables()){
			if (rdy){
				shinyjs::enable("modTNT.overwrite")
			}
			if (!input$modTNT.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["tracks_and_tables"])), "View report"))))))
		} else {
			shinyjs::disable("modTNT.overwrite")
		}
		if (rdy){
			shinyjs::enable("modTNT.run")
		} else {
			shinyjs::disable("modTNT.run")
		}
		tagList(res)
	})
	observeEvent(input$modTNT.run, {
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.tracks_and_tables() && input$modTNT.overwrite){
						unlink(file.path(reportDir(), paste0("tracks_and_tables","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modTNT.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.tnt(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Tracks and Tables) completed")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Tracks and Tables) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Tracks and Tables")
	})

	# COVARIATE INFERENCE
	reportExists.covariate_inference <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["covariate_inference", "Report"])
	})
	output$modInference.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.covariate_inference()){
			if (rdy){
				shinyjs::enable("modInference.overwrite")
			}
			if (!input$modInference.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["covariate_inference"])), "View report"))))))
		} else {
			shinyjs::disable("modInference.overwrite")
		}
		if (rdy){
			shinyjs::enable("modInference.run")
		} else {
			shinyjs::disable("modInference.run")
		}
		tagList(res)
	})
	observeEvent(input$modInference.run, {
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.covariate_inference() && input$modInference.overwrite){
						unlink(file.path(reportDir(), paste0("covariate_inference","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modInference.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.inference(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					rnbData$rnbSet <- res$rnb.set
					if (input$modInference.save) {
						rnbsPath <- file.path(reportDir(), "rnbSet_inference")
						if (dir.exists(rnbsPath)){
							unlink(rnbsPath, recursive=TRUE)
							showNotification(tags$span(icon("trash"), paste0("Deleted previous RnBSet object from report directory")))
						}
						save.rnb.set(rnbData$rnbSet, rnbsPath, archive=FALSE)
					}
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Covariate Inference) completed")))
					showNotification(tags$span(style="color:green", icon("warning"), paste0("Replaced loaded RnBSet object with inference object")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Covariate Inference) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Covariate Inference")
	})

	# EXPLORATORY ANALYSIS
	reportExists.exploratory_analysis <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["exploratory_analysis", "Report"])
	})
	output$modExploratory.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.exploratory_analysis()){
			if (rdy){
				shinyjs::enable("modExploratory.overwrite")
			}
			if (!input$modExploratory.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["exploratory_analysis"])), "View report"))))))
		} else {
			shinyjs::disable("modExploratory.overwrite")
		}
		if (rdy){
			shinyjs::enable("modExploratory.run")
		} else {
			shinyjs::disable("modExploratory.run")
		}
		tagList(res)
	})
	observeEvent(input$modExploratory.run, {
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.exploratory_analysis() && input$modExploratory.overwrite){
						unlink(file.path(reportDir(), paste0("exploratory_analysis","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modExploratory.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.exploratory(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Explorarory Analysis) completed")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Explorarory Analysis) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Explorarory Analysis")
	})

	# DIFFERENTIAL METHYLATION
	reportExists.differential_methylation <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["differential_methylation", "Report"])
	})
	output$modDifferential.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.differential_methylation()){
			if (rdy){
				shinyjs::enable("modDifferential.overwrite")
			}
			if (!input$modDifferential.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["differential_methylation"])), "View report"))))))
		} else {
			shinyjs::disable("modDifferential.overwrite")
		}
		if (rdy){
			shinyjs::enable("modDifferential.run")
		} else {
			shinyjs::disable("modDifferential.run")
		}
		tagList(res)
	})
	observeEvent(input$modDifferential.run, {
		withProgress({
			tryCatch({
					setNumCoresFromSlider()
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.differential_methylation() && input$modDifferential.overwrite){
						unlink(file.path(reportDir(), paste0("differential_methylation","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modDifferential.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.differential(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Differential Methylation) completed")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Differential Methylation) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Differential Methylation")
	})
}

################################################################################
# Main
################################################################################
shinyApp(ui = ui, server = server)

################################################################################
# Sandbox
################################################################################
# useful link for themes/layout
# https://shiny.rstudio.com/gallery/shiny-theme-selector.html

################################################################################
# TODOs:
# - sannot.cols.grps updates last during XML file loading --> SVA column option has to be loaded twice
################################################################################

