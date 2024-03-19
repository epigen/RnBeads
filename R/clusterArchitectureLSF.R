################################################################################
# Cluster Architecture Descriptions
################################################################################
################################################################################
# Concrete implementations for the LSF environment
################################################################################

#' ClusterArchitectureLSF Class
#'
#' A child class of \code{\linkS4class{ClusterArchitecture}} implementing specifications of IBM LSF architectures.
#'
#' @details
#' Follow this template if you want to create your own ClusterArchitecture class.
#'
#' @section Slots:
#' see \code{\linkS4class{ClusterArchitecture}}
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getSubCmdTokens,ClusterArchitectureLSF-method}}}{Returns a vector of command line tokens corresponding to submitting
#'   a job with the given command to the cluster}
#' }
#'
#' @name ClusterArchitectureLSF-class
#' @rdname ClusterArchitectureLSF-class
#' @author Michael Scherer
#' @exportClass ClusterArchitecture
setClass("ClusterArchitectureLSF",
	contains = "ClusterArchitecture"
)

#' initialize.ClusterArchitectureLSF
#'
#' Initialize an ClusterArchitecture object for a LSF
#' 
#' @param .Object New instance of \code{ClusterArchitectureLSF}.
#' @param name A name or identifier
#' @param ... arguments passed on to the constructor of \code{\linkS4class{ClusterArchitecture}} (the parent class)
#'
#' @export
#' @author Michael Scherer
#' @docType methods
setMethod("initialize","ClusterArchitectureLSF",
	function(
		.Object,
		name="ClusterArchitectureLSF",
		...
	) {
		.Object <- callNextMethod(.Object=.Object, name=name, ...)
		.Object <- setExecutable(.Object,"R","R")
		.Object <- setExecutable(.Object,"Rscript","Rscript")
		.Object <- setExecutable(.Object,"python","python")
		.Object
	}
)

#' getSubCmdTokens-methods
#'
#' Returns a string for the of command line corresponding to submitting
#' a job with the given command to the cluster.
#' @details
#' For a concrete child class implementation for a LSF architecture specification see \code{\linkS4class{ClusterArchitectureLSF}}
#'
#' @param object \code{\linkS4class{ClusterArchitectureLSF}} object
#' @param cmd.tokens a character vector specifying the executable command that should be wrapped in the cluster submission command
#' @param log file name and path of the log file that the submitted job writes to
#' @param job.name name of the submitted job
#' @param res.req named vector of requested resources. Two options are available: \code{"clock.limit"} and \code{"memory.size"}
#' @param depend.jobs character vector containg names or ids of jobs the submitted job will depend on.
#' @return A character vector containing the submission command tokens
#'
#' @rdname getSubCmdTokens-ClusterArchitectureLSF-methods
#' @docType methods
#' @aliases getSubCmdTokens,ClusterArchitectureLSF-method
#' @author Michael Scherer
#' @export
#' @examples
#' \donttest{
#' arch <- new("ClusterArchitectureLSF",
#' 	name="my_lsf_architecture"
#' )
#' getSubCmdTokens(arch,c("Rscript","my_great_script.R"),"my_logfile.log")
#' }
setMethod("getSubCmdTokens",
	signature(
		object="ClusterArchitectureLSF"
	),
	function(
	  object,
	  cmd.tokens,
	  log,
	  job.name = "",
	  res.req = character(0),
	  depend.jobs = character(0)
	) {
	  res.req.token <- NULL
	  queue <- "short"
		if(length(res.req)>0){
		  if("clock.limit" %in% names(res.req)){
		    res.req.token <- paste(res.req.token,"-W ",res.req["clock.limit"]," ",collapse = "")
		  }
		  if("mem.size" %in% names(res.req)){
		    res.req.token <- paste0(res.req.token,'-R "rusage[mem=', res.req["mem.size"], ']"', collapse="")
		  }
		  if("queue" %in% names(res.req)){
		    queue <-  res.req["queue"]
		  }
		}
		log.token <- NULL
		if (nchar(log)>0) {
			log.token <- c("-e", log, "-o", log)
		}
		job.name.token <- NULL
		if (nchar(job.name)>0) {
			job.name.token <- paste0(job.name.token,"-J ",job.name,collapse = "")
		}
		dependency.token <- NULL
		if (length(depend.jobs)>0){
			dependency.token <- paste0(dependency.token, "-w ", paste0('"', paste(paste0('done(', depend.jobs, ')'),collapse="&&"), '"'),collapse = "")
		}

		res <- c(
			"bsub",
			paste("-q", queue),
			res.req.token,
			log.token,
			job.name.token,
			dependency.token,
			cmd.tokens
		)
		print(paste(res, collapse=" "))
		return(res)
	}
)
