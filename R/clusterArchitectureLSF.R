################################################################################
# Cluster Architecture Descriptions
################################################################################
################################################################################
# Concrete implementations for the SLURM environment
################################################################################

#' ClusterArchitectureSLURM Class
#'
#' A child class of \code{\linkS4class{ClusterArchitecture}} implementing specifications of Simple Linux Utility for Resource Management (SLURM) architectures.
#'
#' @details
#' Follow this template if you want to create your own ClusterArchitecture class.
#'
#' @section Slots:
#' see \code{\linkS4class{ClusterArchitecture}}
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getSubCmdTokens,ClusterArchitectureSGE-method}}}{Returns a vector of command line tokens corresponding to submitting
#'   a job with the given command to the cluster}
#' }
#'
#' @name ClusterArchitectureSLURM-class
#' @rdname ClusterArchitectureSLURM-class
#' @author Michael Scherer
#' @exportClass ClusterArchitecture
setClass("ClusterArchitectureSLURM",
	contains = "ClusterArchitecture"
)

#' initialize.ClusterArchitectureSLURM
#'
#' Initialize an ClusterArchitecture object for a SLURM
#' 
#' @param .Object New instance of \code{ClusterArchitectureSLURM}.
#' @param name A name or identifier
#' @param ... arguments passed on to the constructor of \code{\linkS4class{ClusterArchitecture}} (the parent class)
#'
#' @export
#' @author Michael Scherer
#' @docType methods
setMethod("initialize","ClusterArchitectureSLURM",
	function(
		.Object,
		name="ClusterArchitectureSLURM",
		...
	) {
		.Object <- callNextMethod(.Object=.Object, name=name, ...)
		.Object <- setExecutable(.Object,"R","R")
		.Object <- setExecutable(.Object,"Rscript","Rscript")
		.Object <- setExecutable(.Object,"python","python")
		.Object@getSubCmdTokens.optional.args <- c("sub.binary","quote.cmd")
		.Object
	}
)

#' getSubCmdTokens-methods
#'
#' Returns a string for the of command line corresponding to submitting
#' a job with the given command to the cluster.
#' @details
#' For a concrete child class implementation for a SLURM architecture specification see \code{\linkS4class{ClusterArchitectureSLURM}}
#'
#' @param object \code{\linkS4class{ClusterArchitectureSLURM}} object
#' @param cmd.tokens a character vector specifying the executable command that should be wrapped in the cluster submission command
#' @param log file name and path of the log file that the submitted job writes to
#' @param job.name name of the submitted job
#' @param res.req named vector of requested resources. Two options are available: \code{"clock.limit"} and \code{"memory.size"}
#' @param sub.binary flag indicating if the command is to be submitted using the \code{"wrap"} option of SLURM
#' @param depend.jobs character vector containg names or ids of jobs the submitted job will depend on.
#' @param quote.cmd Flag indicating whether the submitted cammed should also be wrapped in quotes
#' @return A character vector containing the submission command tokens
#'
#' @rdname getSubCmdTokens-ClusterArchitectureSLURM-methods
#' @docType methods
#' @aliases getSubCmdTokens,ClusterArchitectureSLURM-method
#' @author Michael Scherer
#' @export
#' @examples
#' \donttest{
#' arch <- new("ClusterArchitectureSLURM",
#' 	name="my_slurm_architecture"
#' )
#' getSubCmdTokens(arch,c("Rscript","my_great_script.R"),"my_logfile.log")
#' }
setMethod("getSubCmdTokens",
	signature(
		object="ClusterArchitectureSLURM"
	),
	function(
	  object,
	  cmd.tokens,
	  log,
	  job.name = "",
	  res.req = character(0),
	  depend.jobs = character(0),
	  sub.binary = TRUE,
	  quote.cmd = TRUE
	) {
	  res.req.token <- NULL
		if(length(res.req)>0){
		  if("clock.limit" %in% names(res.req)){
		    res.req.token <- paste(res.req.token,"-t",res.req["clock.limit"]," ",collapse = "")
		  }
		  if("mem.size" %in% names(res.req)){
		    res.req.token <- paste0(res.req.token,"--mem=",res.req["mem.size"],collapse="")
		    
		  }
		}
		log.token <- NULL
		if (nchar(log)>0) {
			log.token <- c("-o",log)
		}
		job.name.token <- NULL
		if (nchar(job.name)>0) {
			job.name.token <- paste0(job.name.token,"--job-name=",job.name,collapse = "")
		}
		dependency.token <- NULL
		if (length(depend.jobs)>0){
		  get.job.id <- function(x){
		    cmd <- paste0("echo $(squeue --noheader --format %i --name ",x,")")
		    system(cmd,intern = T)
		  }
		  depend.jobs <- sapply(depend.jobs,get.job.id)
			dependency.token <- paste0(dependency.token, "--depend=", paste0(paste(depend.jobs,collapse=",")),collapse = "")
		}
		wrap.token <- NULL
		if(sub.binary){
		  if (quote.cmd){
		    cmd.tokens <- paste0("--wrap=",paste0("'",paste(cmd.tokens,collapse=" "),"'"),collapse="")
		  }else{
		    cmd.tokens <- paste0("--wrap=",paste(cmd.tokens,collapse=" "),collapse="")
		  }
		}else{
		  if (quote.cmd){
		    cmd.tokens <- paste0("'",paste(cmd.tokens,collapse=" "),"'")
		  }else{
		    cmd.tokens <- paste(cmd.tokens,collapse=" ")
		  }
		}

		res <- c(
			"sbatch",
			"--export=ALL",
			res.req.token,
			log.token,
			job.name.token,
			dependency.token,
			wrap.token,
			cmd.tokens
		)
		return(res)
	}
)
