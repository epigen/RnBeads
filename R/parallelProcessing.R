########################################################################################################################
## parallelProcessing.R
## created: 2012-05-04
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for parallel processing.
########################################################################################################################

## G L O B A L S #######################################################################################################

.parallel <- new.env()
.parallel[["do.par"]] <- FALSE
.parallel[["num.cores"]] <- -1L

## F U N C T I O N S ###################################################################################################

#' parallel.setup
#'
#' Sets up parallel processing. Requires the \pkg{foreach} and \pkg{doParallel} packages
#' @param ... Parameters for \code{registerDoParallel} from the \pkg{doParallel} package.
#' 			  This allows, for instance, for specificying the number of workers.
#' @return \code{TRUE} (invisible) to indicate that parallelization is set up.
#' @note Requires the packages \pkg{foreach} and \pkg{doParallel}.
#'
#' @author Fabian Mueller
#' @export parallel.setup
#' @examples
#' \donttest{
#' parallel.setup(2)
#' parallel.teardown()
#' }
parallel.setup <- function(...){
	logger.start("Setting up Multicore")
	require(foreach)
	require(doParallel)
	# .parallel[["cl"]] <- makeCluster(...)
	# registerDoParallel(.parallel[["cl"]])
	registerDoParallel(...)
	.parallel[["num.cores"]] <- getDoParWorkers()
	.parallel[["do.par"]] <- TRUE
	logger.info(c("Using",.parallel[["num.cores"]],"cores"))
	logger.completed()
	invisible(.parallel[["do.par"]])
}

########################################################################################################################

#' parallel.teardown
#'
#' Disables parallel processing.
#'
#' @return \code{TRUE}, invisibly.
#'
#' @author Fabian Mueller
#' @export parallel.teardown
#' @examples
#' \donttest{
#' parallel.getNumWorkers()
#' parallel.setup(2)
#' parallel.getNumWorkers()
#' parallel.teardown()
#' parallel.getNumWorkers()
#' }
parallel.teardown <- function(){
	.parallel[["do.par"]] <- FALSE
	.parallel[["num.cores"]] <- -1
	# stopCluster(.parallel[["cl"]])
	.stopImplicitCluster()
	invisible(TRUE)
}

########################################################################################################################

#' parallel.getNumWorkers
#'
#' Gets the number of workers used for parallel processing.
#' 
#' @return Number of workers used for parallel processing; \code{-1} if parallel processing is not enabled.
#'
#' @author Fabian Mueller
#' @export parallel.getNumWorkers
#' @examples
#' \donttest{
#' parallel.getNumWorkers()
#' parallel.setup(2)
#' parallel.getNumWorkers()
#' parallel.teardown()
#' parallel.getNumWorkers()
#' }
parallel.getNumWorkers <- function(){
	return(.parallel[["num.cores"]])
}

########################################################################################################################

#' parallel.isEnabled
#'
#' Checks if whether parallel processing is enabled.
#'
#' @return \code{TRUE} if multicore processing is enabled, \code{FALSE} otherwise.
#'
#' @author Fabian Mueller
#' @export parallel.isEnabled
#' @examples
#' \donttest{
#' parallel.isEnabled()
#' parallel.setup(2)
#' parallel.isEnabled()
#' parallel.teardown()
#' parallel.isEnabled()
#' }
parallel.isEnabled <- function(){
	return(.parallel[["do.par"]])
}
