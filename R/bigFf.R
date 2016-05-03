################################################################################
# Wrapper for ff objects with more than INT_MAX slots
################################################################################
setClassUnion("characterOrNULL", c("character", "NULL"))

#' BigFfMat Class
#'
#' A class implementing matrices as lists of ff objects by column
#' 
#' @details
#' useful to get rid of the INT_MAX limitatation of ff which is easily reached when specifying large matrices
#'
#' @section Slots:
#' \describe{
#'   \item{\code{cols}}{List of ff objects to store the columns of the matrix}
#'   \item{\code{colNames}}{vector of column names or NULL}
#'   \item{\code{rowNames}}{vector of row names or NULL} 
#'   \item{\code{colN}}{number of matrix columns} 
#'   \item{\code{rowN}}{number of matrix rows} 
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{show}}{Display the object}
#'   \item{\code{nrow}}{get the number of rows in the matrix}
#'   \item{\code{ncol}}{get the number of columns in the matrix}
#'   \item{\code{dim}}{vector of the number of rows and columns in the matrix}
#'   \item{\code{rownames}}{get row names}
#'   \item{\code{rownames<-}}{assign row names}
#'   \item{\code{colnames}}{get column names}
#'   \item{\code{colnames<-}}{assign column names}
#'   \item{\code{[}}{array indexing of the stored matrix}
#'   \item{\code{[<-}}{assign values by array indexing in the stored matrix}
#' }
#'
#' @name BigFfMat-class
#' @noRd
#' @author Fabian Mueller
#' @examples
#' myDf <- data.frame(normal=rnorm(10), integers=sample(1:100,10), index=1:10)
#' tt <- BigFfMat(row.n=4, col.n=3, col.names=LETTERS[1:3])
#' ttr <- BigFfMat(myDf)
#' nrow(tt)
#' ncol(tt)
#' dim(tt)
#' colnames(tt)
#' colnames(tt) <- c("blabb","bli","blubb")
#' colnames(tt)[1:2] <- c("bla","padoink")
#' rownames(tt) <- paste0("r",1:nrow(tt))
#' ttr[6,2]
#' ttr[5:7,c(TRUE,TRUE,FALSE)]
#' ttr[,2]
#' ttr[1:2,]
#' rownames(ttr) <- paste0("r",1:nrow(ttr))
#' ttr[c("r5","r7"),c("integers","normal")]
#' tt[4,2] <- 555
#' tt[1,] <- 1:3
#' tt[2,] <- 6
#' tt[3:4,2:3] <- matrix(10,nrow=2,ncol=2)
#' tt[c("r2","r3"),c("padoink","blubb")] <- -8
#' tt[3:4] <- -666 # currently assigns the value to the entire rows instead of using matrix indices
#' # tt[3:5] <- 666 #does not work
#' tt[,1][3:4] <- -1
#' tt[2:3,][,2] <- c(77,777)
#' tt[2:3,2:3][1:4] <- 1:4
#' mask <- matrix(c(TRUE,FALSE,FALSE,TRUE),ncol=2)
#' tt[2:3,2:3][mask] <- c(1000,9999)
#' tt
#' ttr
setClass("BigFfMat",
	slots=list(cols="list", colNames="characterOrNULL", rowNames="characterOrNULL", colN="integer", rowN="integer", activeClosing="logical")
)

#' Construct BigFfMat objects
#' 
#' @param df       	a matrix or data.frame to be stored
#' @param row.n     specify the number of rows. only plays a role if \code{df} is missing.
#' @param col.n     specify the number of columns. only plays a role if \code{df} is missing.
#' @param row.names specify row names. only plays a role if \code{df} is missing.
#' @param col.names specify column names. only plays a role if \code{df} is missing.
#' @param na.prototype type of NA object that is used for filling empty matrices. only plays a role if \code{df} is missing.
#' @param ...		paramters for \code{ff}
#' 
#' @return an object of class BigFfMat
#' 
#' @name BigFfMat
#' @noRd
BigFfMat <- function(df, row.n, col.n, row.names=NULL, col.names=NULL, na.prototype=as.numeric(NA), active.closing=TRUE, ...){
	#initialize empty matrix if the df object is missing
	if (missing(df)){
		if (missing(row.n) | missing(col.n)){
			stop("Must specify row.n and col.n")
		}
		# df <- matrix(na.prototype, nrow=row.n, ncol=col.n)
		ffList <- lapply(1:col.n, FUN=function(j){
			res <- ff(rep(na.prototype, row.n), ...)
			if (active.closing) close.ff(res)
			res
		})
		row.n <- as.integer(row.n)
		col.n <- as.integer(col.n)
		if (!is.null(row.names)){
			if (is.factor(row.names)){
				row.names <- as.character(row.names)
			}
		}
		if (!is.null(col.names)){
			if (is.factor(col.names)){
				col.names <- as.character(col.names)
			}
		}
	} else {	
		if (is.data.frame(df)){
			df <- as.matrix(df)
		}
		if (is.character(df)){
			stop("Character matrices currently not supported")
		}
		row.n <- nrow(df)
		col.n <- ncol(df)
		ffList <- lapply(1:col.n, FUN=function(j){
			res <- ff(df[,j], ...)
			if (active.closing) close.ff(res)
			res
		})
		row.names <- rownames(df)
		col.names <- colnames(df)
	}
	res <- new("BigFfMat", cols=ffList, colNames=col.names, rowNames=row.names, colN=col.n, rowN=row.n, activeClosing=active.closing)
}

setMethod("show", "BigFfMat",
	function(object){
		cat("BigFfDf object of dimension ", object@rowN, " x ", object@colN, "\n", sep="")
		if (!is.null(object@rowNames)){
			if (object@rowN > 10) {
				rowNstr <- paste(paste(rep(object@rowNames, length.out=10), collapse=", "),", ...")
			} else {
				rowNstr <- paste(object@rowNames, collapse=", ")
			}
			cat("  row names: ", rowNstr, "\n", sep="")
		}
		if (!is.null(object@colNames)){
			if (object@colN > 10) {
				colNstr <- paste(paste(rep(object@colNames, length.out=10), collapse=", "),"...")
			} else {
				colNstr <- paste(object@colNames, collapse=", ")
			}
			cat("  col names: ", colNstr, "\n", sep="")
		}
		# show(object[,])
		invisible(NULL)	
	}
)
setMethod("nrow", signature(x="BigFfMat"),
	function(x){
		x@rowN
	}
)
setMethod("ncol", signature(x="BigFfMat"),
	function(x){
		x@colN
	}
)
setMethod("dim", signature(x="BigFfMat"),
	function(x){
		c(x@rowN, x@colN)
	}
)
setMethod("rownames", signature(x="BigFfMat"),
	function(x){
		x@rowNames
	}
)
setReplaceMethod("rownames", signature(x="BigFfMat"),
	function(x, value){
		if (!is.null(value) && length(value)!=x@rowN){
				stop("Number of replacement items does not fit")
		}
		x@rowNames <- value
		invisible(x)
	}
)
setMethod("colnames", signature(x="BigFfMat"),
	function(x){
		x@colNames
	}
)
setReplaceMethod("colnames", signature(x="BigFfMat"),
	function(x, value){
		if (!is.null(value) && length(value)!=x@colN){
			stop("Number of replacement items does not fit")
		}
		x@colNames <- value
		invisible(x)
	}
)

#TODO: -indices
setMethod("[", "BigFfMat",
	function(x, i, j, drop=TRUE){
		if (missing(i)) i <- 1:x@rowN
		if (missing(j)) j <- 1:x@colN
		if (is.logical(i)) i <- which(i)
		if (is.logical(j)) j <- which(j)
		if (is.character(i)){
			if (all(i %in% x@rowNames)){
				i <- match(i, x@rowNames)
			} else {
				stop("Invalid row names specified")
			}
		}
		if (is.character(j)){
			if (all(j %in% x@colNames)){
				j <- match(j, x@colNames)
			} else {
				stop("Invalid col names specified")
			}
		}
		if (!is.vector(i)){
			stop("Invalid row selection (non-vector)")
		}
		if (!is.vector(j)){
			stop("Invalid column selection (non-vector)")
		}
		# print(paste("get: i:",paste(i,collapse=","),"j:",paste(j,collapse=",")))
		res <- do.call("cbind", lapply(x@cols[j],FUN=function(cc){
			if (x@activeClosing) open.ff(cc)
			res <- suppressMessages(cc[i])
			if (x@activeClosing) close.ff(cc)
			res
		}))
		colnames(res) <- x@colNames[j]
		rownames(res) <- x@rowNames[i]
		if (drop) res <- drop(res)
		return(res)
	}
)

setReplaceMethod("[", "BigFfMat",
	function(x, i, j, value){
		if (missing(i)) i <- 1:x@rowN
		if (missing(j)) j <- 1:x@colN
		if (is.logical(i)) i <- which(i)
		if (is.logical(j)) j <- which(j)
		if (is.character(i)){
			if (all(i %in% x@rowNames)){
				i <- match(i, x@rowNames)
			} else {
				stop("Invalid row names specified")
			}
		}
		if (is.character(j)){
			if (all(j %in% x@colNames)){
				j <- match(j, x@colNames)
			} else {
				stop("Invalid col names specified")
			}
		}
		if (!is.vector(i)){
			stop("Invalid row selection (non-vector)")
		}
		if (!is.vector(j)){
			stop("Invalid column selection (non-vector)")
		}
		# print(paste("set: i:",paste(i,collapse=","),"j:",paste(j,collapse=",")))

		if (is.vector(value)){
			if (length(value)==1){
				if (x@activeClosing) {
					for (jjj in j) open.ff(x@cols[[jjj]])
				}
				x@cols[j] <- lapply(1:length(j), FUN=function(jj){
					ccMod <- x@cols[j][[jj]]
					ccMod[i] <- value
					return(ccMod)
				})
				if (x@activeClosing) {
					for (jjj in j) close.ff(x@cols[[jjj]])
				}
			} else if (length(i)==length(value) && length(j)==1){
				if (x@activeClosing) open.ff(x@cols[[j]])
				x@cols[[j]][i] <- value 
				if (x@activeClosing) close.ff(x@cols[[j]])
			} else if (length(j)==length(value) && length(i)==1){
				if (x@activeClosing) {
					for (jjj in j) open.ff(x@cols[[jjj]])
				}
				x@cols[j] <- lapply(1:length(value), FUN=function(jj){
					ccMod <- x@cols[j][[jj]]
					ccMod[i] <- value[jj]
					return(ccMod)
				})
				if (x@activeClosing) {
					for (jjj in j) close.ff(x@cols[[jjj]])
				}
			} else {
				stop("Invalid specification for replacement object of type vector.")
			}
		} else {
			if (nrow(value)!=length(i) || ncol(value)!=length(j)){
				stop("Incompatible dimensions of replacement object and indices")
			}
			# x@cols[j] <- lapply(1:ncol(value), FUN=function(jj){
			# 	ccMod <- x@cols[j][[jj]]
			# 	ccMod[i] <- value[,jj]
			# 	return(ccMod)
			# })
			dummy <- lapply(1:ncol(value), FUN=function(jj){
				if (x@activeClosing) open.ff(x@cols[j][[jj]])
				x@cols[j][[jj]][i] <<- value[,jj]
				if (x@activeClosing) close.ff(x@cols[j][[jj]])
				return(NULL)
			})
		}
		invisible(x)
	}
)

if (!isGeneric("delete")) setGeneric("delete", function(x) standardGeneric("delete"))
setMethod("delete", signature(x="BigFfMat"),
	function(x){
		# message("removing BigFfMat object") #TODO: remove me
		for (j in 1:x@colN){
			ff::delete(x@cols[[j]])
		}
		invisible(NULL)
	}
)
# adapt vmode from ff package
setMethod("vmode", signature(x="BigFfMat"),
	function(x){
		if (ncol(x)<1) return(NA)
		return(vmode(x@cols[[1]]))
	}
)

#' save.bigFfMat
#'
#' save an \code{\linkS4class{BigFfMat}} object to disk
#'
#' @param bff \code{\linkS4class{BigFfMat}} object
#' @param file path on the disk to save to.
#' @param ...  arguments passed on to \code{ffsave}
#' @return Nothing of particular interest
#' @author Fabian Mueller
#' @aliases save.bigFfMat
#' @noRd
save.bigFfMat <- function(bff, file, ...){
	if (dir.exists(file) || file.exists(file)){
		stop("Destination directory or file already exists")
	}
	dir.create(file)
	objFn <- file.path(file,"bff.RData")
	save(bff, file=objFn)
	ee <- new.env()
	for (j in 1:bff@colN){
		cnCur <- paste0("col",j)
		ee[[cnCur]] <- bff@cols[[j]]
	}
	ffFn <- file.path(file,"bff_data")
	ffsave(list=ls(ee), envir=ee, file=ffFn, ...)
	invisible(NULL)
}

#' load.bigFfMat
#'
#' loads a saved \code{\linkS4class{BigFfMat}} object from disk
#'
#' @param path path of the saved object (a directory containing a corresponding \code{bff.RData} file and possibly \code{bff_data} files)
#' @param ...  arguments passed on to \code{ffload}
#' @return the loaded \code{\linkS4class{BigFfMat}} object
#' @author Fabian Mueller
#' @aliases load.bigFfMat
#' @noRd
load.bigFfMat <- function(path, ...){
	if (!file.exists(path)){
		stop("Loading unsuccessfull. Path does not exist")
	}
	load.env <- new.env(parent=emptyenv())
	load(file.path(path,"bff.RData"),envir=load.env)
	bff <- get("bff", load.env)
	eec <- new.env(parent=emptyenv())
	save.file <- file.path(path,"bff_data")
	ffload(save.file, envir=eec, ...)
	bff@cols <- lapply(1:bff@colN, FUN=function(j){
		cnCur <- paste0("col",j)
		get(cnCur,eec)
	})
	return(bff)
}
