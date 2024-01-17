########################################################################################################################
## readGEO.R
## created: 2012-03-27
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Loading array-based methylation datasets from GEO.
########################################################################################################################

## G L O B A L S #######################################################################################################

GEO.PLATFORMS <- c(
	"GPL8490" = "probes27",
	"GPL13534" = "probes450",
	"GPL16304" = "probes450",
	"GPL21145" = "probesEPIC",
	"GPL33022" = "probesEPICv2") ## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL33022

## F U N C T I O N S ###################################################################################################

#' Preprocesse sample characteristics
#'
#' Preprocesses sample characteristics lines from a series matrix definition file.
#'
#' @param txt               Lines in the file to be parsed.
#' @param regex.sample.info Regular expression capturing a definition of sample characteristics.
#' @return \code{character} vector storing sequence of preprocessed lines defining sample characteristics; an empty
#'         \code{vector} if none of the lines in \code{txt} matches the given regular expression.
#'
#' @author Yassen Assenov
#' @noRd
rnb.geo.parse.sample.info <- function(txt, regex.sample.info = "^!Sample_([^\t]+)\\t(.+)$") {
	s.info.lines <- grep(regex.sample.info, txt, value = TRUE)
	if (length(s.info.lines) == 0) {
		return(character())
	}
	result <- gsub(regex.sample.info, "\\2", s.info.lines)
	names(result) <- gsub(regex.sample.info, "\\1", s.info.lines)
	result
}

########################################################################################################################

#' Parses sample annotation
#'
#' Parses the sample annotation table from the corresponding lines in a series matrix definition file.
#'
#' @param txt Lines of the series matrix definition file that define sample characteristics. These should be
#'            preprocessed and not contain the initial sample specification (the token starting with \code{"!Sample_"}).
#' @return Sample annotation table as a \code{data.frame}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.geo.parse.sample.anno <- function(txt) {

	## Validate that the sample annotation forms a table
	regex.sample.value <- "^\"([^\"]*)\"$"
	s.infos <- strsplit(txt, "\t", fixed = TRUE)
	s.valid <- mapply(grepl, x = s.infos, MoreArgs = list(pattern = regex.sample.value), SIMPLIFY = FALSE)
	if (!all(unlist(s.valid))) {
		stop("Cannot parse sample information; quotation marks and/or TAB characters are probably escaped")
	}
	s.infos <- lapply(s.infos, function(x) { gsub(regex.sample.value, "\\1", x) })
	sample.count <- range(sapply(s.infos, length))
	if (sample.count[1] != sample.count[2]) {
		stop("Incosistent number of samples in the annotation")
	}
	sample.count <- sample.count[1]

	## Parse characteristics_ch1
	for (i in which(names(s.infos) == "characteristics_ch1")) {
		ch1.values <- strsplit(s.infos[[i]], ": ", fixed = TRUE)
		if (all(sapply(ch1.values, length) %in% c(0, 2))) {
			j <- which(sapply(ch1.values, length) == 2)
			if (length(j) == 0) { next }
			column.title <- unique(sapply(ch1.values[j], '[', 1L))
			if (length(column.title) == 1) {
				s.infos[[i]] <- sapply(ch1.values, '[', 2L)
				names(s.infos)[i] <- column.title
			}
		}
	}
	suppressWarnings(rm(i, ch1.values, j, column.title))

	## Detect column types
	s.infos <- as.data.frame(s.infos, check.names = FALSE, stringsAsFactors = FALSE)
	buffer <- character()
	b.connection <- textConnection('buffer', open = "w", local = TRUE)
	write.table(s.infos, file = b.connection, quote = FALSE, sep = "\t", row.names = FALSE)
	close(b.connection)
	b.connection <- textConnection(buffer)
	s.infos <- read.delim(b.connection, quote = "", check.names = FALSE)
	close(b.connection)
	for (i in 1:ncol(s.infos)) {
		if (is.factor(s.infos[, i]) && nlevels(s.infos[, i]) == sample.count) {
			s.infos[, i] <- as.character(s.infos[, i])
		}
	}
	## Auto-rename columns in case there are duplicates
	if (ncol(s.infos) > 1) {
		s.infos <- s.infos[, rep(TRUE, ncol(s.infos))]
	}
	s.infos
}

########################################################################################################################

#' Parse sample identifiers
#'
#' Parses column names from the first line of series matrix table definition.
#' @param txt Line in the series matrix file that defines table columns; this is the line following the matrix
#'            definition \code{"series_matrix_table_begin"}.
#' @return All column names in the form of a \code{character} vector.
#'
#' @author Yassen Assenov
#' @noRd
rnb.geo.parse.ids <- function(txt) {
	txt <- strsplit(txt, "\t", fixed = TRUE)[[1]]
	if (length(txt) < 2 || txt[1] != '"ID_REF"') {
		stop("Unexpected header of series matrix table")
	}
	regex.id <- "^\\\"(GSM\\d+)\\\"$"
	if (!all(grepl(regex.id, txt[-1]))) {
		stop("Unexpected header of series matrix table")
	}
	gsub(regex.id, "\\1", txt[-1])
}

########################################################################################################################

#' Initialize methylation matrix
#'
#' Initializes the beta value matrix from a series matrix file.
#' @param txt        Line in the series matrix file that defines table columns; this is the line following the matrix
#'                   definition \code{"series_matrix_table_begin"}.
#' @param N.expected Number of columns expected, based on the sample definition lines above.
#' @param assay.type Expected assay. This must be one of \code{"probes27"}, \code{"probes450"} or \code{"probesEPIC"}.
#' @return Newly initialized \code{matrix} of type \code{double}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.geo.init.matrix <- function(txt, N.expected, assay.type) {
	## TODO: Does it need implementation for EPICv2?
	sample.ids <- rnb.geo.parse.ids(txt)
	if (length(sample.ids) != N.expected) {
		stop("Inconsistent sample characteristics and series matrix table")
	}
	genome.assembly <- rnb.getOption("assembly")
	probe.ids <- rnb.get.annotation(assay.type, assembly = genome.assembly) ## TODO: Improve genome build selection
	probe.ids <- unlist(lapply(probe.ids, names), use.names = FALSE)

	## TODO: Consider rnb.getOption("disk.dump.big.matrices") and initialize different things
	matrix(NA_real_, length(probe.ids), length(sample.ids), dimnames = list(probe.ids, sample.ids))
}

########################################################################################################################

#' Parse a series matrix file
#'
#' Parses the sample annotation data and methylation values from a series matrix file.
#'
#' @param fname   Name of the file that contains the series matrix definition. This can also be gzipped.
#' @param verbose Flag indicating if messages are to be sent to the logger.
#' @return \code{list} with three elements: sample annotation table (\code{data.frame}), methylation beta value matrix,
#'         and a platform, specified as one of \code{"probes27"}, \code{"probes450"} or \code{"probesEPIC"}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.geo.parse.series.matrix <- function(fname, verbose) {

	if (grepl("\\.gz$", fname)) {
		con <- gzfile(fname, "rt")
	} else {
		con <- file(fname, "rt")
	}
	if (verbose) {
		logger.start("Parsing Series Matrix File")
	}

	N.buffer <- 100000L # size of buffer, in number of lines
	regex.platform <- "^\\!Series_platform_id\t\\\"(GPL\\d+)\\\"$"
	data.platform <- NULL
	s.infos <- character()
	state.data.matrix <- 0L # 0: not started, 1: parsing header, 2: parsing data, 3: ready
	data.matrix <- NULL
	regex.data <- "^\\\"([^\"]+)\\\"\t(.*)$"
	skipped.records <- 0L

	repeat {
		txt <- NULL; invisible(gc())
		txt <- scan(con, "", n = N.buffer, sep = "\n", quiet = TRUE)
		if (length(txt) == 0) {
			break
		}

		## Parse platform
		i <- grep(regex.platform, txt)
		if (length(i) == 1) {
			if (!is.null(data.platform)) {
				stop("Multiple platforms defined")
			}
			data.platform <- gsub(regex.platform, "\\1", txt[i])
			if (!(data.platform %in% names(GEO.PLATFORMS))) {
				stop(paste("Platform", data.platform, "is not supported"))
			}
			if (verbose) {
				i <- paste0("Identified platform ", data.platform, " (", GEO.PLATFORMS[data.platform], ")")
				logger.status(i)
			}
		} else if (length(i) > 1) {
			stop("Multiple platforms defined; RnBeads supports a single platform only")
		}

		i <- which(txt == "!series_matrix_table_begin")
		if (length(i) == 0) {
			if (state.data.matrix == 0L) {
				s.infos <- c(s.infos, rnb.geo.parse.sample.info(txt))
				next
			}
			if (state.data.matrix == 1L) {
				data.matrix <- rnb.geo.init.matrix(txt[1], nrow(s.infos), GEO.PLATFORMS[data.platform])
				i.first <- 2L
				state.data.matrix <- 2L
			} else if (state.data.matrix == 2L) {
				i.first <- 1L
			}
		} else if (length(i) == 1 && state.data.matrix == 0L) {
			## Parse sample annotation table from the identified lines
			if (i > 1) {
				s.infos <- c(s.infos, rnb.geo.parse.sample.info(head(txt, i - 1)))
			}
			s.infos <- rnb.geo.parse.sample.anno(s.infos)
			## Initialize data matrix
			if (i + 1L <= length(txt)) {
				data.matrix <- rnb.geo.init.matrix(txt[i + 1], nrow(s.infos), GEO.PLATFORMS[data.platform])
				state.data.matrix <- 2L
				if (i + 1L == length(txt)) {
					next
				}
				i.first <- i + 2L
			} else {
				state.data.matrix <- 1L
				next
			}
		} else {
			stop("Multiple series matrix defined; RnBeads supports a single series matrix only")
		}

		i <- which(txt == "!series_matrix_table_end")
		if (length(i) == 0) {
			if (state.data.matrix != 2L) {
				next
			}
			i.last <- length(txt)
		} else if (length(i) == 1 && state.data.matrix == 2L) {
			i.last <- i - 1L
			if (i.last < i.first) {
				stop("Empty series matrix definition")
			}
		} else {
			stop("Multiple or misplaced end definition(s) of a series matrix")
		}

		## Parse lines to fill in the beta value matrix
		if (state.data.matrix == 2L) {
			txt <- txt[i.first:i.last]
			i.probes <- grep(regex.data, txt, invert = TRUE)
			if (length(i.probes) != 0) {
				stop("Invalid format of data table")
			}
			i.probes <- match(gsub(regex.data, "\\1", txt), rownames(data.matrix))
			i.found <- which(i.probes != 0L)
			skipped.records <- skipped.records + length(i.probes) - length(i.found)
			if (length(i.found) != 0) {
				txt <- gsub(regex.data, "\\2", txt[i.found])
				i.probes <- i.probes[i.found]
				txt <- gsub("\t$", "\t\t", txt)
				mm <- suppressWarnings(sapply(strsplit(txt, "\t", fixed = TRUE), as.double))
				if (!(is.matrix(mm) && nrow(mm) == ncol(data.matrix))) {
					stop("Invalid format of data table")
				}
				data.matrix[i.probes, ] <- t(mm)
			}
			suppressWarnings(rm(i.probes, i.found, mm))

			if (length(i) == 1) {
				state.data.matrix <- 3L
			}
		}
	}
	close(con)

	if (!(is.data.frame(s.infos) && state.data.matrix == 3L)) {
		stop("Missing or incomplete data table")
	}
	if (skipped.records != 0 && verbose) {
		i <- ifelse(skipped.records == 1, "record", "records")
		logger.warning(paste("Skipped", skipped.records, i, "referring to unsupported probes"))
	}
	data.matrix <- data.matrix[!apply(is.na(data.matrix), 1, all), , drop = FALSE]
	colnames(data.matrix) <- NULL
	if (verbose) {
		logger.status(paste("Loaded non-missing data for", nrow(data.matrix), "probes"))
		logger.completed()
	}
	list("annotation" = s.infos, "betas" = data.matrix, "platform" = GEO.PLATFORMS[data.platform])
}

########################################################################################################################

#' Import methylation data from GEO
#'
#' Imports Infinium 450K, MethylationEPIC or MethylationEPIC2 data series from the Gene Expression Omnibus. This function uses the
#' series matrix file.
#'
#' @param accession  Character string, starting with \code{"GSE"}, representing the GEO series for download and parsing.
#'                   Alternatively, this parameter can specify the file name of a previously downloaded GEO series
#'                   matrix file or its gzipped representation (in which case the filename must end in \code{".gz"}).
#'                   Other file formats, such as SOFT files, are not supported.
#' @param verbose    Flag indicating if messages should be created informing about the progress. If the logger is
#'                   initialized prior to calling this function, the informative messages are sent to the logger.
#'                   Warnings and errors are not affected by this parameters, the function always outputs them.
#' @param destdir    The destination directory for any downloads. Defaults to the (architecture-dependent) temporary
#'                   directory. Keep in mind that GEO series can be demanding in terms of storage space.
#' @return \code{\linkS4class{RnBeadSet}} object with phenotypic and beta value information.
#'
#' @author Yassen Assenov, modified by Baris Kalem
#' @export
rnb.read.geo <- function(accession = NULL, verbose = logger.isinitialized(), destdir = tempdir()) {
	## TODO: Make it compatible with EPIC v2
	if (verbose) {
		rnb.logger.start("Loading GEO Data Series")
	}
	if (!(is.character(accession) && length(accession) == 1 && isTRUE(accession != ""))) {
		stop("Invalid value for accession")
	}
	if (isTRUE(file.info(accession)[1, "isdir"] == FALSE)) {
		## Specified existing file
		fname <- accession
	} else if (!grepl("^GSE\\d\\d+$", accession)) {
		stop("Invalid value for accession; expected GEO data series identifier or an existing file")
	} else {
		## Specified GSE data series, download from GEO
		fname <- tempfile("gse", destdir, ".txt.gz")
		gse.url <- gsub("^(GSE\\d\\d).*$", "\\1", accession)
		gse.url <- sprintf("%snnn/%s/matrix/%s", gse.url, accession, accession)
		gse.url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", gse.url, "_series_matrix.txt.gz")
		if (verbose) {
			rnb.info(paste("Using URL", gse.url))
			rnb.info(paste("Using download method", getOption("download.file.method")))
		}
		result <- tryCatch( # quiet = TRUE would be cleaner, but it seems to fail on MacOS
			download.file(gse.url, fname, quiet = FALSE, mode = "wb"),
			error = function(err) { return(-99L) })
		if (result != 0) {
			invisible(file.remove(fname))
			stop(paste("Could not download GEO Data Series, error code", result))
		}
		if (verbose) {
			rnb.status(paste("Downloaded GSE matrix to", fname))
		}
		rm(gse.url, result)
	}

	## Parse the series matrix file
	result <- rnb.geo.parse.series.matrix(fname, verbose)
	x <- c("probes27" = "27k", "probes450" = "450k", "probesEPIC" = "EPIC", "probesEPICv2" = "EPICv2")[result[[3]]]
	RnBeadSet(pheno = result[[1]], betas = result[[2]], platform = x)
}
