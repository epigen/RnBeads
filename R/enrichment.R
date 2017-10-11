########################################################################################################################
## differentialMethylation.R
## created: 2017-07-04
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for enrichment analysis.
########################################################################################################################

################################################################################
# GO enrichment analyses
################################################################################
#' performGOenrichment.diffMeth.entrez
#'
#' performs Gene Ontology (GO) enrichment analysis for a list of Entrez identifiers
#' @param gids gene ids to test (entrez IDs)
#' @param uids ids to test against (universe)
#' @param ontology which ontology should be used (see \code{GOHyperGParams} from the \code{GOstats} package for details)
#' @param assembly Genome to be used. One of the following: hg19, mm9, mm10 or rn5
#' @param ... arguments passed on to the parameters of \code{GOHyperGParams} from the \code{GOstats} package
#' @return a \code{GOHyperGresult} object (see the \code{GOstats} package for further details)
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' dmt <- get.table(dm,get.comparisons(dm)[1],"promoters")
#' annot <- annotation(rnb.set.example,"promoters")
#' all.promoters <- annot$entrezID
#' #get the hypermethylated promoters
#' hyper.promoters <- annot$entrezID[dmt[,"mean.mean.diff"]>0]
#' result <- performGOenrichment.diffMeth.entrez(hyper.promoters,all.promoters,"BP",assembly="hg19")
#' }
performGOenrichment.diffMeth.entrez <- function(gids,uids,ontology,assembly="hg19",...){
	rnb.require("GOstats")
#	#testing
#	dmt <- diffmeth$region[[1]][["promoters"]]
#	annot <- rnb.get.annotations(type = "promoters")
#	uids <- unique(unlist(sapply(annot[rownames(dmt),]$entrezid,FUN=function(idstr){
#		unlist(strsplit(idstr, ";", fixed = TRUE))
#	})))
#	dmrs <- dmt[dmt[,"combinedRank"] < 1000,]
#	gids <- unique(unlist(sapply(annot[rownames(dmrs),]$entrezid,FUN=function(idstr){
#		unlist(strsplit(idstr, ";", fixed = TRUE))
#	})))
#	#/testing
	if (is.element(assembly,c("hg19","hg38"))){
		ass <- "org.Hs.eg.db"
	} else if (is.element(assembly,c("mm9","mm10"))){
		ass <- "org.Mm.eg.db"
	} else if (is.element(assembly,c("rn5"))){
		ass <- "org.Rn.eg.db"
	} else if (is.element(assembly,c("zv9"))){
		ass <- "org.Dr.eg.db"
	} else {
		stop("Unsupported assembly")
	}
	params <- new("GOHyperGParams",annotation=ass,geneIds = gids, universeGeneIds = uids, 
			ontology = ontology,conditional = TRUE, testDirection = "over",...)
	nonMatchIdErrorMsg <- c(
		"The genes you are testing do not have any corresponding GO terms for the ontology you are searching.",
		"genes being tested do not have corresponding GO terms"
	)
	hgOver <- tryCatch({
			GOstats::hyperGTest(params)
		}, error = function(ee) {
			if (is.element(ee$message, nonMatchIdErrorMsg)){
				logger.info("Could not conduct enrichment analysis as associated genes are not in GO database.")
				NULL
			} else {
				logger.warning(c("Could not conduct enrichment analysis:", ee$message))
				NULL
			}
		}
	)
	return(hgOver)
}

#' performGoEnrichment.diffMeth
#'
#' performs Geno Ontology (GO) enrichment analysis for a given differential methylation table.
#' @author Fabian Mueller
#' @param rnbSet RnBSet object for which dirrential methylation was computed
#' @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param ontologies GO ontologies to use for enrichment analysis
#' @param rank.cuts.region Cutoffs for combined ranking that are used to determine differentially methylated regions
#' @param add.auto.rank.cut flag indicating whether an automatically computed cut-off should also be considered.
#' @param rerank For deterimining differential methylation: should the ranks be ranked again or should the absolute ranks be used.
#' @param verbose Enable for detailed status report
#' @param ... arguments passed on to the parameters of \code{GOHyperGParams} from the \code{GOstats} package
#' @return a DiffMeth.go.enrich object (S3) containing the following attributes
#' \item{region}{Enrichment information for differential methylation on the region level. See \code{GOHyperGresult} from the \code{GOstats} package for furthert details}
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' res <- performGoEnrichment.diffMeth(rnb.set.example,dm)
#' }
performGoEnrichment.diffMeth <- function(rnbSet,diffmeth,ontologies=c("BP","MF"),rank.cuts.region=c(100,500,1000),add.auto.rank.cut=TRUE,rerank=TRUE,verbose=TRUE,...){
	gene.id.col <- "entrezID"
	logger.start("Differential Methylation GO Enrichment Analysis")
	if (!is.element(assembly(rnbSet),c("hg19","hg38","mm9","mm10","rn5","zv9"))){
		logger.warning(c("Enrichment analysis currently not supported for genome assembly:",assembly(rnbSet),"--> skipping enrichment analysis"))
		logger.completed()
		return(NULL)
	}
	comps <- get.comparisons(diffmeth)
	region.types <- get.region.types(diffmeth)
	skipSites <- !includes.sites(diffmeth)
	diff.col.reg <- "mean.mean.diff"
	if (skipSites) diff.col.reg <- "mean.diff"
	dm.go.enrich <- list(probe=list(),region=list())
	for (cc in comps){
		logger.start(c("Comparison: ",cc))
		dm.go.enrich$probe[[cc]] <- list()
		dm.go.enrich$region[[cc]] <- list()
		for (oo in ontologies){
			logger.start(c("Ontology: ",oo))
			logger.start("Region Level")
			dm.en.region.list.list <- list()
			for (rt in region.types){
				logger.start(c("Region Type:",rt))
				#rt <- "promoters"
				annot <- annotation(rnbSet,type = rt)
				if (!(gene.id.col %in% colnames(annot))){
					logger.info(c("Not annotated with",gene.id.col,"--> Skipped"))
					logger.completed()
					next
				}
#				dmt <- diffmeth$region[[cc]][[rt]]
				dmt <- get.table(diffmeth,cc,rt,undump=TRUE,return.data.frame=TRUE)
				rrs <- dmt[,"combinedRank"]
				rrs.hyper <- rrs
				rrs.hypo <- rrs
				rrs.hyper[dmt[,diff.col.reg] <= 0] <- NA
				rrs.hypo[dmt[,diff.col.reg] >= 0] <- NA
				#automatically select rank cutoff
				rc.auto <- 0L
				if (add.auto.rank.cut){
					rc.auto <- as.integer(auto.select.rank.cut(dmt$comb.p.adj.fdr,dmt$combinedRank))
				}

				if (rerank){
					rrs <- rank(rrs,na.last="keep",ties.method="min")
					rrs.hyper <- rank(rrs.hyper,na.last="keep",ties.method="min")
					rrs.hypo <- rank(rrs.hypo,na.last="keep",ties.method="min")
					if (add.auto.rank.cut && rc.auto > 0L){
						rc.auto <- na.omit(rrs[dmt[,"combinedRank"]==rc.auto])[1] #arbitrarily select the first rerank if the rank cut threshold is occupied by multiple ranks
					}
				}
				
				gene.ids.chr <- as.character(annot[,gene.id.col])
				uids <- unique(unlist(sapply(gene.ids.chr,FUN=function(idstr){
					unlist(strsplit(idstr, ";", fixed = TRUE))
				})))
				uids <- na.omit(uids)
				dm.en.region.list <- list()
				rank.cuts <- rank.cuts.region
				if (add.auto.rank.cut){
					rank.cuts <- c(rank.cuts,"autoSelect")
				}
				for (i in 1:length(rank.cuts)){
					rc <- rank.cuts[i]
					if (rc == "autoSelect") {
						rc <- rc.auto
						if (verbose) logger.info(c("Rank cutoff:",rc,"(auto-select)"))
					} else {
						rc <- as.integer(rc)
						if (verbose) logger.info(c("Rank cutoff:",rc))
					}
					is.dmr.hyper <- rrs.hyper <= rc
					is.dmr.hypo <- rrs.hypo <= rc
					gids.hyper <- unique(unlist(sapply(gene.ids.chr[is.dmr.hyper],FUN=function(idstr){
						unlist(strsplit(idstr, ";", fixed = TRUE))
					})))
					gids.hyper <- na.omit(gids.hyper)
					if (length(gids.hyper)>0 && length(uids)>0){
						dm.en.region.hyper <- performGOenrichment.diffMeth.entrez(gids.hyper,uids,ontology=oo,assembly=assembly(rnbSet),...)
					} else {
						dm.en.region.hyper <- NULL
					}
					
					gids.hypo <- unique(unlist(sapply(gene.ids.chr[is.dmr.hypo],FUN=function(idstr){
						unlist(strsplit(idstr, ";", fixed = TRUE))
					})))
					gids.hypo <- na.omit(gids.hypo)
					if (length(gids.hypo)>0 && length(uids)>0){
						dm.en.region.hypo <- performGOenrichment.diffMeth.entrez(gids.hypo,uids,ontology=oo,assembly=assembly(rnbSet),...)
					} else {
						dm.en.region.hypo <- NULL
					}
					dm.en.region.list <- c(dm.en.region.list,list(list(hyper=dm.en.region.hyper,hypo=dm.en.region.hypo)))
				}
				names(dm.en.region.list) <- paste("rankCut_",rank.cuts,sep="")
				dm.en.region.list.list <- c(dm.en.region.list.list,list(dm.en.region.list))
				names(dm.en.region.list.list)[length(dm.en.region.list.list)] <- rt
				logger.completed()
			}
			dm.go.enrich$region[[cc]][[oo]] <- dm.en.region.list.list
			logger.completed()
			logger.completed()
		}
		logger.completed()
	}
	class(dm.go.enrich) <- "DiffMeth.go.enrich"
	logger.completed()
	return(dm.go.enrich)
}

################################################################################
# LOLA enrichment analyses
################################################################################
#' downloadLolaDbs
#' 
#' Downloading prepared LOLA DBs from server
#' @param dest    destination directory
#' @param dbs     vector of names of LOLA DBs to be downloaded. Currently 'LOLACore' and 'LOLAExt'
#'                are supported
#' @details
#' Requires a stable internet connection. Could take a while depending on the size of the database and the internet connection
#' @return a list containing vectors of directory names for each available genome assembly
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' }
downloadLolaDbs <- function(dest, dbs=c("LOLACore")){
	res <- list()
	for (dbName in dbs){
		dbDownloadLink <- NULL
		if (dbName == "LOLACore"){
			dbDownloadLink <- "http://regiondb.s3.amazonaws.com/LOLACoreCaches_170206.tgz"
		} else if (dbName == "LOLAExt"){
			dbDownloadLink <- "http://regiondb.s3.amazonaws.com/LOLAExtCaches_170206.tgz"
		}
		if (is.null(dbDownloadLink)){
			logger.error(c("Invalid DB name:", dbName))
		}

		dbDir <- file.path(dest, dbName)
		if (dir.exists(dbDir)){
			logger.warning(c("Directory", dbDir, "already existing --> skipping download"))
		} else {
			tmpFn <- tempfile(fileext=".tgz")
			logger.start(c("Downloading from", dbDownloadLink))
			download.file(dbDownloadLink, destfile=tmpFn, mode = "wb")
			exDir <- tempfile()
			untar(tmpFn, exdir=exDir)
			baseDir <- grep(paste0(dbName, "$"), list.dirs(exDir, recursive=TRUE), value=TRUE)
			if (length(baseDir) < 1) logger.error(c("Downloaded archive does not contain a subdirectory", dbName))
			dir.create(dbDir)
			file.copy(baseDir, dest, recursive=TRUE)
			unlink(c(tmpFn, exDir), recursive=TRUE) # remove temp files
			logger.completed()
		}
		assemblies <- list.dirs(dbDir, full.names=FALSE, recursive=FALSE)
		for (aa in assemblies){
			res[[aa]] <- c(res[[aa]], file.path(dbDir, aa))
		}
	}
	return(res)
}

#' prepLolaDbPaths
#' 
#' prepare LOLA Database paths from options
#' @param assembly name of the assembly/subdirectory of the DB path
#' @param dbs  LOLA database paths, potentially containing placeholders
#' @param downloadDir directory where default LOLA DBs corresponding to placeholders will be downloaded to
#'             (provided that placeholders are contained in \code{dbs})
#' 
#' @return vector of LOLA DB paths
#' @author Fabian Mueller
#' @noRd
prepLolaDbPaths <- function(assembly, dbs=rnb.getOption("differential.enrichment.lola.dbs"), downloadDir=tempdir()){
	#download DBs for placeholders
	phPattern <- "^\\$\\{(.+)\\}$"
	isPh <- grepl(phPattern, dbs)
	res <- file.path(dbs[!isPh], assembly)
	if (any(isPh)){
		logger.start("Downloading LOLA databases")
			if (!dir.exists(downloadDir)) dir.create(downloadDir)
			downloadDbs <- gsub(phPattern, "\\1", dbs[isPh])
			downloadedPaths <- downloadLolaDbs(downloadDir, dbs=downloadDbs)
		logger.completed()
		res <- c(res, downloadedPaths[[assembly]])
	}
	nonExisting <- !dir.exists(res)
	if (any(nonExisting)){
		logger.error(c("The following LOLA DB directories are missing:", paste(res[nonExisting], collapse="; ")))
	}
	return(res)
}

#' loadLolaDbs
#' 
#' Load LOLA databases from disk and merge them
#' @param lolaDbPaths  vector of names of LOLA DB paths to be loaded
#' 
#' @return LOLA DB list as returned by \code{LOLA::loadRegionDB}
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' }
loadLolaDbs <- function(lolaDbPaths){
	require(data.table) #explicitely load data.table to adress LOLA namespace issues
	require(LOLA)
	require(simpleCache) # TODO: include requirement in dependencies
	logger.start("Loading LOLA DBs")
		lolaDb <- loadRegionDB(lolaDbPaths[1])
		if (length(lolaDbPaths)>1){
			for (i in 2:length(lolaDbPaths)){
				lolaDb <- mergeRegionDBs(lolaDb, loadRegionDB(lolaDbPaths[i]))
			}
		}
	logger.completed()
	return(lolaDb)
}

#' performLolaEnrichment.diffMeth
#'
#' performs LOLA enrichment analysis for a given differential methylation table.
#' @author Fabian Mueller
#' @param rnbSet RnBSet object for which dirrential methylation was computed
#' @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param lolaDbPaths LOLA database paths
#' @param rank.cuts.region Cutoffs for combined ranking that are used to determine differentially methylated regions
#' @param add.auto.rank.cut flag indicating whether an automatically computed cut-off should also be considered.
#' @param rerank For deterimining differential methylation: should the ranks be ranked again or should the absolute ranks be used.
#' @param verbose Enable for detailed status report
#' @return a DiffMeth.lola.enrich object (S3) containing the following attributes
#' \item{region}{Enrichment information for differential methylation on the region level. A \code{data.table} object
#' as returned by the \code{runLOLA} function from the \code{LOLA} package for furthert details. Each element will contain different
#' user sets for different rank cutoffs and hyper/hypomethylation events(\code{userSet} column)}
#' \item{lolaDb}{The loaded \code{lolaDb} object containing the merged databases as returned by \code{\link{loadLolaDbs}}}
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
#' }
performLolaEnrichment.diffMeth <- function(rnbSet, diffmeth, lolaDbPaths, rank.cuts.region=c(100,500,1000), add.auto.rank.cut=TRUE, rerank=TRUE, verbose=TRUE){
	nCores <- ifelse(parallel.isEnabled(), parallel.getNumWorkers(), 1)
	lolaDb <- loadLolaDbs(lolaDbPaths)
	comps <- get.comparisons(diffmeth)
	region.types <- get.region.types(diffmeth)
	diff.col.reg <- "mean.mean.diff"
	dm.lola.enrich <- list(probe=list(),region=list())
	for (cc in comps){
		dm.lola.enrich$probe[[cc]] <- list()
		dm.lola.enrich$region[[cc]] <- list()
	}
	logger.start("Differential Methylation LOLA Enrichment Analysis")
	for (rt in region.types){
		logger.start(c("Region Type:",rt))
			annot <- annotation(rnbSet, type=rt)
			# convert annotaiton to GRanges object
			univGr <- data.frame2GRanges(annot, chrom.column="Chromosome", start.column="Start", end.column="End", assembly=assembly(rnbSet), sort.result=FALSE)

			lolaUserSets <- list()

			for (i in 1:length(comps)){
				cc <- comps[i]
				compKey <- paste0("comp", i)
				logger.start(c("Comparison: ",cc))
					dmt <- get.table(diffmeth,cc,rt,undump=TRUE,return.data.frame=TRUE)

					rrs <- dmt[,"combinedRank"]
					rrs.hyper <- rrs
					rrs.hypo <- rrs
					rrs.hyper[dmt[,diff.col.reg] <= 0] <- NA
					rrs.hypo[dmt[,diff.col.reg] >= 0] <- NA
					#automatically select rank cutoff
					rc.auto <- 0L
					if (add.auto.rank.cut){
						rc.auto <- as.integer(auto.select.rank.cut(dmt$comb.p.adj.fdr,dmt$combinedRank))
					}

					if (rerank){
						rrs <- rank(rrs,na.last="keep",ties.method="min")
						rrs.hyper <- rank(rrs.hyper,na.last="keep",ties.method="min")
						rrs.hypo <- rank(rrs.hypo,na.last="keep",ties.method="min")
						if (add.auto.rank.cut && rc.auto > 0L){
							rc.auto <- na.omit(rrs[dmt[,"combinedRank"]==rc.auto])[1] #arbitrarily select the first rerank if the rank cut threshold is occupied by multiple ranks
						}
					}

					rank.cuts <- rank.cuts.region
					if (add.auto.rank.cut){
						rank.cuts <- c(rank.cuts,"autoSelect")
					}
					for (j in 1:length(rank.cuts)){
						rcn <- paste("rankCut_", rank.cuts[j] ,sep="")

						rc <- rank.cuts[j]
						if (rc == "autoSelect") {
							rc <- rc.auto
							if (verbose) logger.info(c("Rank cutoff:",rc,"(auto-select)"))
						} else {
							rc <- as.integer(rc)
							if (verbose) logger.info(c("Rank cutoff:",rc))
						}
						is.dmr.hyper <- rrs.hyper <= rc
						is.dmr.hyper[is.na(is.dmr.hyper)] <- FALSE
						is.dmr.hypo <- rrs.hypo <= rc
						is.dmr.hypo[is.na(is.dmr.hypo)] <- FALSE

						lolaUserSets <- c(lolaUserSets, list(univGr[is.dmr.hyper]), list(univGr[is.dmr.hypo]))
						N <- length(lolaUserSets)
						names(lolaUserSets)[c(N-1,N)] <- paste(compKey, rcn, c("hyper", "hypo"), sep="_")
					}
				logger.completed()
			}

			lolaUserSets <- GRangesList(lolaUserSets)
			logger.start("Running LOLA")
				lolaRes <- runLOLA(lolaUserSets, univGr, lolaDb, cores=nCores)
			logger.completed()

			for (i in 1:length(comps)){
				cc <- comps[i]
				curPat <- paste0("^comp",i,"_")
				curSubset <- grepl(curPat, lolaRes[["userSet"]])
				lolaRes.cur <- lolaRes[curSubset,]
				lolaRes.cur[["userSet"]] <- gsub(curPat, "", lolaRes.cur[["userSet"]])
				dm.lola.enrich$region[[cc]][[rt]] <- lolaRes.cur
			}

		logger.completed()
	}
	dm.lola.enrich[["lolaDb"]] <- lolaDb
	class(dm.lola.enrich) <- "DiffMeth.lola.enrich"
	logger.completed()
	return(dm.lola.enrich)
}

#' performGOEnrichment.diffVar
#'
#' performs Geno Ontology (GO) enrichment analysis for a given differential variability table.
#' @author Fabian Mueller and Michael Scherer
#' @param rnbSet RnBSet object for which dirrential variability was computed
#' @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param enrich.diffMeth Result of \code{performGOEnrichment.diffMeth}. NULL, if enrichment should only be performed for differential variability.
#' @param ontologies GO ontologies to use for enrichment analysis
#' @param rank.cuts.region Cutoffs for combined ranking that are used to determine differentially variable regions
#' @param add.auto.rank.cut flag indicating whether an automatically computed cut-off should also be considered.
#' @param rerank For deterimining differential variability: should the ranks be ranked again or should the absolute ranks be used.
#' @param verbose Enable for detailed status report
#' @param ... arguments passed on to the parameters of \code{GOHyperGParams} from the \code{GOstats} package
#' @return a DiffMeth.enrich object (S3) containing the following attributes
#' \item{region}{Enrichment information for differential variability on the region level. See \code{GOHyperGresult} from the \code{GOstats} package for furthert details}
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.diffVar(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' res <- performEnrichment.diffVar(rnb.set.example,dm)
#' }
performGOEnrichment.diffVar <- function(rnbSet,diffmeth,enrich.diffMeth=NULL,ontologies=c("BP","MF"),rank.cuts.region=c(100,500,1000),add.auto.rank.cut=TRUE,rerank=TRUE,verbose=TRUE,...){
  gene.id.col <- "entrezID"
  logger.start("Differential Variability Enrichment Analysis")
  if (!is.element(assembly(rnbSet),c("hg19","hg38","mm9","mm10","rn5","zv9"))){
    logger.warning(c("Enrichment analysis currently not supported for genome assembly:",assembly(rnbSet),"--> skipping enrichment analysis"))
    logger.completed()
    return(NULL)
  }
  comps <- get.comparisons(diffmeth)
  region.types <- get.region.types(diffmeth)
  skipSites <- !includes.sites(diffmeth)
  diff.col.reg <- "mean.var.diff"
  if (skipSites) diff.col.reg <- "var.diff"
  if(!is.null(enrich.diffMeth)){
    enrich.diffMeth$probe_var <- list()
    enrich.diffMeth$region_var <- list()
    dv.enrich <- enrich.diffMeth
    rm(enrich.diffMeth)
  }else{
    dv.enrich <- list(probe_var=list(),region_var=list())
  }
  for (cc in comps){
    logger.start(c("Comparison: ",cc))
    dv.enrich$probe_var[[cc]] <- list()
    dv.enrich$region_var[[cc]] <- list()
    for (oo in ontologies){
      logger.start(c("Ontology: ",oo))
      logger.start("Region Level")
      dm.en.region.list.list <- list()
      for (rt in region.types){
        logger.start(c("Region Type:",rt))
        #rt <- "promoters"
        annot <- annotation(rnbSet,type = rt)
        if (!(gene.id.col %in% colnames(annot))){
          logger.info(c("Not annotated with",gene.id.col,"--> Skipped"))
          logger.completed()
          next
        }
        dmt <- get.table(diffmeth,cc,rt,undump=TRUE,return.data.frame=TRUE)
        rrs <- dmt[,"combinedRank.var"]
        rrs.hyper <- rrs
        rrs.hypo <- rrs
        rrs.hyper[dmt[,diff.col.reg] <= 0] <- NA
        rrs.hypo[dmt[,diff.col.reg] >= 0] <- NA
        #automatically select rank cutoff
        rc.auto <- 0L
        if (add.auto.rank.cut){
          rc.auto <- as.integer(auto.select.rank.cut(dmt$comb.p.adj.var.fdr,dmt$combinedRank.var))
        }
        
        if (rerank){
          rrs <- rank(rrs,na.last="keep",ties.method="min")
          rrs.hyper <- rank(rrs.hyper,na.last="keep",ties.method="min")
          rrs.hypo <- rank(rrs.hypo,na.last="keep",ties.method="min")
          if (add.auto.rank.cut && rc.auto > 0L){
            rc.auto <- na.omit(rrs[dmt[,"combinedRank.var"]==rc.auto])[1] #arbitrarily select the first rerank if the rank cut threshold is occupied by multiple ranks
          }
        }
        
        gene.ids.chr <- as.character(annot[,gene.id.col])
        uids <- unique(unlist(sapply(gene.ids.chr,FUN=function(idstr){
          unlist(strsplit(idstr, ";", fixed = TRUE))
        })))
        uids <- na.omit(uids)
        dm.en.region.list <- list()
        rank.cuts <- rank.cuts.region
        if (add.auto.rank.cut){
          rank.cuts <- c(rank.cuts,"autoSelect")
        }
        for (i in 1:length(rank.cuts)){
          rc <- rank.cuts[i]
          if (rc == "autoSelect") {
            rc <- rc.auto
            if (verbose) logger.info(c("Rank cutoff:",rc,"(auto-select)"))
          } else {
            rc <- as.integer(rc)
            if (verbose) logger.info(c("Rank cutoff:",rc))
          }
          is.dmr.hyper <- rrs.hyper <= rc
          is.dmr.hypo <- rrs.hypo <= rc
          gids.hyper <- unique(unlist(sapply(gene.ids.chr[is.dmr.hyper],FUN=function(idstr){
            unlist(strsplit(idstr, ";", fixed = TRUE))
          })))
          gids.hyper <- na.omit(gids.hyper)
          if (length(gids.hyper)>0 && length(uids)>0){
            dm.en.region.hyper <- performGOenrichment.diffMeth.entrez(gids.hyper,uids,ontology=oo,assembly=assembly(rnbSet),...)
          } else {
            dm.en.region.hyper <- NULL
          }
          
          gids.hypo <- unique(unlist(sapply(gene.ids.chr[is.dmr.hypo],FUN=function(idstr){
            unlist(strsplit(idstr, ";", fixed = TRUE))
          })))
          gids.hypo <- na.omit(gids.hypo)
          if (length(gids.hypo)>0 && length(uids)>0){
            dm.en.region.hypo <- performGOenrichment.diffMeth.entrez(gids.hypo,uids,ontology=oo,assembly=assembly(rnbSet),...)
          } else {
            dm.en.region.hypo <- NULL
          }
          dm.en.region.list <- c(dm.en.region.list,list(list(hyper=dm.en.region.hyper,hypo=dm.en.region.hypo)))
        }
        names(dm.en.region.list) <- paste("rankCut_",rank.cuts,sep="")
        dm.en.region.list.list <- c(dm.en.region.list.list,list(dm.en.region.list))
        names(dm.en.region.list.list)[length(dm.en.region.list.list)] <- rt
        logger.completed()
      }
      dv.enrich$region_var[[cc]][[oo]] <- dm.en.region.list.list
      logger.completed()
      logger.completed()
    }
    logger.completed()
  }
  class(dv.enrich) <- "DiffMeth.go.enrich"
  logger.completed()
  return(dv.enrich)
}

#' performLolaEnrichment.diffVar
#'
#' performs LOLA enrichment analysis for a given differential variability table.
#' @author Michael Scherer and Fabian Mueller
#' @param rnbSet RnBSet object for which differential variability was computed
#' @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param enrich.diffMeth Enrichment object as obtained from \code{\link{performLolaEnrichment.diffMeth}}. If it is not provided
#'          a new object is created.
#' @param lolaDbPaths LOLA database paths
#' @param rank.cuts.region Cutoffs for combined ranking that are used to determine differentially variable regions
#' @param add.auto.rank.cut flag indicating whether an automatically computed cut-off should also be considered.
#' @param rerank For deterimining differential variability: should the ranks be ranked again or should the absolute ranks be used.
#' @param verbose Enable for detailed status report
#' @return a DiffMeth.lola.enrich object (S3) containing the following attributes
#' \item{region}{Enrichment information for differential variability on the region level. A \code{data.table} object
#' as returned by the \code{runLOLA} function from the \code{LOLA} package for further details. Each element will contain different
#' user sets for different rank cutoffs and hyper/hypomethylation events(\code{userSet} column)}
#' \item{lolaDb}{The loaded \code{lolaDb} object containing the merged databases as returned by \code{\link{loadLolaDbs}}}
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.diffVar(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffVar(rnb.set.example,dm,lolaDirs[["hg19"]])
#' }
performLolaEnrichment.diffVar <- function(rnbSet, diffmeth, enrich.diffMeth=NULL, lolaDbPaths, rank.cuts.region=c(100,500,1000), add.auto.rank.cut=TRUE, rerank=TRUE, verbose=TRUE){
  nCores <- ifelse(parallel.isEnabled(), parallel.getNumWorkers(), 1)
  lolaDb <- loadLolaDbs(lolaDbPaths)
  comps <- get.comparisons(diffmeth)
  region.types <- get.region.types(diffmeth)
  diff.col.reg <- "mean.var.diff"
  if(!is.null(enrich.diffMeth)){
    enrich.diffMeth$region_var <- list()
    dv.lola.enrich <- enrich.diffMeth
    rm(enrich.diffMeth)
  }else{
    dv.lola.enrich <- list(region_var=list())
  }
  for (cc in comps){
    dv.lola.enrich$probe_var[[cc]] <- list()
    dv.lola.enrich$region_var[[cc]] <- list()
  }
  logger.start("Differential Variability LOLA Enrichment Analysis")
  for (rt in region.types){
    logger.start(c("Region Type:",rt))
    annot <- annotation(rnbSet, type=rt)
    # convert annotaiton to GRanges object
    univGr <- data.frame2GRanges(annot, chrom.column="Chromosome", start.column="Start", end.column="End", assembly=assembly(rnbSet), sort.result=FALSE)
    
    lolaUserSets <- list()
    
    for (i in 1:length(comps)){
      cc <- comps[i]
      compKey <- paste0("comp", i)
      logger.start(c("Comparison: ",cc))
      dmt <- get.table(diffmeth,cc,rt,undump=TRUE,return.data.frame=TRUE)
      
      if(!is.element("combinedRank.var",colnames(dmt))){
        stop("No differential variability information available. Be sure to execute differential variability module beforehand.")
      }
      rrs <- dmt[,"combinedRank.var"]
      rrs.hyper <- rrs
      rrs.hypo <- rrs
      rrs.hyper[dmt[,diff.col.reg] <= 0] <- NA
      rrs.hypo[dmt[,diff.col.reg] >= 0] <- NA
      #automatically select rank cutoff
      rc.auto <- 0L
      if (add.auto.rank.cut){
        rc.auto <- as.integer(auto.select.rank.cut(dmt$comb.p.adj.var.fdr,dmt$combinedRank.var))
      }
      
      if (rerank){
        rrs <- rank(rrs,na.last="keep",ties.method="min")
        rrs.hyper <- rank(rrs.hyper,na.last="keep",ties.method="min")
        rrs.hypo <- rank(rrs.hypo,na.last="keep",ties.method="min")
        if (add.auto.rank.cut && rc.auto > 0L){
          rc.auto <- na.omit(rrs[dmt[,"combinedRank.var"]==rc.auto])[1] #arbitrarily select the first rerank if the rank cut threshold is occupied by multiple ranks
        }
      }
      
      rank.cuts <- rank.cuts.region
      if (add.auto.rank.cut){
        rank.cuts <- c(rank.cuts,"autoSelect")
      }
      for (j in 1:length(rank.cuts)){
        rcn <- paste("rankCut_", rank.cuts[j] ,sep="")
        
        rc <- rank.cuts[j]
        if (rc == "autoSelect") {
          rc <- rc.auto
          if (verbose) logger.info(c("Rank cutoff:",rc,"(auto-select)"))
        } else {
          rc <- as.integer(rc)
          if (verbose) logger.info(c("Rank cutoff:",rc))
        }
        is.dmr.hyper <- rrs.hyper <= rc
        is.dmr.hyper[is.na(is.dmr.hyper)] <- FALSE
        is.dmr.hypo <- rrs.hypo <= rc
        is.dmr.hypo[is.na(is.dmr.hypo)] <- FALSE
        
        lolaUserSets <- c(lolaUserSets, list(univGr[is.dmr.hyper]), list(univGr[is.dmr.hypo]))
        N <- length(lolaUserSets)
        names(lolaUserSets)[c(N-1,N)] <- paste(compKey, rcn, c("hyper", "hypo"), sep="_")
      }
      logger.completed()
    }
    
    lolaUserSets <- GRangesList(lolaUserSets)
    logger.start("Running LOLA")
    lolaRes <- runLOLA(lolaUserSets, univGr, lolaDb, cores=nCores)
    logger.completed()
    
    for (i in 1:length(comps)){
      cc <- comps[i]
      curPat <- paste0("^comp",i,"_")
      curSubset <- grepl(curPat, lolaRes[["userSet"]])
      lolaRes.cur <- lolaRes[curSubset,]
      lolaRes.cur[["userSet"]] <- gsub(curPat, "", lolaRes.cur[["userSet"]])
      dv.lola.enrich$region_var[[cc]][[rt]] <- lolaRes.cur
    }
    
    logger.completed()
  }
  if(!is.element("lolaDb",names(dv.lola.enrich))){
    dv.lola.enrich[["lolaDb"]] <- lolaDb
  }
  class(dv.lola.enrich) <- "DiffMeth.lola.enrich"
  logger.completed()
  return(dv.lola.enrich)
}
