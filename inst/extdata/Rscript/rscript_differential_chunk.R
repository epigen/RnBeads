library(argparse)
suppressPackageStartupMessages(library(RnBeads))

ap <- ArgumentParser()
ap$add_argument("-x", "--xml", action="store", help="Configuration xml file")
ap$add_argument("-s", "--rnbSet", action="store", dest="rnbSetLoc", help="Location of stored RnBSet")
ap$add_argument("-i", "--id", action="store", help="An analysis id that will be used to name the output")
ap$add_argument("-o", "--output", action="store", help="Output directory")
ap$add_argument("-c", "--cores", action="store", type="integer", default=1, help="Number of cores used for the analysis")
ap$add_argument("-r", "--region", action="append", dest="regionTypes", help="region type")
ap$add_argument("-p", "--pheno", action="append", dest="comparisons", help="Phenotype/comparison column")
cmdArgs <- ap$parse_args()
module.name <- "differential_chunk"

logger.start(fname=NA)
logger.status(c("...Started module:",module.name))

logger.start("Configuring Analysis")
	rnb.settings <- rnb.xml2options(cmdArgs$xml,return.full.structure=TRUE)

	report.dir <- rnb.settings$analysis.params[["dir.reports"]]
	analysis.options <- rnb.settings$options

	cmp.cols <- cmdArgs$comparisons
	reg.types <- cmdArgs$regionTypes
	chunk.id <- cmdArgs$id

	if ("preanalysis.script" %in% names(rnb.settings)){
		source(rnb.settings$preanalysis.script)
	} 
	## Set options
	if (length(analysis.options) != 0) {
		do.call(rnb.options, analysis.options)
	}

	if(rnb.getOption("region.subsegments")>1L){
		reg.types.subseg <- grep("\\.subsegments$",reg.types,value=TRUE)
		if (length(reg.types.subseg)>0){
			subseg.annot.dir <- file.path(report.dir,"preprocessing_data","subsegments")
			logger.start("Loading region sub-segmentation annotation")
			for (rt in reg.types.subseg){
				logger.status(c("Loading annotation",rt))
				fn <- file.path(subseg.annot.dir,paste0("regions_",rt,".RData"))
				rnb.load.annotation(fn, rt)
			}
			logger.completed()
		}
	}

	logger.machine.name()

	if (cmdArgs$cores > 1) {
		parallel.setup(cmdArgs$cores)
	}

	aname <- rnb.getOption("analysis.name")
	if (!(is.null(aname) || is.na(aname) || nchar(aname) == 0)) {
		logger.info(c("Analysis Title:", aname))
	}
	ncores <- parallel.getNumWorkers()
	if (ncores == -1) {
		ncores <- 1L
	}
	logger.info(c("Number of cores:", ncores))
	rm(aname, ncores)

	lolaDbPaths <- NULL
	savedLolaDbPath <- file.path(cmdArgs$output, "lolaDbPaths.rds")
	if(rnb.getOption("differential.enrichment.lola")){
		if (file.exists(savedLolaDbPath)){
			lolaDbPaths <- readRDS(savedLolaDbPath)
			logger.info(c("Successfully loaded LOLA DB paths"))
		} else{
			logger.info(c("No LOLA DB file found"))
		}
		
	}
logger.completed()

logger.start("Loading RnBSet")
	rnb.set <- load.rnb.set(cmdArgs$rnbSetLoc)
logger.completed()

logger.start(fname=c(file.path(report.dir,paste0("analysis_",module.name,"_",chunk.id,".log")),NA))

################################################################################
# main script
################################################################################

permutations <- rnb.getOption("differential.permutations")

logger.start("Differential Methylation")
	logger.info(c("Using",permutations,"permutation tests"))

	if (is.null(cmp.cols)){
		cmp.cols <- colnames(pheno(rnb.set))
	}
	logger.info(c("Using columns:",paste(cmp.cols,collapse=",")))

	if (length(reg.types)>0){
		logger.info(c("Using region types:",paste(reg.types,collapse=",")))
	} else {
		logger.info(c("Skipping region level analysis"))
	}
	disk.dump <- rnb.getOption("disk.dump.big.matrices")
	diffmeth <- rnb.execute.computeDiffMeth(rnb.set,cmp.cols,region.types=reg.types,n.perm=permutations,covg.thres=rnb.getOption("filtering.coverage.threshold"),
			pheno.cols.all.pairwise=rnb.getOption("differential.comparison.columns.all.pairwise"),columns.pairs=rnb.getOption("columns.pairing"),
			skip.sites=!rnb.getOption("analyze.sites"),
			# disk.dump=disk.dump,disk.dump.dir=paste0(cmdArgs$output,"_diffMethTableDir"))
			disk.dump=disk.dump,disk.dump.dir=paste0(tempfile(pattern=""),"_diffMethTableDir"))
	if (rnb.getOption("differential.enrichment.go")){
		dm.go.enrich <- performGoEnrichment.diffMeth(rnb.set,diffmeth,verbose=FALSE)
		if(rnb.getOption("differential.variability")){
		  dm.go.enrich <- performGOEnrichment.diffVar(rnb.set,diffmeth,dm.go.enrich,verbose=FALSE)
		}
	} else {
		dm.go.enrich <- NULL
		logger.info(c("Skipping GO enrichment analysis of differentially methylated regions"))
	}
	if (rnb.getOption("differential.enrichment.lola") && !is.null(lolaDbPaths)){
		dm.lola.enrich <- performLolaEnrichment.diffMeth(rnb.set,diffmeth,lolaDbPaths,verbose=FALSE)
		if(rnb.getOption("differential.variability")){
		  dm.lola.enrich <- performLolaEnrichment.diffVar(rnb.set,diffmeth,enrich.diffMeth = dm.lola.enrich, lolaDbPaths,verbose=FALSE)
		}
	} else {
		dm.lola.enrich <- NULL
		logger.info(c("Skipping LOLA enrichment analysis of differentially methylated regions"))
	}


	logger.start("Saving")
		diffmeth.path <- file.path(cmdArgs$output,paste0(module.name,"_",chunk.id,"_rnbDiffMeth"))
		save.rnb.diffmeth(diffmeth, diffmeth.path)
		diffmeth.go.enrichment <- dm.go.enrich
		if (!is.null(diffmeth.go.enrichment)){
			save(diffmeth.go.enrichment, file=file.path(diffmeth.path, "enrichment_go.RData"))
		}
		diffmeth.lola.enrichment <- dm.lola.enrich
		if (!is.null(diffmeth.lola.enrichment)){
			save(diffmeth.lola.enrichment, file=file.path(diffmeth.path, "enrichment_lola.RData"))
		}
	logger.completed()
logger.completed()

logger.start("Check result")
#check whether the output regions match the annotation regions
comps <- get.comparisons(diffmeth)
rts <- reg.types
if (rnb.getOption("analyze.sites")) rts <- c("sites",rts)
for (rr in rts){
	aa <- annotation(rnb.set,rr)
	n.regs <- nrow(aa)
	for (cc in comps){
		dmt <- get.table(diffmeth,cc,rr,return.data.frame=TRUE)
		if (n.regs != nrow(dmt)){
			msg <- paste0("Number of regions in the DMR table does not match the number in the annotation table [",cc,"--",rr,"]. This could be a memory issue. Try assigning more main memory to the task or reducing parallelity.")
			logger.error(msg)
		}
	}
}
logger.completed()

logger.status(c("...Completed module:",module.name))
quit(save='no')
