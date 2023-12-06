########################################################################################################################
## combining.R
## created: 2017-04-13
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions related to combining datasets of the same or different platforms, differential methylation results, etc.
########################################################################################################################

## TODO: Move the "combine" method here; use the function rnb.combine.arrays in it

## F U N C T I O N S ###################################################################################################

#' Combine dataset slots
#' 
#' Combines the data matrices of two datasets into a single one.
#' 
#' @param m1    Matrix of the first dataset. Set this to \code{NULL} if the first dataset does not contain this slot.
#' @param m2    Matrix of the second dataset. Set this to \code{NULL} if the second dataset does not contain this slot.
#' @param ii    Two-column \code{integer} \code{matrix} of indices in \code{m1} and \code{m2} to be considered for the
#'              combined matrix.
#' @param nn    Two-element \code{integer} \code{vector} listing the number of samples in the first and second dataset,
#'              respectively. If \code{m1} is provided, it should contain \code{nn[1]} columns.
#' @param useff Flag indicating if disk dumping for the output matrix must be enabled.
#' @return When \code{m1} and \code{m2} are both \code{NULL}, the return value is also \code{NULL}. Otherwise, the
#'         constructed combined matrix as an \code{ff} object (if \code{useff} is \code{TRUE}), or as an instance of
#'         \code{matrix} (if \code{useff} is \code{FALSE}). The dimensions of the matrix are \code{nrow(ii)} x
#'         \code{sum(nn)}.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.combine.matrices <- function(m1, m2, ii, nn, useff = rnb.getOption("disk.dump.big.matrices")) {
	if (is.null(m1)) {
		if (is.null(m2)) {
			return(NULL)
		}
		if ("ff" %in% class(m2)) {
			d.type <- vmode(m2)
		} else {
			d.type <- typeof(m2)
		}
	} else if ("ff" %in% class(m1)) {
		d.type <- vmode(m1)
	} else {
		d.type <- typeof(m1)
	}
	na.value <- do.call(paste0("as.", d.type), list(x = NA))
	if (useff) {
		mm <- ff(na.value, dim = c(nrow(ii), sum(nn)), dimnames = list(NULL, NULL), vmode = d.type)
	} else {
		mm <- matrix(na.value, nrow = nrow(ii), ncol = sum(nn))
	}
	if (!is.null(m1)) {
        present<-which(!is.na(ii[,1]))
		for (i in 1:(nn[1])) {
			mm[present, i] <- m1[ii[present, 1], i]
		}
	}
	if (!is.null(m2)) {
        present<-which(!is.na(ii[,2]))
		for (i in 1:(nn[2])) {
			mm[present, nn[1] + i] <- m2[ii[present, 2], i]
		}
	}
	mm
}

########################################################################################################################

#' Combine pheno tables
#' 
#' Combines two sample annotation tables.
#' 
#' @param dataset1 First input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @param dataset2 Second input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @return Combined sample annotation table as a \code{data.frame}.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.combine.pheno <- function(dataset1, dataset2) {
	tbl1 <- dataset1@pheno
	tbl2 <- dataset2@pheno
	if (!identical(sapply(tbl1, class), sapply(tbl2, class))) {
		stop("Mismatch in sample annotation tables")
	}
	for (i in which(sapply(tbl1, class) == "factor")) {
		if (!identical(levels(tbl1[, i]), levels(tbl2[, i]))) {
			stop("Mismatch in sample annotation tables")
		}
	}
	rbind(tbl1, tbl2)
}

########################################################################################################################

#' Combine array-based datasets
#'
#' Concatenates two array-based datasets focusing on the common probes.
#'
#' @param dataset1 First input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @param dataset2 Second input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @param type Type of the combine operation as a character singleton, one of "common", "all.x", "all.y" and "all".
#' @return Combined dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#'
#' @details \describe{
#'   \item{Sample annotation tables}{This method expects that the sample annotation tables of the two datasets have
#'      identical structures.}
#'   \item{Genome assembly}{This method expects that the two datasets target the same genome assembly.}
#'   \item{Platform}{The platform of the combined dataset is the most recent among the platforms of the input datasets.}
#'   \item{Intensity values}{The combined dataset is of type \code{\linkS4class{RnBeadRawSet}} only when both input
#'      datasets are of this type. Otherwise, any intensity value data is ignored.}
#'   \item{Probes}{Only the common probes are included in the resulting dataset.}
#'   \item{Regions}{Regions summarized in any of the input datasets are ignored. In the resulting dataset, regions are
#'      summarized as specified in the analysis option \code{"region.types"}.}
#'   \item{Quality control data}{QC data in the input datasets is ignored. The combined dataset includes no data on QC
#'      probe intensities.}
#'   \item{Infered covariates}{Inferred covariates in the input datasets are ignored. The combined dataset includes no
#'      data on inferred covariates.}
#'   \item{Disk dumping}{The combined dataset stores big tables on disk when the analysis option
#'      \code{"disk.dump.big.matrices"} is enabled.}
#' }
#' 
#' @author Yassen Assenov
#' @export
rnb.combine.arrays <- function(dataset1, dataset2, type="common") {
	if (!inherits(dataset1, "RnBeadSet")) {
		stop("Invalid value for dataset1")
	}
	if (!inherits(dataset2, "RnBeadSet")) {
		stop("Invalid value for dataset2")
	}
	if (dataset1@assembly != dataset2@assembly) {
		stop("Incompatible assemblies")
	}
	i <- c(dataset1@target, dataset2@target)
	common.platform <- c("probesEPICv2" = "EPICv2", "probesEPIC" = "EPIC", "probes450" = "450k", "probes27" = "27k") ## TODO: Validate for EPICv2 (on EPICv1 ?)
	if (!(i[1] %in% names(common.platform))) {
		stop("Unsupported platform for dataset1")
	}
	if (!(i[2] %in% names(common.platform))) {
		stop("Unsupported platform for dataset2")
	}
	common.platform <- common.platform[intersect(names(common.platform), unique(i))]
	common.platform <- unname(common.platform[1])
	is.raw <- (inherits(dataset1, "RnBeadRawSet") && inherits(dataset2, "RnBeadRawSet"))

	## Combine sample annotation tables
	#tbl <- rnb.combine.pheno(dataset1, dataset2)
    tbl<-plyr::rbind.fill(pheno(dataset1),pheno(dataset2))
	nn <- c(nrow(dataset1@pheno), nrow(dataset2@pheno))

	## Identify common probes
    if(type == "common"){
	    common.sites <- intersect(rownames(dataset1@sites), rownames(dataset2@sites))
    }else if(type == "all.x"){
        common.sites <- rownames(dataset1@sites)
    }else if(type == "all.y"){
        common.sites <- rownames(dataset2@sites)
    }else if(type == "all"){
        common.sites <- unique(c(rownames(dataset1@sites), rownames(dataset2@sites)))
    }else{
        rnb.error("Unsupported value for type")
    }
	if (length(common.sites) == 0) {
		stop("No common sites identified")
	}
#	ii <- cbind( # probe names must be sorted in both datasets!
#		which(rownames(dataset1@sites) %in% common.sites),
#		which(rownames(dataset2@sites) %in% common.sites))
    ii <- cbind(
        match(common.sites, rownames(dataset1@sites)),
        match(common.sites, rownames(dataset2@sites))
            )
	useff <- rnb.getOption("disk.dump.big.matrices")
	
	## Combine the data matrices
	result <- list(
		"pheno" = tbl,
		"probes" = common.sites,
		"platform" = common.platform,
		"region.types" = rnb.region.types.for.analysis(dataset1@assembly),
		"qc" = NULL,
		"useff" = useff,
		"ffcleanup" = rnb.getOption("enforce.destroy.disk.dumps"))
	slot.names <- c("p.values")
	if (is.raw) {
		slot.names <- c(slot.names, "M", "U", "M0", "U0", "bead.counts.M", "bead.counts.U")
	} else {
		result[["ffcleanup"]] <- NULL
		slot.names <- c(slot.names, "betas")
	}
	for (slot.name in slot.names) {
		if (slot.name == "p.values") {
			s.name <- "pval.sites"
		} else if (slot.name == "betas") {
			s.name <- "meth.sites"
		} else {
			s.name <- slot.name
		}
		result[[slot.name]] <- rnb.combine.matrices(slot(dataset1, s.name), slot(dataset2, s.name), ii, nn, useff)
	}
	
	## Construct the resulting object
	result <- do.call(ifelse(is.raw, "RnBeadRawSet", "RnBeadSet"), result)
	rm(tbl, i, common.platform, is.raw, nn, common.sites, ii, useff, slot.names, slot.name, s.name)

	## Set normalization and background subtraction methods
	if (dataset1@status$normalized == dataset2@status$normalized) {
		result@status$normalized <- dataset1@status$normalized
	} else {
		result@status$normalized <- "none"
		warning("Inconsistent normalization methods; setting to none")
	}
	if (dataset1@status$background == dataset2@status$background) {
		result@status$background <- dataset1@status$background
	} else {
		result@status$background <- "none"
		warning("Inconsistent background subtraction methods; setting to none")
	}
	result
}
###############################################################################
#' rnb.combine.seq
#'
#' Initial implementation of the combine method for sequencing datasets.
#' 
#' @param x An object of type \code{\link{RnBiseqSet-class}} used for concatenation
#' @param y Another object of type \code{\link{RnBiseqSet-class}} used for concatenation
#' @param type A character representing the type of combination. Needs to be one of \code{"common"},
#'          \code{all.x}, \code{all.y} or \code{all}.
#' @return An \code{\link{RnBiseqSet-class}} object with combined information
#' @details The type parameters determines the mode of combination:
#'           \itemize{
#'             \item{\code{"common"}}{ The intersection between the sites present in the two 
#'                             datasets is used for the new dataset.}
#'             \item{\code{"all.x"}}{ All sites present in x are used.}
#'             \item{\code{"all.y"}}{ All sites present in y are used.}
#'             \item{\code{"all"}}{ The union between the sites of both datasets is used.}
#'           }
#'
#' @export
rnb.combine.seq<-function(x,y,type="common"){
#    if (class(x) != class(y)){
#        stop("Could not combine RnBiseqSet objects: incompatible classes")
#    }
    if (assembly(x) != assembly(y)){
        stop("Could not combine RnBiseqSet objects: incompatible assemblies")
    }
#    if (x@target != y@target){
#        stop("Could not combine RnBiseqSet objects: incompatible platforms")
#    }
    common.samples <- intersect(samples(x),samples(y))
    if (length(common.samples)>0){
        stop(paste0("Could not combine RnBSet objects: the following samples overlap in both objects: ", paste(common.samples,collapse=",")))
    }
    useff <- x@status$disk.dump
    if (x@status$disk.dump != y@status$disk.dump){
        warning(paste0("disk dump status of the two objects to combine disagree. Using disk dump: ", useff))
    }
    usebigff <- useff
    if (usebigff) usebigff <- !is.null(x@status$disk.dump.bigff)
    if (usebigff) usebigff <- x@status$disk.dump.bigff
    if(usebigff){
        bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")
    }
    # prepare a new object
    new.set<-x
    # if(nrow(pheno(x))>=nrow(pheno(y))){
    # 	new.set<-y
    # }else{
    # 	new.set<-x
    # }
    
    new.set@pheno <- plyr::rbind.fill(pheno(x),pheno(y))
    
    # combine sites
    sites1<-x@sites
    sites2<-y@sites
    
    ## Identify common probes
    if(type == "common"){
        common.chr<-intersect(unique(sites1[,2]), unique(sites2[,2]))
    }else if(type == "all.x"){
        common.chr<-unique(sites1[,2])
    }else if(type == "all.y"){
        common.chr<-unique(sites2[,2])
    }else if(type == "all"){
        common.chr<-union(unique(sites1[,2]), unique(sites2[,2]))
    }else{
        rnb.error("Unsupported value for type")
    }
    
    subs1<-list()
    subs2<-list()
    map1<-list()
    map2<-list()
    common.sites<-list()
    
    for(chr in common.chr){
        if(type == "common"){
            sts<-sort(intersect(sites1[sites1[,2]==chr,3],sites2[sites2[,2]==chr,3]))
        }else if(type == "all.x"){
            sts<-sites1[sites1[,2]==chr,3]
        }else if(type == "all.y"){
            sts<-sites2[sites2[,2]==chr,3]
        }else if(type == "all"){
            sts<-sort(union(sites1[sites1[,2]==chr,3],sites2[sites2[,2]==chr,3]))
        }
        mtch1<-match(sites1[sites1[,2]==chr,3], sts)
        valid1<-which(!is.na(mtch1))
        map1[[chr]]<-valid1
        subs1[[chr]]<-mtch1[valid1]
        mtch2<-match(sites2[sites2[,2]==chr,3], sts)
        valid2<-which(!is.na(mtch2))
        map2[[chr]]<-valid2
        subs2[[chr]]<-mtch2[valid2]
        common.sites[[chr]]<-cbind(rep(1,length(sts)), rep(chr,length(sts)), sts)
    }
    
    total.sites<-sum(sapply(common.sites, nrow))
    
    if("ff_matrix" %in% c(class(sites1), class(sites2))){
        new.sites <- ff(vmode="integer", dim=c(total.sites,3))
        ixx<-1
        for(sts in common.sites){
            new.sites[ixx:(ixx+nrow(sts)),]<-sts
            ixx<-ixx+nrow(sts)+1
        }
    }else{
        new.sites<-do.call("rbind", common.sites)
    }
    
    colnames(new.sites)<-NULL
    new.set@sites<-new.sites
    
    slot.names<-RNBSET.SLOTNAMES
    
    if(inherits(x, "RnBeadSet")){
        slot.names<-c(slot.names, RNBEADSET.SLOTNAMES)
    }
    if(inherits(x, "RnBeadRawSet")){
        slot.names<-c(slot.names, RNBRAWSET.SLOTNAMES)
    }
    
    for(sl in slot.names){
        if(all(!is.null(slot(x,sl)),!is.null(slot(y,sl)))){
            if(useff){
                #new.matrix<-ff(vmode=vmode(slot(x,sl)), dim=c(total.sites,nrow(pheno(new.set))))
                if (usebigff){
                    new.matrix <- BigFfMat(row.n=total.sites, col.n=nrow(pheno(new.set)), vmode=vmode(slot(x,sl)), finalizer=bff.finalizer)
                } else {
                    new.matrix <- create.empty.ff.matrix.tmp(vm=vmode(slot(x,sl)), dim=c(total.sites,nrow(pheno(new.set))))
                }
            }else{
                new.matrix<-matrix(NA, nrow=total.sites, ncol=nrow(pheno(new.set)))
            }
            for(chr in common.chr){
                #new.matrix[new.sites[,2]==chr,1:nrow(pheno(x))]<-slot(x,sl)[sites1[map1[[chr]],2]==chr,][subs1[[chr]],]
                #new.matrix[new.sites[,2]==chr,(nrow(pheno(x))+1):nrow(pheno(new.set))]<-slot(y,sl)[sites2[,2]==chr,][subs2[[chr]],]
                ix<-which(new.sites[,2]==chr)
                new.matrix[ix[subs1[[chr]]],1:nrow(pheno(x))]<-slot(x,sl)[which(sites1[,2]==chr)[map1[[chr]]],,drop=FALSE]
                new.matrix[ix[subs2[[chr]]],(nrow(pheno(x))+1):nrow(pheno(new.set))]<-slot(y,sl)[which(sites2[,2]==chr)[map2[[chr]]],,drop=FALSE]
            }
            #colnames(new.matrix)<-c(colnames(slot(x,sl)), colnames(slot(y,sl)))
            slot(new.set, sl)<-new.matrix
            rm(new.matrix)
            rnb.cleanMem()
        }else{
            slot(new.set, sl)<-NULL
        }
        
        if(x@status$disk.dump && isTRUE(x@status$discard.ff.matrices)){
            delete(slot(x, sl))
        }
        if(y@status$disk.dump && isTRUE(y@status$discard.ff.matrices)){
            delete(slot(y, sl))
        }
    }
    
    if(inherits(x,"RnBeadSet")){
        if(all(!is.null(qc(x)),!is.null(qc(y)))){
            cpn<-intersect(rownames(qc(x)$Cy3), rownames(qc(x)$Cy3))
            cy3.new<-cbind(qc(x)$Cy3[cpn,], qc(y)$Cy3[cpn,])
            cy5.new<-cbind(qc(x)$Cy5[cpn,], qc(y)$Cy5[cpn,])
            colnames(cy3.new)<-NULL
            colnames(cy5.new)<-NULL
            new.set@qc<-list(Cy3=cy3.new, Cy5=cy5.new)
        }else{
            new.set@qc<-NULL
        }
    }
    
    new.set@status<-list()
    if(inherits(new.set, "RnBeadSet")){
        new.set@status$normalized<-"none"
        new.set@status$background<-"none"
    }
    new.set@status$disk.dump<-useff
    new.set@status$disk.dump.bigff<-usebigff
    
    for (region.type in union(summarized.regions(x), summarized.regions(y))) {
        if (region.type %in% rnb.region.types(assembly(new.set))) {
            new.set <- summarize.regions(new.set, region.type)
        }
    }
    new.set@inferred.covariates<-list()
    
    rm(common.sites, sites1, sites2, subs1, subs2, x, y)
    rnb.cleanMem()
    new.set
}

###############################################################################
# basic_combine
#
# Initial implementation of the combine method
#
#
basic_combine<-function(x,y){
    if (class(x) != class(y)){
        stop("Could not combine RnBSet objects: incompatible classes")
    }
    if (assembly(x) != assembly(y)){
        stop("Could not combine RnBSet objects: incompatible assemblies")
    }
    if (x@target != y@target){
        stop("Could not combine RnBSet objects: incompatible platforms")
    }
    common.samples <- intersect(samples(x),samples(y))
    if (length(common.samples)>0){
        stop(paste0("Could not combine RnBSet objects: the following samples overlap in both objects: ", paste(common.samples,collapse=",")))
    }
    useff <- x@status$disk.dump
    if (x@status$disk.dump != y@status$disk.dump){
        warning(paste0("disk dump status of the two objects to combine disagree. Using disk dump: ", useff))
    }
    usebigff <- useff
    if (usebigff) usebigff <- !is.null(x@status$disk.dump.bigff)
    if (usebigff) usebigff <- x@status$disk.dump.bigff
    if(usebigff){
        bff.finalizer <- rnb.getOption("disk.dump.bigff.finalizer")
    }
    # prepare a new object
    new.set<-x
    # if(nrow(pheno(x))>=nrow(pheno(y))){
    # 	new.set<-y
    # }else{
    # 	new.set<-x
    # }
    
    new.set@pheno <- plyr::rbind.fill(pheno(x),pheno(y))
    
    # combine sites
    sites1<-x@sites
    sites2<-y@sites
    
    common.chr<-union(unique(sites1[,2]), unique(sites2[,2]))
    
    subs1<-list()
    subs2<-list()
    common.sites<-list()
    
    for(chr in common.chr){
        sts<-sort(union(sites1[sites1[,2]==chr,3],sites2[sites2[,2]==chr,3]))
        subs1[[chr]]<-match(sites1[sites1[,2]==chr,3], sts)
        subs2[[chr]]<-match(sites2[sites2[,2]==chr,3], sts)
        common.sites[[chr]]<-cbind(rep(1,length(sts)), rep(chr,length(sts)), sts)
    }
    
    total.sites<-sum(sapply(common.sites, nrow))
    
    if("ff_matrix" %in% c(class(sites1), class(sites2))){
        new.sites <- ff(vmode="integer", dim=c(total.sites,3))
        ixx<-1
        for(sts in common.sites){
            new.sites[ixx:(ixx+nrow(sts)),]<-sts
            ixx<-ixx+nrow(sts)+1
        }
    }else{
        new.sites<-do.call("rbind", common.sites)
    }
    
    colnames(new.sites)<-NULL
    new.set@sites<-new.sites
    
    slot.names<-RNBSET.SLOTNAMES
    
    if(inherits(x, "RnBeadSet")){
        slot.names<-c(slot.names, RNBEADSET.SLOTNAMES)
    }
    if(inherits(x, "RnBeadRawSet")){
        slot.names<-c(slot.names, RNBRAWSET.SLOTNAMES)
    }
    
    for(sl in slot.names){
        if(all(!is.null(slot(x,sl)),!is.null(slot(y,sl)))){
            if(useff){
                #new.matrix<-ff(vmode=vmode(slot(x,sl)), dim=c(total.sites,nrow(pheno(new.set))))
                if (usebigff){
                    new.matrix <- BigFfMat(row.n=total.sites, col.n=nrow(pheno(new.set)), vmode=vmode(slot(x,sl)), finalizer=bff.finalizer)
                } else {
                    new.matrix <- create.empty.ff.matrix.tmp(vm=vmode(slot(x,sl)), dim=c(total.sites,nrow(pheno(new.set))))
                }
            }else{
                new.matrix<-matrix(NA, nrow=total.sites, ncol=nrow(pheno(new.set)))
            }
            for(chr in common.chr){
                #new.matrix[new.sites[,2]==chr,1:nrow(pheno(x))]<-slot(x,sl)[sites1[,2]==chr,][subs1[[chr]],]
                #new.matrix[new.sites[,2]==chr,(nrow(pheno(x))+1):nrow(pheno(new.set))]<-slot(y,sl)[sites2[,2]==chr,][subs2[[chr]],]
                ix<-which(new.sites[,2]==chr)
                new.matrix[ix[subs1[[chr]]],1:nrow(pheno(x))]<-slot(x,sl)[sites1[,2]==chr,,drop=FALSE]
                new.matrix[ix[subs2[[chr]]],(nrow(pheno(x))+1):nrow(pheno(new.set))]<-slot(y,sl)[sites2[,2]==chr,,drop=FALSE]
            }
            #colnames(new.matrix)<-c(colnames(slot(x,sl)), colnames(slot(y,sl)))
            slot(new.set, sl)<-new.matrix
            rm(new.matrix)
            rnb.cleanMem()
        }else{
            slot(new.set, sl)<-NULL
        }
        
        if(x@status$disk.dump && isTRUE(x@status$discard.ff.matrices)){
            delete(slot(x, sl))
        }
        if(y@status$disk.dump && isTRUE(y@status$discard.ff.matrices)){
            delete(slot(y, sl))
        }
    }
    
    if(inherits(x,"RnBeadSet")){
        if(all(!is.null(qc(x)),!is.null(qc(y)))){
            cpn<-intersect(rownames(qc(x)$Cy3), rownames(qc(x)$Cy3))
            cy3.new<-cbind(qc(x)$Cy3[cpn,], qc(y)$Cy3[cpn,])
            cy5.new<-cbind(qc(x)$Cy5[cpn,], qc(y)$Cy5[cpn,])
            colnames(cy3.new)<-NULL
            colnames(cy5.new)<-NULL
            new.set@qc<-list(Cy3=cy3.new, Cy5=cy5.new)
        }else{
            new.set@qc<-NULL
        }
    }
    
    new.set@status<-list()
    if(inherits(new.set, "RnBeadSet")){
        new.set@status$normalized<-"none"
        new.set@status$background<-"none"
    }
    new.set@status$disk.dump<-useff
    new.set@status$disk.dump.bigff<-usebigff
    
    for (region.type in union(summarized.regions(x), summarized.regions(y))) {
        if (region.type %in% rnb.region.types(assembly(new.set))) {
            new.set <- summarize.regions(new.set, region.type)
        }
    }
    new.set@inferred.covariates<-list()
    
    rm(common.sites, sites1, sites2, subs1, subs2, x, y)
    rnb.cleanMem()
    new.set
}