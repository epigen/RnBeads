################################################################################

#' rnb.execute.pOOBAH
#' 
#' Probe signal intensities are masked based on their out-of-band signal intensities to counter hybridization failure.
#' 
#' rnb.execute.pOOBAH is used to apply the method \emph{pOOBAH} (P-value with OOB probes for Array Hybridization), 
#' which was conceived by Zhou, Triche, Laird and Shen to mask probes associated with hybridization failures. 
#' pOOBAH has been implemented in the R-package \emph{"sesame"}, a dependency needed for this function (see Zhou et al, 2018 and the respective Bioconductor/github pages). 
#' \emph{pOOBAH} computes the detection p-values by constructing 2 empirical cumulative density functions (eCDFs) based 
#' on the out-of-band signal intensities of the red and the green channel, respectively, to detect hybridization failures. 
#' The (in-band) green and red channel signal intensities of the probes are passed to the eCDFs and the probes with a 
#' p-value higher than the given threshold (\code{pval.thresh}) are masked, as they are considered background. 
#' \emph{pOOBAH} is applied \emph{separately to each sample}. 
#' Hybridization failures might occur due to somatic or germline deletions. 
#' In addition, unreliable low-intensity probes might also be masked.
#'  
#' @param raw.set Methylation dataset as an instance of \code{RnBeadRawSet}.
#' @param anno.table Annotation for \code{raw.set}.
#' @param pval.thresh Computed detection p-values above this threshold are masked. Default value is 0.05.
#' @param verbose If set to true, a short information is printed on how many probes are masked by the method. 
#' @return Returns a modified \code{RnBeadRawSet}, in which signal intensities are masked, if their computed p-value  
#'         was greater than \code{pval.thresh}. Note, in datasets with several samples, signal intensities of a specific probe 
#'         might be masked in sample A, but not in sample B, as \emph{pOOBAH} is applied separately to each sample.
#'         For example: the signal intensities of probe cg24488772 might be masked in sample 1, but not in sample 12. 
#' @author \emph{pOOBAH} method: Wanding Zhou. Adapted by Nathan Steenbuck. 
#' @examples
#' library(RnBeads.hg19) 
#' data(small.example.object)
#' 
#'
#' @export 

rnb.execute.pOOBAH <- function(raw.set, anno.table = NULL, pval.thresh = 0.05, verbose = FALSE){
  rnb.require("sesame")
  
  if(!(is.numeric(pval.thresh) && pval.thresh <= 1 && pval.thresh >= 0)){
    stop("Invalid value for pval.thresh. Please specify a numeric in the range of [0, 1].")
  }
  if(is.numeric(anno.table)){
    stop("Invalid value for anno.table. Wanted to specify the p-value threshold?")
  }
  if(!inherits(raw.set, "RnBeadRawSet")){
    stop("Please provide input inhereting from RnBeadRawSet.")
  }
  else{
    if(raw.set@target == "probes450"){
      platform = "HM450"
    }else if(raw.set@target == "probesEPIC"){
      platform = "EPIC"
    }else if(raw.set@target == "probesEPICv2"){
      platform = "EPICv2"
    }else if(raw.set@target == "probesMSA"){
      platform = "MSA"
    }else{
      stop("Invalid value for platform")
    }

    if(is.null(anno.table) || !("ID" %in% colnames(anno.table))){
      anno.table <- annotation(raw.set, add.names=TRUE) 
    }else if(length(anno.table[["ID"]]) != nrow(raw.set@sites)){
      stop("The annotation and dataset are not compatible.")
    }
    probeIDs <- anno.table[["ID"]]
    
    intensities.by.channel <- intensities.by.color(raw.set, address.rownames = FALSE, add.oob = TRUE, 
                                                  add.controls =  FALSE, add.missing = FALSE, re.separate = TRUE)
    
    grn <- intensities.by.channel$Cy3.I
    red <- intensities.by.channel$Cy5.I
    grn.oob <- intensities.by.channel$Cy3.I.oob
    red.oob <- intensities.by.channel$Cy5.I.oob
    tII <- intensities.by.channel$II
    rm(intensities.by.channel, anno.table)
    
    if(sum((rownames(grn$M) %in% rownames(grn$U)) == FALSE) > 0 ||
       sum((rownames(red$M) %in% rownames(red$U)) == FALSE) > 0 ||
       sum((rownames(grn.oob$M) %in% rownames(grn.oob$U)) == FALSE) > 0 ||
       sum((rownames(red.oob$M) %in% rownames(red.oob$U)) == FALSE) > 0 ||
       sum((rownames(tII$M) %in% rownames(tII$U)) == FALSE) > 0 ||
       length(grn$M) != length(red.oob$M) ||
       length(red$M) != length(grn.oob$M)){
      stop("Equal dimensions and IDs are expected.")
    }
    
    nsamples = length(samples(raw.set)) 
    if(nsamples == 0){
      stop("Dataset contains no samples.")
    }
    
    if(is.null(raw.set@pval.sites) || nrow(raw.set@pval.sites) != length(probeIDs) 
       || ncol(raw.set@pval.sites) != nsamples){
      raw.set@pval.sites <- matrix(data = NA, nrow = length(probeIDs), ncol = nsamples)
    }
    
    nmasked = 0
    
    for (i in 1:nsamples){
      sdf = SigDF(rbind(
          data.frame(
              Probe_ID = rownames(grn$M),
              MG = grn$M[,i], MR = red.oob$M[,i],
              UG = grn$U[,i], UR = red.oob$U[,i], col="G", mask=FALSE),
          data.frame(
              Probe_ID = rownames(red$M),
              MG = grn.oob$M[,i], MR = red$M[,i],
              UG = grn.oob$U[,i], UR = red$U[,i], col="R", mask=FALSE),
          data.frame(
              Probe_ID = rownames(tII$M),
              MG = NA, MR = NA,
              UG = tII$M[,i], UR = tII$U[,i], col="2", mask=FALSE)), platform)
      
      pvalues = pOOBAH(sdf, return.pval = TRUE)

      mask <- names(pvalues)[pvalues > pval.thresh]
      mask <- mask[!is.na(mask)]
      raw.set@pval.sites[, i] <- pvalues[match(probeIDs, names(pvalues))] 
      
      if(!(length(mask) == 0)){
        nmasked = nmasked + length(mask)
        maskedIDs <- match(mask, probeIDs) 
        raw.set@pval.sites[maskedIDs, i] <- NA
        raw.set@M[maskedIDs, i] <- NA
        raw.set@U[maskedIDs, i] <- NA
        raw.set@M0[maskedIDs, i]<- NA
        raw.set@U0[maskedIDs, i]<- NA
      }
    }
    rm(sigset.l, maskedIDs, mask, pvalues, grn, grn.oob, red, red.oob, tII, pIG, pII, pIR, platform, i)
    if (verbose){
      nprobes = length(probeIDs)    
      ntotal = nsamples*nprobes
      
      logger.info(paste('=======================',
        '=    pOOBAH           =',
        '=======================',
        paste0('No. probes: ', nprobes),
        paste0('No. samples: ', nsamples),
        paste0('No. probes times samples: ', ntotal),
        paste0('No. of masked probes: ', nmasked),
        paste0('Fraction of masked probes: ', round(nmasked/ntotal, digits = 3)),sep='\n'))
    }
    
    return(raw.set)
  }
}

################################################################################
