########################################################################################################################
## differentialvariability.R
## created: 2017-07-28
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## The differential variability analysis methods between sample groups.
########################################################################################################################

#' rnb.execute.diffVar
#' 
#' This routine computes sites that are differentially variable between two sample groups specified as the column name 
#' in the phenotypic table.
#' 
#' @author Michael Scherer
#' @export
#' @param rnb.set Object of type \code{\linkS4class{RnBSet}} on which differential variability analysis should be conducted
#' @param pheno.cols Column names used to define the classes, whose methylation variability should be compared with each other
#' @param region.types Regions types to be used for the analysis. Defaults to the results given by rnb.region.types.for.analysis of the given RnBSet.
#' @param columns.adj Column names or indices in the table of phenotypic information to be used for confounder adjustment in the differential variability analysis.
#' @param adjust.celltype Flag indicating whether the resulting table should also contain estimated celltype contributions. See \code{\link{rnb.execute.ct.estimation}} for details.
#' @param disk.dump Flag indicating whether the resulting differential methylation object should be file backed, ie.e the matrices dumped to disk
#' @param disk.dump.dir disk location for file backing of the resulting differential methylation object. Only meaningful if \code{disk.dump=TRUE}.
#' @return Object of type \code{\linkS4class{RnBDiffMeth}} containing information about the differential variability analysis.
rnb.execute.diffVar <- function(rnb.set,pheno.cols=rnb.getOption("differential.comparison.columns"),
                                region.types=rnb.region.types.for.analysis(rnb.set),
                                columns.adj=rnb.getOption("covariate.adjustment.columns"),
                                adjust.celltype=rnb.getOption("differential.adjustment.celltype"),
                                disk.dump=rnb.getOption("disk.dump.big.matrices"),
                                disk.dump.dir=tempfile(pattern="diffMethTables_")){
  logger.start("Differential Variability")
  logger.start("Retrieving comparison info")
  cmp.info <- get.comparison.info(rnb.set, pheno.cols=pheno.cols,region.types = region.types, 
                                  columns.adj=columns.adj,
                                  adjust.celltype=adjust.celltype)
  logger.completed()
  if (is.null(cmp.info)) {
    return(NULL)
  }
  
  variability.method <- rnb.getOption("differential.variability.method")
  diff.meth <- new("RnBDiffMeth",variability.method=variability.method,disk.dump=disk.dump,disk.path=disk.dump.dir)
  if(!isImputed(rnb.set)) rnb.set <- rnb.execute.imputation(rnb.set,update.ff=FALSE)

  for(i in 1:length(cmp.info)){
    cmp.info.cur <- cmp.info[[i]]
    logger.start(c("Comparing ",cmp.info.cur$comparison))
    if(cmp.info.cur$paired){
      logger.status("Conducting PAIRED analysis")
    }
    meth.matrix <- meth(rnb.set)
    diffVar <- computeDiffVar.bin.site(X=meth.matrix,inds.g1=cmp.info.cur$group.inds$group1,
                                       inds.g2=cmp.info.cur$group.inds$group2, 
                                       paired=cmp.info.cur$paired,adjustment.table=cmp.info.cur$adjustment.table)

    diff.meth <- addDiffMethTable(diff.meth,diffVar,comparison=cmp.info.cur$comparison,region.type="sites",grp.labs=cmp.info.cur$group.names)
    rnb.cleanMem()
    if (length(cmp.info.cur$region.types)>0){
      diffVar.region <- computeDiffVar.bin.region(rnb.set,diffVar,
                                                  inds.g1=cmp.info.cur$group.inds$group1,inds.g2=cmp.info.cur$group.inds$group2,
                                                  region.types=cmp.info.cur$region.types,
                                                  paired = cmp.info.cur$paired
      )
      for(rt in region.types){
        diff.meth <- addDiffMethTable(diff.meth,diffVar.region[[rt]],comparison=cmp.info.cur$comparison, 
                                      region.type=rt, grp.labs=cmp.info.cur$group.names
        )
      }
    }
    logger.completed()
  }
  if(!all(get.comparisons(diff.meth)%in%names(cmp.info)))
    diff.meth <- addComparisonInfo(diff.meth,cmp.info)
  logger.completed()
  return(diff.meth)
}

#' diffVar
#' 
#' This routine applies the diffVar method from the \code{missMethyl} package that determines sites exhibiting
#' differential variability between two sample groups
#' 
#' @author Michael Scherer
#' @param meth.matrix Matrix containing the methylation information used to calculate differentially variable sites 
#'                    between the two groups
#' @param inds.g1 Indices in the phenotypic table corresponding to the first group.
#' @param inds.g2 Indices in the phenotypic table corresponding to the second group.
#' @param adjustment.table A \code{data.frame} containing variables to adjust for in the testing
#' @param paired Should the analysis be performed in a paired fashion. If yes, the first index in \code{inds.g1} must 
#'           correspond to the first in \code{inds.g2} and so on.
#' @return P-values as the result of the diffVar method not adjusted for multiple hypothesis testing.
#' 
#' @export
#' @references Phipson, Belinda, Oshlack, Alicia (2014)
#'             DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging
#'             Genome Biology 15(9):465.
diffVar <- function(meth.matrix,inds.g1,inds.g2,adjustment.table=NULL,paired=FALSE){
  logger.start("diffVar method")
  rnb.require('missMethyl')
  if(is.logical(inds.g1)) inds.g1 <- which(inds.g1)
  if(is.logical(inds.g2)) inds.g2 <- which(inds.g2)
  n.g1 <- length(inds.g1)
  n.g2 <- length(inds.g2)
  n <- n.g1 + n.g2
  if (!is.null(adjustment.table)){
    if (!(is.data.frame(adjustment.table) && nrow(adjustment.table)==n && (!any(is.na(adjustment.table))))) {
      stop("invalid value for adjustment.table")
    }
    m <- ncol(adjustment.table)
    if (m == 0) {
      adjustment.table <- NULL
    } else {
      colnames(adjustment.table) <- paste0("x",1:m,"x")
    }
  }
  
  ind.vec <- c(inds.g1,inds.g2)
  if (length(ind.vec) < 2) stop("need at least two samples indices to compare")
  meth.matrix <- rnb.beta2mval(meth.matrix[,ind.vec,drop=FALSE])
  df <- data.frame(xg = factor(rep(c("group1","group2"), c(n.g1,n.g2)), levels=c("group1","group2")))
  if (!is.null(adjustment.table)){
    df <- cbind(df,adjustment.table)
  }
  if(paired){
    if(n.g1 != n.g2){
      stop("Could not conduct paired diffVar analysis: unequal groupsizes")
    }
    df$xp <- as.factor(rep(1:n.g1,2))
  }
  formula.text <- paste0(c("~0",colnames(df)),collapse="+")
  design.m <- model.matrix(as.formula(formula.text),data=df)
  colnames(design.m) <- make.names(colnames(design.m),unique=TRUE)
  colnames(design.m)[1:2] <- c("group1","group2")
  var.fit <- tryCatch(varFit(data=meth.matrix,design=design.m),error=function(e){
    logger.warning("Could not compute p-values with diffVar, returning NA")
    return(NULL)
  })
  if(is.null(var.fit)){
    return(rep(NA,nrow(meth.matrix)))
  }
  contrasts.m <- makeContrasts(group1vs2=group1-group2,levels=design.m)
  var.fit <- contrasts.fit(var.fit,contrasts.m)
  var.fit <- eBayes(var.fit)
  logger.completed()
  return(var.fit$p.value[,"group1vs2"])
}

#' apply.iEVORA
#' 
#' This routine applies the iEVORA method created by Teschendorff et.al. to the supplied methylation matrix in a similar way
#' as the diffVar method.
#' 
#' @author Michael Scherer
#' @param meth.matrix Matrix containing the methylation information used to calculate differentially variable sites 
#'                    between the two groups
#' @param inds.g1 Indices in the phenotypic table corresponding to the first group.
#' @param inds.g2 Indices in the phenotypic table corresponding to the second group.
#' @return Q-values as the result of applying the iEVORA method and then correct for multiple testing.
#' @export
apply.iEVORA <- function(meth.matrix,inds.g1,inds.g2){
  if(is.logical(inds.g1)) inds.g1 <- which(inds.g1)
  if(is.logical(inds.g2)) inds.g2 <- which(inds.g2)
  n <- length(inds.g1)+length(inds.g2)
  pheno.binary <- rep(0,n)
  pheno.binary[inds.g2] <- 1
  result <- suppressWarnings(iEVORA(data.m=meth.matrix,pheno.v=pheno.binary))
  if(is.null(result)){
    p.values <- rep(1,dim(meth.matrix)[1])
    return(p.values)
  }
  p.values <- rep(1,dim(meth.matrix)[1])
  p.values[result[,"index"]] <- result[,"P(BT)"]
  return(p.values)
}

#' iEVORA
#' 
#' This method performs the iEVORA algorithm introduced by Teschendorff et.al. to detect sites showing differential
#' variability between two groups. The routine has been modified from the R script associated with the publication.
#' 
#' @author Andrew E. Teschendorff, with minor modifications by Michael Scherer
#' @param data.m Methylation matrix with rows labeling features, and columns labeling samples. 
#'               Rownames should be feature/probe IDs and should be provided.
#' @param pheno.v Phenotype vector with entries either 0 or 1 correpsonding to the group assignment.
#' @param thDV q-value threshold for the differntial variability test. Default is 0.001.
#' @param thDM p-value threshold for differntial methylation means. Default is 0.05.
#' 
#' @return A matrix of ranked differentially variable (DV) and differentially methylated CpGs (DVMCs), 
#'          ranked according to the t-statistic p-value, but selected using the Bartlett's DV test. 
#'          Columns labels are the t-statistic, its P-value, the mean of phenotype-1, the mean of phenotype-0, 
#'          the log-ratio of the variances of phenotype-1 to phenotype-0,as well as the Bartlett's test P-value 
#'          and q-value.
#' @noRd
#' @references Teschendorff, Andrew E, Gao, Yang, Jones, Allison, Ruebner, Matthias, Beckmann, Matthias W., Wachter, David L.
#'              ,Fasching, Peter A., Widschwendter, Martin
#'             DNA methylation outliers in normal breast tissue identify field defects that are enriched in cancer
#'             Nature Communications 7:10478.
iEVORA <- function(data.m,pheno.v,thDV=0.001,thDM=0.05){
  logger.start("iEVORA method")
  rnb.require("qvalue")
  rnb.require("stats")
  statDVC.m <- t(apply(data.m,1,doDV,pheno.v));
  qvDVC.v <- qvalue(statDVC.m[,2])$qval;
  dvc.idx <- which(qvDVC.v < thDV);
  nDVC <- length(dvc.idx);
  if( nDVC > 0 ){
    statDMC.m <- t(apply(data.m[dvc.idx,],1,doTT,pheno.v))
    tmp.s <- sort(statDMC.m[,2],decreasing=FALSE,index.return=TRUE);
    pvDMC.v <- tmp.s$x;
    ntop <- length(which(pvDMC.v < 0.05));
    if(ntop > 0){
      topDVMC.m <- cbind(dvc.idx[tmp.s$ix[1:ntop]],statDMC.m[tmp.s$ix[1:ntop],],statDVC.m[dvc.idx[tmp.s$ix[1:ntop]],c(3:4,1:2)],qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]]);
      colnames(topDVMC.m) <- c("index","t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)");
      rownames(topDVMC.m) <- rownames(statDMC.m)[tmp.s$ix[1:ntop]];
    }
    else {
      logger.info("No DVCs found.")
      return(NULL)
    }
  }
  else {
    logger.info("No DVCs detected. All p-values set to 1.");
    return(NULL)
  }
  logger.completed()
  return(topDVMC.m);
}

#' doDV
#' 
#' Helper function for the iEVORA method also supplied with the publication used to determine differentially variable
#' sites.
#' 
#' @author Andrew Teschendorff, minor modifications by Michael Scherer
#' @param tmp.v A row from the methylation matrix corresponding to a single site.
#' @param pheno.v Vector containing the splitting into the two groups for which differential variability should be detected.
#' 
#' @return Vector containing the log odds of variances, the p-value from the Bartlett's test and the mean methylation 
#'          levels for the two groups.
#' @noRd
doDV <- function(tmp.v,pheno.v){
  rnb.require('stats')
  co.idx <- which(pheno.v==0);
  ca.idx <- which(pheno.v==1);
  bt.o <- bartlett.test(x=tmp.v,g=pheno.v);
  pv <- bt.o$p.value;
  logR <- log2(var(tmp.v[ca.idx])/var(tmp.v[co.idx]));
  avCA <- mean(tmp.v[ca.idx]);
  avCO <- mean(tmp.v[co.idx]);
  out.v <- c(logR,pv,avCA,avCO);
  names(out.v) <- c("log(V1/V0)","P(BT)","Av1","Av0");
  return(out.v);
}

#' doTT
#' Helper function for the iEVORA method also supplied with the publication used to calculate t-statistics of differential
#' methylation between two sample groups.
#' 
#' 
#' @author Andrew Teschendorff, with minor modifications by Michael Scherer
#' @param tmp.v A row from the methylation matrix corresponding to a single site.
#' @param pheno.v Vector containing the splitting into the two groups for which differential variability should be detected.
#' 
#' @return Vector containing the t-statistics and the t-statistics p-value.
#' @noRd
doTT <- function(tmp.v,pheno.v){
  tt.o <- t.test(tmp.v ~ pheno.v);
  out.v <- c(-tt.o$stat,tt.o$p.val);
  names(out.v) <- c("t","P");
  return(out.v);
}

#' Returns the colum names of the differential variability table.
#' @param includeCovg Flag indicating if dataset contains coverage information
#' @return Column names of the differential variability table
get.table.ids <- function(includeCovg=FALSE){
  ret <- c("var.g1","var.g2","var.diff","var.log.ratio","diffVar.p.val")
  if(includeCovg){
    ret <- c(ret,"mean.covg.g1","covg.thresh.nsamples.g1","mean.covg.g2","covg.thresh.nsamples.g2")
  }
  return(ret)
}

#' addReportPlot.diffVar.scatter.site
#' 
#' This function creates a scatterplot for the given comparison comparing the difference in variances with the
#' usted) p-value as the result of the differentially variability analysis.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added
#' @param var.table Variability table containing the relevant information to be put into the scatterplot
#' @param comparison.name Comparison for which the scatterplot should be creaed
#' @param rank.cutoffs.numbers Numbers to be used for plotting the best ranked sites
#' @param auto.cutoff Auto selected cutoffs for the given comparison
#' @param group.name1 Name of the first group of the comparison
#' @param group.name2 Name of the second group of the comparison
#' @return list of scatterplots
#' @noRd
addReportPlot.diffVar.scatter.site <- function (report, var.table, comparison.name,
                              rank.cutoffs.numbers, auto.cutoff, rerank=TRUE,
                              group.name1="Group1",group.name2="Group2"){
  ret <- list()
  sparse.points <- DENS.SCATTER.SPARSE.POINTS.PERC
  if (DENS.SCATTER.SPARSE.POINTS.MAX < sparse.points*nrow(var.table)){
    sparse.points <- DENS.SCATTER.SPARSE.POINTS.MAX
  }
  #dens.subsample <- FALSE
  dens.subsample <- FALSE
  if(nrow(var.table)>10000){
    dens.subsample <- 10000
  }
  #if (nrow(var.table) > DENS.SCATTER.SUBSAMPLE.THRES){
  #  dens.subsample <- DENS.SCATTER.SUBSAMPLE.THRES
  #}
  
  if("diffVar.p.adj.fdr" %in% colnames(var.table)){
    var.sites <- var.table[,"diffVar.p.adj.fdr"] < P.VAL.CUT
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")],is.special=var.sites,
                                  dens.subsample=dens.subsample, add.text.cor=TRUE) +
        xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2)) + coord_fixed(xlim = c(0,max(var.table[,c("var.g1","var.g2")])),ylim = c(0,max(var.table[,c("var.g1","var.g2")])))
    comp.type <- "fdrAdjPval"
    fig.name <- paste("diffVar_site",comparison.name,comp.type,sep = "_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE, high.png = 200)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  ranks <- var.table[,"combinedRank.var"]
  if (rerank)	ranks <- rank(ranks,na.last="keep",ties.method="min")
  for(i in 1:length(rank.cutoffs.numbers)){
    number <- rank.cutoffs.numbers[i]
    cutoff.name <- paste0("rc",i)
    var.table$isDVC <- ranks < number
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")], is.special = var.table$isDVC,
                                  dens.subsample=dens.subsample, add.text.cor=TRUE) +
      xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2)) + coord_fixed(xlim = c(0,max(var.table[,c("var.g1","var.g2")])),ylim = c(0,max(var.table[,c("var.g1","var.g2")])))
    fig.name <- paste("diffVar_site",comparison.name,cutoff.name,sep="_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE, high.png = 200)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  if(is.integer(auto.cutoff)){
    var.table$isDVC <- ranks < auto.cutoff
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")], is.special = var.table$isDVC,
                                  dens.subsample=dens.subsample, add.text.cor=TRUE) +
      xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2))  + coord_fixed(xlim = c(0,max(var.table[,c("var.g1","var.g2")])),ylim = c(0,max(var.table[,c("var.g1","var.g2")])))
    fig.name <- paste("diffVar_site",comparison.name,"rcAuto",sep="_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE, high.png = 200)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  return(ret)
}
#' @noRd
addReportPlot.diffVar.scatter.region <- function (report, var.table, comparison.name, region.name, rerank=TRUE,
                                           ranking.cutoffs, auto.cutoff=NULL, group.name1="Group1",group.name2="Group2",useSiteCols=FALSE){
  ret <- list()
  
  cn.x <- "mean.var.g1"
  cn.y <- "mean.var.g2"
  cn.pa <- "comb.p.adj.var.fdr"
  if (useSiteCols){
    cn.x <- "var.g1"
    cn.y <- "var.g2"
    cn.pa <- "diffVar.p.adj.fdr"
  }
  al.x <- paste("Mean Variance",group.name1,sep=".")
  al.y <- paste("Mean Variance",group.name2,sep=".")
  
  dens.subsample <- FALSE
  if(nrow(var.table) > DENS.SCATTER.SUBSAMPLE.THRES){
    dens.subsample <- DENS.SCATTER.SUBSAMPLE.THRES
  }
  
  if(cn.pa %in% colnames(var.table)){
    var.sites <- var.table[,cn.pa] < P.VAL.CUT
    plot <- create.densityScatter(var.table[,c(cn.x,cn.y)],is.special=var.sites,
                                  dens.subsample = dens.subsample, add.text.cor=TRUE) +
      xlab(al.x) + ylab(al.y) + coord_fixed(xlim = c(0,max(var.table[,c("mean.var.g1","mean.var.g2")])),ylim = c(0,max(var.table[,c("mean.var.g1","mean.var.g2")])))
    comp.type <- "fdrAdjPval"
    fig.name <- paste("diffVar_region",comparison.name,region.name,comp.type,sep = "_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  
  rrs <- var.table[,"combinedRank.var"]
  if (rerank)	rrs <- rank(rrs,na.last="keep",ties.method="min")
  for (i in 1:length(ranking.cutoffs)){
    rc <- ranking.cutoffs[i]
    cur.cut.name <- paste("rc",i,sep="")
    var.table$isDVR <- rrs < rc
    
    pp <- create.densityScatter(var.table[,c(cn.x,cn.y)],is.special=var.table$isDVR,
                                dens.subsample = dens.subsample, add.text.cor=TRUE) +
      labs(x=al.x, y=al.y) + coord_fixed(xlim = c(0,max(var.table[,c("mean.var.g1","mean.var.g2")])),ylim = c(0,max(var.table[,c("mean.var.g1","mean.var.g2")])))
    
    figName <- paste("diffVar_region",comparison.name,region.name,cur.cut.name,sep="_")
    report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,list(report.plot))
  }
  
  if (is.integer(auto.cutoff)){
    var.table$isDVR <- var.table[,"combinedRank.var"] <= auto.cutoff
    pp <- create.densityScatter(var.table[,c(cn.x,cn.y)],is.special=var.table$isDVR,
                                dens.subsample = dens.subsample, add.text.cor=TRUE) +
      labs(x=al.x, y=al.y) + coord_fixed(xlim = c(0,max(var.table[,c("mean.var.g1","mean.var.g2")])),ylim = c(0,max(var.table[,c("mean.var.g1","mean.var.g2")])))
    figName <- paste("diffVar_region",comparison.name,region.name,"rcAuto",sep="_")
    report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,list(report.plot))
  }

  return(ret)
}

#' addReportPlot.diffVar.volcano
#' 
#' This function creates volcano plots for the given comparison comparing the difference in variances with the
#' negative decadic logarithm of the (adjusted) p-value as the result of the differentially variability analysis.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added
#' @param var.table Variability table containing the relevant information to be put into the scatterplot
#' @param comparison.name Comparison for which the scatterplot should be creaed
#' @param group.name1 Name of the first group of the comparison
#' @param group.name2 Name of the second group of the comparison
#' @aliases addReportPlot.diffVar.volcano
#' @aliases addReportPlot.diffVar.volcano.region
#' @return list of volcano plots
#' @noRd
addReportPlot.diffVar.volcano <- function (report, var.table, comparison.name,
                                           group.name1="Group1",group.name2="Group2"){
  ret <- list()
  dont.plot.p.val <- all(is.na(var.table[,"diffVar.p.val"]))
  
  fig.name <- paste("diffVar_volcano_site",comparison.name,"pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes(x=var.diff,y=log10P,color=log10(combinedRank.var)) + 
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point()
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  fig.name <- paste("diffVar_volcano_site",comparison.name,"pValAdj",sep="_")
  pp <- ggplot(var.table) + aes(x=var.diff,y=log10FDR,color=log10(combinedRank.var)) +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient")))+
    geom_point()
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  return(ret)
}
#' addReportPlot.diffVar.volcano.region
#' @noRd
addReportPlot.diffVar.volcano.region <- function(report, var.table,comparison.name,region.type,
                                                  group.name1,group.name2, useSiteCols=FALSE){
  cn.d <- "mean.var.diff"
  cn.q <- "mean.var.log.ratio"
  cn.p <- "comb.p.val.var"
  cn.pa <- "comb.p.adj.var.fdr"
  if (useSiteCols){
    cn.d <- "var.diff"
    cn.q <- "var.log.ratio"
    cn.p <- "diffVar.p.val"
    cn.pa <- "diffVar.p.adj.fdr"
  }

  figPlots <- list()
  dont.plot.p.val <- all(is.na(var.table[,cn.p]))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"diff","pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes_string(cn.d, paste0("-log10(",cn.p,")"), color="log10(combinedRank.var)") +
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point()#(alpha=0.3)
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"diff","pValAdj",sep="_")
  pp <- ggplot(var.table) + aes_string(cn.d, paste0("-log10(",cn.pa,")"), color="log10(combinedRank.var)") +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
    geom_point()#(alpha=0.3)
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"quot","pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes_string(cn.q, paste0("-log10(",cn.p,")"), color="log10(combinedRank.var)") +
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point()#(alpha=0.3)
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"quot","pValAdj",sep="_")
  pp <- ggplot(var.table) + aes_string(cn.q, paste0("-log10(",cn.pa,")"), color="log10(combinedRank.var)") +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
    geom_point()#(alpha=0.3)
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"diff","quotSig",sep="_")
  pp <- ggplot(var.table) + aes_string(cn.d, cn.q, color="log10(combinedRank.var)") +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
    geom_point()#(alpha=0.3)
  report.plot <- createReportGgPlot(pp, figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"quot","quotSig",sep="_")
  pp <- rnb.message.plot("Quotient--Quotient scatterplot not available")
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  return(figPlots)
}

#' addReportPlot.diffVar.meth
#' 
#' This function compares differentially methylated sites (DMCs) with differentially variable sites (DVCs) in a scatterplot
#' and colors significantly differential sites with different colors.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added.
#' @param var.table Table containing the differentially variable sites.
#' @param comparison.name Name of the comparison to be conducted.
#' @param group.name1 Name of the first group in the comparison
#' @param group.name2 Name of the second group in the comparison
#' @param auto.diffVar Automatically selected rank cutoff for differential varibility
#' @param auto.diffMeth Automatically selected rank cutoff for differential methylation
#' @param rank.cuts Cutoffs used for differentiality detection
#' @return The modified report
#' @noRd
addReportPlot.diffVar.meth <- function(report, var.table, comparison.name,
                           group.name1,group.name2,auto.diffVar=NULL,auto.diffMeth=NULL,rank.cuts=c(100,1000,10000,100000),rerank=TRUE){
  ret <- list()
  dens.subsample <- FALSE
  if (nrow(var.table) > 10000){
    dens.subsample <- 10000
  }

  ranks.diffVar <- var.table[,"combinedRank.var"]
  if(rerank) ranks.diffVar <- rank(ranks.diffVar,na.last = "keep",ties.method = "min")
  ranks.diffMeth <- var.table[,"combinedRank"]
  if(rerank) ranks.diffMeth <- rank(ranks.diffMeth,na.last = "keep",ties.method = "min")
  
  fig.name <- paste("diffVarMeth_site",comparison.name,sep="_")
  for(i in 1:length(rank.cuts)){
    cutoff <- rank.cuts[i]
    cut.id <- paste0("rc",i)
    toPlot <- data.frame(Rank.meth=var.table[,"combinedRank"],Rank.var=var.table[,"combinedRank.var"])
    is.diff.meth <- ranks.diffMeth < cutoff
    is.diff.var <- ranks.diffVar < cutoff
    rank.cut.diffMeth <- max(toPlot[is.diff.meth,1])
    rank.cut.diffVar <- max(toPlot[is.diff.var,2])
    color.diff <- rep('Not Significant',dim(toPlot)[1])
    color.diff[is.diff.meth] <- "DMC"
    color.diff[is.diff.var] <- "DVC"
    color.diff[is.diff.meth&is.diff.var] <- "Both"
    toPlot <- data.frame(toPlot,Color.diff=color.diff)
    is.special <- is.diff.meth|is.diff.var
    pp <- create.diffMeth.diffVar.subsample(toPlot,dens.subsample, is.special=is.special, rank.cut.diffMeth = rank.cut.diffMeth, rank.cut.diffVar = rank.cut.diffVar)
    fig.name <- paste("diffMethVar",comparison.name,cut.id,sep="_")
    report.plot <- createReportGgPlot(pp,fig.name,report,create.pdf = FALSE,high.png = 200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,report.plot)
  }
  if(is.integer(auto.diffVar)&&is.integer(auto.diffMeth)){
    toPlot <- data.frame(Rank.meth=var.table[,"combinedRank"],Rank.var=var.table[,"combinedRank.var"])
    is.diff.meth <- toPlot$Rank.meth < auto.diffMeth
    is.diff.var <- toPlot$Rank.var < auto.diffVar
    color.diff <- rep('Not Significant',dim(toPlot)[1])
    color.diff[is.diff.meth] <- "DMC"
    color.diff[is.diff.var] <- "DVC"
    color.diff[is.diff.meth&is.diff.var] <- "Both"
    toPlot <- data.frame(toPlot,Color.diff=color.diff)
    is.special <- is.diff.meth|is.diff.var
    pp <- create.diffMeth.diffVar.subsample(toPlot,dens.subsample,is.special=is.special, rank.cut.diffMeth = auto.diffMeth, rank.cut.diffVar = auto.diffVar)
    fig.name <- paste("diffMethVar",comparison.name,"rcAuto",sep="_")
    report.plot <- createReportGgPlot(pp,fig.name,report,create.pdf = FALSE,high.png = 200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,report.plot)
  }
  return(ret)
}

#' addReportPlot.diffVar.meth.region
#' 
#' This function compares differentially methylated sites (DMCs) with differentially variable sites (DVCs) in a scatterplot
#' and colors significantly differential sites with different colors.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added.
#' @param var.table Table containing the differentially variable sites.
#' @param comparison.name Name of the comparison to be conducted.
#' @param region.type Region type used
#' @param group.name1 Name of the first group in the comparison
#' @param group.name2 Name of the second group in the comparison
#' @param auto.diffVar Automatically selected rank cutoff for differential varibility
#' @param auto.diffMeth Automatically selected rank cutoff for differential methylation
#' @param rank.cuts Cutoffs used for differentiality detection
#' @return The modified report
#' @noRd
addReportPlot.diffVar.meth.region <- function(report, var.table, comparison.name, region.type,
                                       group.name1,group.name2,auto.diffVar=NULL,auto.diffMeth=NULL,rank.cuts=c(100,500,1000),rerank=TRUE){
  ret <- list()
  dens.subsample <- FALSE
  if (nrow(var.table) > 10000){
    dens.subsample <- 10000
  }
  
  ranks.diffVar <- var.table[,"combinedRank.var"]
  if(rerank) ranks.diffVar <- rank(ranks.diffVar,na.last = "keep",ties.method = "min")
  ranks.diffMeth <- var.table[,"combinedRank"]
  if(rerank) ranks.diffMeth <- rank(ranks.diffMeth,na.last = "keep",ties.method = "min")
  
  fig.name <- paste("diffVarMeth_region",comparison.name,region.type,sep="_")
  for(i in 1:length(rank.cuts)){
    cutoff <- rank.cuts[i]
    cut.id <- paste0("rc",i)
    toPlot <- data.frame(Rank.meth=var.table[,"combinedRank"],Rank.var=var.table[,"combinedRank.var"])
    is.diff.meth <- ranks.diffMeth < cutoff
    is.diff.var <- ranks.diffVar < cutoff
    rank.cut.diffMeth <- max(toPlot[is.diff.meth,1])
    rank.cut.diffVar <- max(toPlot[is.diff.var,2])
    color.diff <- rep('Not Significant',dim(toPlot)[1])
    color.diff[is.diff.meth] <- "DMR"
    color.diff[is.diff.var] <- "DVR"
    color.diff[is.diff.meth&is.diff.var] <- "Both"
    toPlot <- data.frame(toPlot,Color.diff=color.diff)
    is.special <- is.diff.meth|is.diff.var
    pp <- create.diffMeth.diffVar.subsample(toPlot,dens.subsample, is.special=is.special, rank.cut.diffMeth = rank.cut.diffMeth, rank.cut.diffVar = rank.cut.diffVar)
    fig.name <- paste("diffMethVar",comparison.name,region.type,cut.id,sep="_")
    report.plot <- createReportGgPlot(pp,fig.name,report,create.pdf = FALSE,high.png = 200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,report.plot)
  }
  if(is.integer(auto.diffVar)&&is.integer(auto.diffMeth)){
    toPlot <- data.frame(Rank.meth=var.table[,"combinedRank"],Rank.var=var.table[,"combinedRank.var"])
    is.diff.meth <- toPlot$Rank.meth < auto.diffMeth
    is.diff.var <- toPlot$Rank.var < auto.diffVar
    color.diff <- rep('Not Significant',dim(toPlot)[1])
    color.diff[is.diff.meth] <- "DMR"
    color.diff[is.diff.var] <- "DVR"
    color.diff[is.diff.meth&is.diff.var] <- "Both"
    toPlot <- data.frame(toPlot,Color.diff=color.diff)
    is.special <- is.diff.meth|is.diff.var
    pp <- create.diffMeth.diffVar.subsample(toPlot,dens.subsample, is.special=is.special, rank.cut.diffMeth = auto.diffMeth, rank.cut.diffVar = auto.diffVar)
    fig.name <- paste("diffMethVar",comparison.name,region.type,"rcAuto",sep="_")
    report.plot <- createReportGgPlot(pp,fig.name,report,create.pdf = FALSE,high.png = 200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,report.plot)
  }
  return(ret)
}

#' create.diffMeth.diffVar.subsample
#' 
#' This routine creates a plot showing both differentially methylated and differentially variable sites for the given rank cutoffs,
#' where the x-axis is the combined rank of the differential methylation and the y-axis the combined rank of the differential
#' variability analysis.
#' @param df2p Data frame containing the ranks and color of the sites
#' @param dens.subsample Number specifiying the subsample that should be taken from the whole population of sites
#' @param rank.cut.diffMeth Rank cutoff for differential methylation (absolute rank)
#' @param rank.cut.diffVar Rank cutoff for differential variability (absolute rank)
#' @return Plot object with the corresponding plot
#' @noRd
create.diffMeth.diffVar.subsample <- function(df2p,dens.subsample,is.special=NULL,rank.cut.diffMeth, rank.cut.diffVar,sparse.points=0.01,dens.n=100){
  if (!(is.numeric(sparse.points) && sparse.points>=0)) {
    stop("Invalid parameter value: sparse.points")
  }
  if (is.null(df2p) || nrow(df2p)<1){
    logger.warning(c("Could not create density scatterplot"))
    pp <- rnb.message.plot("Could not create plot")
    return(pp)
  }
  df2p <- na.omit(df2p)
  if (is.null(df2p) || nrow(df2p)<1){
    logger.warning(c("Could not create density scatterplot (NA omission removed all entries)"))
    pp <- rnb.message.plot("Could not create plot")
    return(pp)
  }
  df2p.sub <- df2p
  dens.ranks <- NULL
  tryCatch(
    dens.ranks <- densRanks(x=df2p[,1],y=df2p[,2]),
    error=function(ee){
      logger.warning(c("Could not assess density ranking:",ee$message))
    }
  )
  if (is.numeric(dens.subsample) && dens.subsample>0){
    ss <- as.integer(dens.subsample)
    if (nrow(df2p) > ss) {
      df2p.sub <- df2p[sample(nrow(df2p),ss),]
    }
  }
  
  #the standard bandwith function of MASS::kde2d is unstable when looking at
  #distributions with very low variance. Here's a more stable version
  stable.bandwidth.fun <- function(x,eps=1e-4){
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2] - r[1])/1.34
    if (h==0) h <- eps
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
  }
  stable.h <- c(stable.bandwidth.fun(df2p.sub[,1]),stable.bandwidth.fun(df2p.sub[,2]))
  
  if (is.null(dens.ranks)){
    pp <- rnb.message.plot("Could not assess density")
  } else {
    pp <- ggplot(df2p.sub) + aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2]) + 
      stat_density2d(geom="tile", fill=DENS.COLORS.LOW[1], aes(alpha=after_stat(density)^0.25), contour=FALSE, n=dens.n, h=stable.h) +
      scale_alpha(range = c(0.0, 1),guide=FALSE) +
      geom_vline(xintercept = rank.cut.diffMeth,linetype='dotted') + geom_hline(yintercept = rank.cut.diffVar,linetype='dotted')+
      scale_color_manual(values=c(rnb.getOption("colors.category")[c(1,2,4)],DENS.COLORS.LOW[1]))
    if (sparse.points > 0){
      if (sparse.points <= 1){
        thres <- ceiling(nrow(df2p)*sparse.points)
      } else {
        thres <- sparse.points
      }
      df2p.loose <- df2p[dens.ranks<=thres,]#the sub data.frame in of the least dens points
      pp <- pp + geom_point(data=df2p.loose,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2],colour=colnames(df2p)[3]),size=0.4)
    }
    if(!is.null(is.special)){
      df2p.special <- df2p[is.special,]
      if(is.numeric(dens.subsample) && dens.subsample>0){
        ss <- as.numeric(dens.subsample)
        if(nrow(df2p.special)>ss){
          pp <- pp + geom_point(data=df2p.special,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2],colour=colnames(df2p)[3]),size=0.4) + guides(fill=FALSE)
        }else{
          pp <- pp + geom_point(data=df2p.special,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2],colour=colnames(df2p)[3]),size=1) + guides(fill=FALSE)
        }
      }else{
        pp <- pp + geom_point(data=df2p.special,aes_string(x=colnames(df2p)[1],y=colnames(df2p)[2],colour=colnames(df2p)[3]),size=1) + guides(fill=FALSE)
      }
    }  
  }
  return(pp)
}

#' rnb.section.diffVar
#'
#' This function adds information about the differential variability analysis to the specified report.
#' 
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}} on which the analysis was conducted.
#' @param diff.meth Object of class \code{\linkS4class{RnBDiffMeth}} as the result of applying \code{rnb.execute.diffVar} 
#'                  containing the differentially variable sites
#' @param report Report object (\code{\linkS4class{Report}}) to which information should be added.
#' @param gzTable Flag indicating if the tables should ne gzipped
#' @param level Which level of section should be created. See \code{rnb.add.section}.
#' @return Modified report with section about result of differential variability analysis
#' @author Michael Scherer
#' @noRd
rnb.section.diffVar <- function(rnb.set,diff.meth,report,gzTable=FALSE,differentiality.method=rnb.getOption("differential.variability.method"),level=1){
  if (length(get.comparisons(diff.meth))<1){
    stop("no valid comparisons")
  }
  diff.var.cutoff <- 0.01
  rank.cutoffs.numbers <- c(100,1000,10000,100000)
  logger.start("Adding Differential Variability Information")
  txt <- c("Differentially variable sites were computed with <code>", differentiality.method, "</code>.",
           " For more information about the method, have a look at the <a href=")
  if(differentiality.method=='diffVar'){
    refText.missMethyl <- c("Phipson, B., & Oshlack, A. (2014). DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging. ",
                      "<i>Genome Biology</i>, <b>15</b>(9), 465")
    report <- rnb.add.reference(report, refText.missMethyl)
    txt <- c(txt,"https://bioconductor.org/packages/release/bioc/html/missMethyl.html>missMethyl</a> Bioconductor package.",
             rnb.get.reference(report,refText.missMethyl))
  }else if (differentiality.method=='iEVORA'){
    refText.iEVORA <- c("Teschendorff, A.E. et.al., DNA methylation outliers in normal breast tissue identify field defects that are enriched in cancer",
                        "<i>Nature Communications</i>, <b>7</b>, 10478")
    txt <- c(txt,"https://www.nature.com/articles/ncomms10478>iEVORA</a> method published in Nature.",
             "For the iEVORA method, adjustment columns are not supported and therefore ignored.",
             rnb.get.reference(report,refText.iEVORA))
  }
  txt <- c(txt," This section contains plots and tables describing the results of this test and further analyses ",
            "of the sites that were selected as differentially variable. Please note that missing methylation values have ",
           "been imputed with ",rnb.getOption("imputation.method"),".")
  report <- rnb.add.section(report, 'Differential Variability', txt,level = level)
  
  diffSitesRankCut <- c(100,1000,10000,100000)
  comps <- get.comparisons(diff.meth)

  logger.start("Selection of rank cutoffs")
  rank.cuts.auto <- lapply(1:length(comps),FUN=function(i){
    comp.name <- names(comps)[i]
    comp <- comps[[i]]
    var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
    res <- auto.select.rank.cut(var.table$diffVar.p.adj.fdr,var.table$combinedRank.var,alpha=0.1)
    return(as.integer(res))
  })
  rank.cuts.auto.meth <- lapply(1:length(comps),FUN=function(i){
    comp.name <- names(comps)[i]
    comp <- comps[[i]]
    var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
    res <- auto.select.rank.cut(var.table$diffmeth.p.adj.fdr,var.table$combinedRank,alpha=0.1)
    return(as.integer(res))
  })
  txt <- paste("The following rank cutoffs have been automatically selected for the analysis of differentially",
               "variable sites:")
  rnb.add.paragraph(report, txt)
  
  tt <- data.frame(matrix(unlist(rank.cuts.auto),ncol=1,nrow=length(comps),byrow=TRUE))
  colnames(tt) <- "Rank Cutoff"
  rownames(tt) <- comps
  rnb.add.table(report,tt)
  logger.completed()
  
  #scatterplots
  logger.start("Adding scatterplots")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      selected.cutoff <- rank.cuts.auto[[i]]
      res <- addReportPlot.diffVar.scatter.site(report,var.table,comp.name.allowed,
                                                      rank.cutoffs.numbers=rank.cutoffs.numbers,
                                                      auto.cutoff=selected.cutoff, rerank = TRUE,
                                                      group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      selected.cutoff <- rank.cuts.auto[[i]]
      res <- addReportPlot.diffVar.scatter.site(report,var.table,comp.name.allowed,
                                                rank.cutoffs.numbers=rank.cutoffs.numbers,
                                                auto.cutoff=selected.cutoff,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  
  comps <- get.comparisons(diff.meth)
  diff.var.type <- c(paste("FDR adjusted p-value &lt;",P.VAL.CUT),
                     paste("combined rank among the",rank.cutoffs.numbers,"best ranking sites"),
                     "automatically selected rank cutoff")
  names(diff.var.type) <-  c("fdrAdjPval",paste0("rc",1:length(rank.cutoffs.numbers)),"rcAuto")
  setting.names <- list(
    'comparison' = comps,
    'differential variability measure' = diff.var.type)
  description <- c('Scatterplot for differential variable sites.',
                   ' The transparency corresponds to point density. If the number of points exceeds ',DENS.SCATTER.SUBSAMPLE.THRES,
                   ' then the number of points for density estimation is reduced to that number by random sampling.',
                   ' The',round(DENS.SCATTER.SPARSE.POINTS.PERC*100),
                   '% of the points in the sparsest populated plot regions are drawn explicitly (up to a maximum of ',DENS.SCATTER.SPARSE.POINTS.MAX,
                   " points).",
                   ' Additionally, the colored points represent differentially variable sites (according to the selected criterion).')
  report <- rnb.add.figure(report, description, plots, setting.names)
  logger.completed()
  
  #volcano plots
  logger.start("Adding volcano plots")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.volcano(report, var.table, comp.name.allowed,
                                                      group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- comp.name
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.volcano(report, var.table, comp.name.allowed,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  
  signif.measure <- c("pVal"="p-value","pValAdj"="adjusted p-value")
  setting.names <- list(
    'comparison' = comps,
    'significance metric' = signif.measure)
  description <- 'Volcano plot for differential variable sites.'
  report <- rnb.add.figure(report, description, plots, setting.names)
  logger.completed()
  
  #Comparing differentially methylated sites and differentially variable
  logger.start("Comparing DMCs and DVCs")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      auto.cutoff <- rank.cuts.auto[[i]]
      auto.cutoff.meth <- rank.cuts.auto.meth[[i]]
      res <- addReportPlot.diffVar.meth(report, var.table, comp.name.allowed,auto.diffVar = auto.cutoff, auto.diffMeth = auto.cutoff.meth,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      auto.cutoff <- rank.cuts.auto[[i]]
      auto.cutoff.meth <- rank.cuts.auto.meth[[i]]
      res <- addReportPlot.diffVar.meth(report, var.table, comp.name.allowed,auto.diffVar = auto.cutoff, auto.diffMeth = auto.cutoff.meth,
                                        group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  diffVarType = c(paste("combined rank among the ",diffSitesRankCut," best ranking sites",sep=""),
                  "automatically selected rank cutoff")
  names(diffVarType) = c(paste("rc",1:length(diffSitesRankCut),sep=""),"rcAuto")
  setting.names <- list(
    'comparison' = comps,
    'rankCutoff' = diffVarType)
  description <- c('Scatterplot comparing differentially methylated (DMCs) and variable sites (DVCs), as well as sites ',
  'that are both differentially methylated and variable. The dotted lines corrspond to the respective rank cutoffs used ',
  'to call a site differentially methylated/variable.')
  report <- rnb.add.figure(report, description, plots, setting.names)
  logger.completed()
  
  logger.completed()
  return(report)
}

#' rnb.section.diffVar.region
#'
#' Adds information for differentially variable regions to the report.
#' @author Michael Scherer
#' @param rnb.set Object of type \code{\linkS4class{RnBSet}} containing methylation information
#' @param diff.meth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param report Report object to which the content is added
#' @param gzTable Flag indicating if tables should be gzipped
#' @param level Which level of section should be created. See \code{rnb.add.section}.
#' @return The modified report object
rnb.section.diffVar.region <- function(rnb.set,diff.meth,report,gzTable=FALSE,level=1){
  if (length(get.comparisons(diff.meth))<1){
    stop("no valid comparisons")
  }
  if (length(get.region.types(diff.meth))<1){
    stop("no valid region types")
  }

  skipSites <- !includes.sites(diff.meth)
  diffRegionRankCut <- c(100,500,1000)
  logger.start("Adding Region Level Information (Differential Variability)")

  sectionText <- c("Differential variability on the region level was computed similar to differential methylation, ", 
                   "but the mean of variances, ",
                   "the log-ratio of the quotient of variances as well as the p-values from the differentiality test were ",
                   "employed. Ranking was performed in line with the ranking of differential methylation.")
  report <- rnb.add.section(report, "Differential Variability", sectionText,level = level)
  
  comps <- get.comparisons(diff.meth)
  reg.types <- get.region.types(diff.meth)
  
  logger.start("Selection of rank cutoffs")
  rank.cuts.auto <- lapply(1:length(comps),FUN=function(i){
    lapply(1:length(reg.types),FUN=function(j){
      var.table <- get.table(diff.meth,comps[i],region.type=reg.types[j],return.data.frame=TRUE)
      res <- auto.select.rank.cut(var.table$comb.p.adj.var.fdr,var.table$combinedRank.var,alpha=0.1)
      return(as.integer(res))
    })
  })
  rank.cuts.auto.meth <- lapply(1:length(comps),FUN=function(i){
    lapply(1:length(reg.types),FUN=function(j){
      var.table <- get.table(diff.meth,comps[i],region.type=reg.types[j],return.data.frame=TRUE)
      res <- auto.select.rank.cut(var.table$comb.p.adj.fdr,var.table$combinedRank,alpha=0.1)
      return(as.integer(res))
    })
  })
  txt <- paste("The following rank cutoffs have been automatically selected for the analysis of differentially",
               "variable regions:")
  rnb.add.paragraph(report, txt)
  
  tt <- data.frame(matrix(unlist(rank.cuts.auto),ncol=length(reg.types),nrow=length(comps),byrow=TRUE))
  colnames(tt) <- reg.types
  rownames(tt) <- comps
  rnb.add.table(report,tt)
  logger.completed()
  
  #scatterplots
  logger.start("Adding scatterplots")
  rnb.cleanMem()
  addedPlots <- list()
  grp.labels <- get.comparison.grouplabels(diff.meth)
  if(parallel.isEnabled()){
    iis <- 1:length(comps)
    jjs <- 1:length(reg.types)
    pps <- expand.grid(iis,jjs)
    
    addedPlots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
      i <- pps[k,1]
      j <- pps[k,2]
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      rr <- reg.types[j]
      rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      auto.rank.cut <- rank.cuts.auto[[i]][[j]]
      var.table <- get.table(diff.meth,ccc,rr,return.data.frame=TRUE)
      res <- addReportPlot.diffVar.scatter.region(report,var.table,comparison.name=ccn,region.name=rrn,
                                                  ranking.cutoffs = diffRegionRankCut, rerank=TRUE,
                                                  auto.cutoff = auto.rank.cut,group.name1 = grp.labels[ccc,1],group.name2 = grp.labels[ccc,2], useSiteCols=skipSites)
      rnb.cleanMem()
      res
    }
    
  } else {
    for (i in 1:length(comps)){
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      for (j in 1:length(reg.types)){
        rr <- reg.types[j]
        rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
        auto.rank.cut <- rank.cuts.auto[[i]][[j]]
        var.table <- get.table(diff.meth,ccc,rr,return.data.frame=TRUE)
        addedPlots <- c(addedPlots,addReportPlot.diffVar.scatter.region(report,var.table,comparison.name=ccn,region.name=rrn,
                                                    ranking.cutoffs = diffRegionRankCut,
                                                    auto.cutoff = auto.rank.cut,group.name1 = grp.labels[ccc,1],group.name2 = grp.labels[ccc,2], useSiteCols=skipSites)
        )
       rnb.cleanMem()
      }
    }
  }
  
  names(reg.types) <- ifelse(is.valid.fname(reg.types),reg.types,paste("reg",1:length(reg.types),sep=""))
  diffVarType = c(paste("FDR adjusted p-value <",P.VAL.CUT),
                   paste("combined rank among the ",diffRegionRankCut," best ranking regions",sep=""),
                   "automatically selected rank cutoff")
  names(diffVarType) = c("fdrAdjPval",paste("rc",1:length(diffRegionRankCut),sep=""),"rcAuto")
  setting.names <- list(
    'comparison' = comps,
    'regions' = reg.types ,
    'differential variability measure' = diffVarType)
  description <- 'Scatterplot for differential variable regions. The transparency corresponds to point density. The 1% of the points in the sparsest populated plot regions are drawn explicitly.
  Additionally, the colored points represent differentially methylated regions (according to the selected criterion).'
  report <- rnb.add.figure(report, description, addedPlots, setting.names)
  logger.completed()
  
  #volcano plots
  logger.start("Adding volcano plots")
  rnb.cleanMem()
  addedPlots <- list()
  if(parallel.isEnabled()){
    #generate pairs of comparison, region combinations with indices
    iis <- 1:length(comps)
    jjs <- 1:length(reg.types)
    pps <- expand.grid(iis,jjs)

    addedPlots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
      i <- pps[k,1]
      j <- pps[k,2]
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      rr <- reg.types[j]
      rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      var.table <- get.table(diff.meth,ccc,region.type=rr,return.data.frame=TRUE)
      res <- addReportPlot.diffVar.volcano.region(report,var.table,comparison.name=ccn,region.type=rrn,
                                                        group.name1=grp.labels[ccc,1],group.name2=grp.labels[ccc,2], useSiteCols=skipSites)
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(comps)){
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      for (j in 1:length(reg.types)){
        rr <- reg.types[j]
        rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      var.table <- get.table(diff.meth,ccc,region.type=rr,return.data.frame=TRUE)
      addedPlots <- c(addedPlots,addReportPlot.diffVar.volcano.region(report,var.table,comparison.name=ccn,region.type=rrn,
                                                        group.name1=grp.labels[ccc,1],group.name2=grp.labels[ccc,2], useSiteCols=skipSites)
      )
      rnb.cleanMem()
      }
    }
  }
  diff.measure <- c("diff"="Difference","quot"="Quotient")
  signif.measure <- c("pVal"="combined p-value","pValAdj"="adjusted combined p-value","quotSig"="Quotient (only meaningful if 'Difference' is selected above)")
  setting.names <- list(
    'comparison' = comps,
    'regions' = reg.types ,
    'difference metric' = diff.measure,
    'significance metric' = signif.measure)
  description <- 'Volcano plot for differential variability quantified by various metrics. Color scale according to
  combined ranking.'
  report <- rnb.add.figure(report, description, addedPlots, setting.names)
  logger.completed()
  
  #Comparing differentially methylated regions and differentially variable
  logger.start("Comparing DMRs and DVRs")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
      i <- pps[k,1]
      j <- pps[k,2]
      comp.name <- names(comps)[i]
      comp <- comps[comp.name]
      ccn <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      rr <- reg.types[j]
      rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      var.table <- get.table(diff.meth,comp,region.type=rr,return.data.frame=TRUE)
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      auto.cutoff <- rank.cuts.auto[[i]][[j]]
      auto.cutoff.meth <- rank.cuts.auto.meth[[i]][[j]]
      res <- addReportPlot.diffVar.meth.region(report, var.table, ccn, region.type=rrn, auto.diffVar = auto.cutoff, auto.diffMeth = auto.cutoff.meth,
                                        group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(comps)){
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      for (j in 1:length(reg.types)){
        rr <- reg.types[j]
        rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
        var.table <- get.table(diff.meth,ccc,region.type=rr,return.data.frame=TRUE)
        grp.names <- get.comparison.grouplabels(diff.meth)[ccc,]
        auto.cutoff <- rank.cuts.auto[[i]][[j]]
        auto.cutoff.meth <- rank.cuts.auto.meth[[i]][[j]]
        plots <- c(plots,addReportPlot.diffVar.meth.region(report,var.table,ccn,region.type=rrn,auto.diffVar = auto.cutoff, auto.diffMeth = auto.cutoff.meth,
                                                                        group.name1=grp.labels[ccc,1],group.name2=grp.labels[ccc,2])
        )
        rnb.cleanMem()
      }
    }
  }
  diffVarType = c(paste("combined rank among the ",diffRegionRankCut," best ranking regions",sep=""),
                  "automatically selected rank cutoff")
  names(diffVarType) = c(paste("rc",1:length(diffRegionRankCut),sep=""),"rcAuto")
  setting.names <- list(
    'comparison' = comps,
    'regions' = reg.types,
    'rankCutoff' = diffVarType)
  description <- c('Scatterplot comparing differentially methylated (DMRs) and variable regions (DVRs), as well as regions ',
  'that are both differentially methylated and variable. The dotted lines corrspond to the respective rank cutoffs used ',
  'to call a region differentially methylated/variable.')
  report <- rnb.add.figure(report, description, plots, setting.names)
  

  logger.completed()
  return(report)
}

#' cols.to.rank.site
#' 
#' Return a matrix containing the negative absolute values of the information used to rank the sites. Those are currently:
#' the variance difference, the log ratio in variances and the p-value from the statistical test.
#' 
#' @param diff.var A differential variability table.
#' @return A matrix with the absolute values of the relevant columns
#' @author Michael Scherer
#' @rdname cols.to.rank
#' @aliases cols.to.rank.site
#' @aliases cols.to.rank.region
cols.to.rank.site <- function(diff.var){
  return(cbind(-abs(diff.var$var.diff), -abs(diff.var$var.log.ratio), diff.var$diffVar.p.val))
}
#' @rdname cols.to.rank
cols.to.rank.region <- function(diff.var){
  return(cbind(-abs(diff.var$mean.var.diff), -abs(diff.var$mean.var.log.ratio), diff.var$comb.p.adj.var.fdr))
}

#' computeDiffVar.bin.site
#' 
#' This function computes the table of differentially variable sites with the corresponding statistics.
#' 
#' @param X methylation matrix with sites as rows and samples as columns
#' @param inds.g1 the column indicies of the first group
#' @param inds.g2 the column indicies of the second group
#' @param variability.method method to use for calculating p-values. One of `diffVar` or `iEVORA`
#' @param paired Should paired analysis be conducted? If yes, the first entry in \code{inds.g1} should correspond to the first
#'        in \code{indg.g2} and so on.
#' @param eps Epsilon for computing quotients (avoid division by 0 by adding this value to denominator and enumerator before calculating the quotient)
#' @param adjustment.table Table specifying the variables to be adjusted for in the analysis.
#' @noRd
computeDiffVar.bin.site <- function(X,inds.g1,inds.g2,
                                    variability.method=rnb.getOption("differential.variability.method"),
                                    paired=FALSE, adjustment.table=NULL, eps=0.01){
  p.vals.var <- rep(as.double(NA),nrow(X))
  do.p.vals <- ncol(X[,inds.g1]) > 1 || ncol(X[inds.g2]) > 1
  if (do.p.vals) {
    if (variability.method == "diffVar"){
      logger.info("Conducting differential variability using diffVar")
      tryCatch(
        p.vals.var <- diffVar(X,inds.g1,inds.g2,adjustment.table=adjustment.table,paired = paired),
        error = function(ee) {
          logger.warning(c("Could not compute p-values using diffVar:",ee$message))
          logger.completed()
        }
      )
    } else if (variability.method == "iEVORA"){
      if(paired){
        logger.warning("Cannot conduct paired variability analysis with iEVORA, changing to diffVar.")
        rnb.options("differential.variability.method"="diffVar")
        tryCatch(
          p.vals.var <- diffVar(X,inds.g1,inds.g2,adjustment.table=adjustment.table,paired=paired),
          error = function(ee) {
            logger.warning(c("Could not compute p-values using diffVar:",ee$message))
            logger.completed()
          }
        )
      }else{
        tryCatch(
          p.vals.var <- apply.iEVORA(X,inds.g1,inds.g2),
          error = function(ee) {
            logger.warning(c("Could not compute p-values using iEVORA:",ee$message))
            logger.completed()
          }
        )
      }
    }
  } else {
    logger.warning("Skipping p-value computation due to insufficient sample numbers")
  }
  p.vals.is.na <- is.na(p.vals.var)
  if (!all(p.vals.is.na)){
    if (any(p.vals.is.na)){
      logger.info(c(sum(p.vals.is.na),"p-values are NA. They are treated as 1 in FDR adjustment"))
      p.vals.var[is.na(p.vals.t.na.adj)] <- 1
    }
    p.vals.var.adj <- p.adjust(p.vals.var, method = "fdr")
  } else {
    p.vals.var.adj <- rep(NA,length(p.vals.var))
  }
  tab.g1 <- X[,inds.g1,drop=FALSE]
  tab.g2 <- X[,inds.g2,drop=FALSE]
  if(length(inds.g1)<2) {
    logger.info("Group 1 has less than 2 members")
    tab.g1 <- as.matrix(tab.g1)
  }
  if(length(inds.g2)<2) {
    logger.info("Group 2 has less than 2 members")
    tab.g2 <- as.matrix(tab.g2)
  }
  
  var.g1 <- apply(tab.g1,1,var)
  var.g2 <- apply(tab.g2,1,var)
  if(paired){
    var.diff <- apply(tab.g1 - tab.g2,1,var)
    var.log.ratio <- apply(X,1,function(x){
      log2(var(x[inds.g1]+eps)/var(x[inds.g2]+eps))
    })
  }else{
    var.diff <- var.g1-var.g2
    var.log.ratio <- ifelse(var.g1==0|var.g2==0,1,log2(var.g1/var.g2))
  }
  neg.log10.p <- -log10(p.vals.var)
  neg.log10.fdr <- -log10(p.vals.var.adj)
  tt <- data.frame(var.g1=var.g1,var.g2=var.g2,var.diff=var.diff,var.log.ratio=var.log.ratio,diffVar.p.val=p.vals.var,diffVar.p.adj.fdr=p.vals.var.adj,log10P=neg.log10.p,
                   log10FDR=neg.log10.fdr)
}

#' computeDiffVar.default.region
#' @noRd
computeDiffVar.default.region <- function(dmtp,regions2sites){
  n.regs.with.sites <- length(regions2sites)
  col.id.var.g1 <- "var.g1"
  col.id.var.g2 <- "var.g2"
  col.id.diff.var <- "var.diff"
  col.id.quot.var <- "var.log.ratio"
  col.id.p.var <- "diffVar.p.val"
  mean.var.g1 <- rep(NA,n.regs.with.sites)
  mean.var.g2 <- rep(NA,n.regs.with.sites)
  diff.var <- rep(NA,n.regs.with.sites)
  quot.var <- rep(NA,n.regs.with.sites)
  p.vals.var <- rep(NA,n.regs.with.sites)
  col.vec <- c(col.id.var.g1, col.id.var.g2, col.id.diff.var, col.id.quot.var, col.id.p.var)
  
  dmt4fastProc <- dmtp[,col.vec]
  if(parallel.isEnabled()){
    dm <- foreach(i=1:n.regs.with.sites, .combine='rbind',.multicombine=TRUE,.maxcombine=200) %dopar% {
      pids <- regions2sites[[i]]
      subtab <- dmt4fastProc[pids,]#these lookups take up most of the time
      
      mean.var.g1   <- mean(subtab[,1],na.rm=TRUE)
      mean.var.g2   <- mean(subtab[,2],na.rm=TRUE)
      diff.var      <- mean(subtab[,3],na.rm=TRUE)
      quot.var      <- mean(subtab[,4],na.rm=TRUE)
      res <- combineTestPvalsMeth(na.omit(subtab[,5]),correlated=TRUE)
      p.vals.var <- NA
      if (length(res)>0) p.vals.var <- res
      c(mean.var.g1,mean.var.g2,diff.var,quot.var,p.vals.var)
    }
    mean.var.g1                  <- dm[,1]
    mean.var.g2                  <- dm[,2]
    diff.var                     <- dm[,3]
    quot.var                     <- dm[,4]
    p.vals.var                  <- dm[,5]
    } else {
      dummy <- sapply(1:n.regs.with.sites,FUN=function(i){
        pids <- regions2sites[[i]]
        subtab <- dmt4fastProc[pids,]#these lookups take up most of the time
        
        mean.var.g1[i] <<- mean(subtab[,1],na.rm=TRUE)
        mean.var.g2[i] <<- mean(subtab[,2],na.rm=TRUE)
        diff.var[i] <<- mean(subtab[,3],na.rm=TRUE)
        quot.var[i] <<- mean(subtab[,4],na.rm=TRUE)
        res <- combineTestPvalsMeth(na.omit(subtab[,5]),correlated=TRUE)
        if (length(res)>0) p.vals.var[i]  <<- res
        return(TRUE)
      })
    }
  p.vals.var.na.adj <- p.vals.var
  p.vals.var.is.na <- is.na(p.vals.var)
  if (any(p.vals.var.is.na)){
    logger.info(c(sum(p.vals.var.is.na),"p-values are NA. They are treated as 1 in FDR adjustment"))
    p.vals.var.na.adj[is.na(p.vals.var.na.adj)] <- 1
  }
  p.vals.var.adj <- p.adjust(p.vals.var.na.adj, method = "fdr")
  tt <- cbind(data.frame(mean.var.g1=mean.var.g1,mean.var.g2=mean.var.g2,mean.var.diff=diff.var,
                            mean.var.log.ratio=quot.var,comb.p.val.var=p.vals.var,comb.p.adj.var.fdr=p.vals.var.adj))
  return(tt)
}

#' computeDiffVar.bin.region
#' @noRd
computeDiffVar.bin.region <- function(rnbSet,dmtp,inds.g1,inds.g2,region.types=rnb.region.types(assembly(rnbSet)), ...){
  if (length(union(inds.g1,inds.g2)) != (length(inds.g1)+length(inds.g2))){
    logger.error("Overlapping sample sets in differential methylation analysis")
  }
  logger.start('Computing Differential Variability Tables (Region Level)')

  diffmeth.tabs <- list()
  for (rt in region.types){
    regions2sites <- regionMapping(rnbSet,rt)
    regions2sites.is.all.na <- sapply(regions2sites,FUN=function(x){all(is.na(x))})
    if (any(regions2sites.is.all.na)) {
      stop(paste("Region mapping of RnBSet from sites to regions is inconsistent (",rt,")"))
    }
    dmtr <- computeDiffVar.default.region(dmtp,regions2sites)
    diff.var.ranks <- cols.to.rank.region(dmtr)
    comb.rank.var <- combinedRanking.tab(diff.var.ranks,rerank=FALSE)
    dmtr$combinedRank.var <- comb.rank.var
    
    diffmeth.tabs <- c(diffmeth.tabs,list(dmtr))
    logger.status(c("Computed table for", rt))
  }
  names(diffmeth.tabs) <- region.types
  logger.completed()
  return(diffmeth.tabs)
}
