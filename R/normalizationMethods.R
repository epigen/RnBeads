########################################################################################################################
## normalizationMethods.R
## created: 2021-02-04
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Implementation of the normalization methods.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

# reimplementation of the basic scaling normalization method for dye bias correction
#
# method='reference', default Illumina scaling normalization to a reference sample, taken from methylumi
# method='internal', singe-sample correlction, taken from minfi
#

rnb.norm.scaling<-function(rnb.set, method=c("internal", "reference")[1L], ref.sample=NULL){
    
    if(!method %in% c("reference", "internal")){
        rnb.error("Unsupported value for method")
    }
    
    qc_int<-qc(rnb.set)
    ctrl_target<-gsub("probes", "controls", rnb.set@target)
    qc_annot<-rnb.get.annotation(ctrl_target, rnb.set@assembly)
    annot<-annotation(rnb.set)  
    
    channels<-c("Cy5", "Cy3")
    tokens<-get.platform.tokens(rnb.set@target)
    
    get_norm_int<-function(ch){
        ids<-qc_annot[[tokens$id_col]][grep(tokens$norm_token[ch], qc_annot[[tokens$trg_col]])]
        ids<-as.character(ids)
        qc_int[[ch]][ids,]
    }
    
    norm_int_list<-lapply(channels, get_norm_int)
    names(norm_int_list)<-channels
    
    #taken from methylumi
    Grn.avg <- colMeans(norm_int_list[["Cy3"]])
    Red.avg <- colMeans(norm_int_list[["Cy5"]])
    R.G.ratio = Red.avg/Grn.avg
    if(method=="reference"){
        if(is.null(ref.sample)) ref.sample <- which.min( abs(R.G.ratio-1) )
        rnb.info(paste('Using sample number', ref.sample, 'as reference level...'))
        
        ref <- (Grn.avg + Red.avg)[ref.sample]/2
        if(is.na(ref)) rnb.error("'reference' refers to an array that is not present")
        Grn.factor <- ref/Grn.avg
        Red.factor <- ref/Red.avg
    }else if(method=="internal"){
        Grn.factor <- rep(1, length(R.G.ratio))
        Red.factor <- 1/R.G.ratio
    }
    
    probe.categories<-INTENSITY.SUMMARIZATION.INFO[[rnb.set@target]]

    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        M_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@M)[1], ncol=dim(rnb.set@M)[2]))
        U_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@U)[1], ncol=dim(rnb.set@U)[2]))
        M0_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@M0)[1], ncol=dim(rnb.set@M0)[2]))
        U0_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@U0)[1], ncol=dim(rnb.set@U0)[2]))
    }else{
        M_temp<-matrix(NA_real_, nrow=nrow(rnb.set@M), ncol=ncol(rnb.set@M))
        U_temp<-matrix(NA_real_, nrow=nrow(rnb.set@U), ncol=ncol(rnb.set@U))
        M0_temp<-matrix(NA_real_, nrow=nrow(rnb.set@M0), ncol=ncol(rnb.set@M0))
        U0_temp<-matrix(NA_real_, nrow=nrow(rnb.set@U0), ncol=ncol(rnb.set@U0))
    }
    
    for(pc in probe.categories){
        
        if(pc$Color=="Both"){
            M.factor<-Grn.factor
            U.factor<-Red.factor
            M0.factor<-Red.factor
            U0.factor<-Grn.factor
        }else if(pc$Color=="Red"){
            M.factor<-Red.factor
            U.factor<-Red.factor
            M0.factor<-Grn.factor
            U0.factor<-Grn.factor
        }else if(pc$Color=="Grn"){
            M.factor<-Grn.factor
            U.factor<-Grn.factor
            M0.factor<-Red.factor
            U0.factor<-Red.factor
        }
        indices<-which(annot$Design==pc$Design & annot$Color==pc$Color)
        
        M_temp[indices,]<-sweep(rnb.set@M[indices,], 2, FUN="*", M.factor)
        U_temp[indices,]<-sweep(rnb.set@U[indices,], 2, FUN="*", U.factor)
        M0_temp[indices,]<-sweep(rnb.set@M0[indices,], 2, FUN="*", M0.factor)
        U0_temp[indices,]<-sweep(rnb.set@U0[indices,], 2, FUN="*", U0.factor)
        
    }
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump && isTRUE(rnb.set@status$discard.ff.matrices)){
        delete(rnb.set@M)
        delete(rnb.set@U)
        delete(rnb.set@M0)
        delete(rnb.set@U0)
    }
    
    rnb.set@M<-M_temp
    rnb.set@U<-U_temp
    rnb.set@M0<-M0_temp
    rnb.set@U0<-U0_temp
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        rm(M_temp); rm(U_temp); rm(M0_temp); rm(U0_temp); rnb.cleanMem()
    }
    
    rnb.set<-update.meth(rnb.set)
    
    return(rnb.set)
}



# reimplementation of the basic subtraction  normalization method

rnb.bgcorr.subtr<-function(rnb.set, bg.quant=0.05, offset=0){
    
    qc_int<-qc(rnb.set)
    ctrl_target<-gsub("probes", "controls", rnb.set@target)
    qc_annot<-rnb.get.annotation(ctrl_target, rnb.set@assembly)
    annot<-annotation(rnb.set)
    
    channels<-c("Cy5", "Cy3")
    tokens<-get.platform.tokens(rnb.set@target)
    
    get_bg_int<-function(ch){
        ids<-qc_annot[[tokens$id_col]][grep(tokens$bg_token, qc_annot[[tokens$trg_col]])]
        ids<-as.character(ids)
        qc_int[[ch]][ids,]
    }
    
    bg_int_list<-lapply(channels, get_bg_int)
    names(bg_int_list)<-channels
    
    bg<-list()
    for(ch in channels) {
        bg[[ch]]<-colQuantiles(bg_int_list[[ch]], probs=bg.quant)
    }
    
    probe.categories<-INTENSITY.SUMMARIZATION.INFO[[rnb.set@target]]
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        M_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@M)[1], ncol=dim(rnb.set@M)[2]))
        U_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@U)[1], ncol=dim(rnb.set@U)[2]))
    }else{
        M_temp<-matrix(NA_real_, nrow=nrow(rnb.set@M), ncol=ncol(rnb.set@M))
        U_temp<-matrix(NA_real_, nrow=nrow(rnb.set@U), ncol=ncol(rnb.set@U))
    }
    
    for(pc in probe.categories){
        
        if(pc$Color=="Both"){
            bg.M<-bg[["Cy5"]]
            bg.U<-bg[["Cy3"]]
        }else if(pc$Color=="Red"){
            bg.M<-bg[["Cy3"]]
            bg.U<-bg[["Cy3"]]
        }else if(pc$Color=="Grn"){
            bg.M<-bg[["Cy5"]]
            bg.U<-bg[["Cy5"]]
        }
        indices<-which(annot$Design==pc$Design & annot$Color==pc$Color)
        for(i in 1:ncol(rnb.set@M)) {
            M_temp[indices,i] <- pmax(rnb.set@M[indices,i] - bg.M[i], 1)
            U_temp[indices,i] <- pmax(rnb.set@U[indices,i] - bg.U[i], 1)
        }
    }
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump && isTRUE(rnb.set@status$discard.ff.matrices)){
        delete(rnb.set@M)
        delete(rnb.set@U)
    }
    
    rnb.set@M<-M_temp
    rnb.set@U<-U_temp
    
    if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
        rm(M_temp); rm(U_temp); rnb.cleanMem()
    }
    
    rnb.set<-update.meth(rnb.set)
    
    return(rnb.set)
    
}

#
# NOOB method as implemented in SeSAME
#
noob.sesame <- function(rnb.set, background=c("internal", "bleed-through")[2L], offset=15) {
    
    annot<-annotation(rnb.set)
    
    ## if no Infinium-I probes
    #if (nrow(InfIG(sdf)) == 0 && nrow(InfIR(sdf)) == 0) { return(sdf) }
    
    ##################### GRN CHANNEL ESTIMATION
    ## background
    oobG <- rbind(
            rnb.set@M0[annot$Design=="I" & annot$Color=="Red",], 
            rnb.set@U0[annot$Design=="I" & annot$Color=="Red",])
    oobG[is.na(oobG)] <-0
    #oobG[oobG == 0] <- 1
    
    ## if not enough out-of-band signal
    if (sum(oobG > 0, na.rm=TRUE) < 100) { return(rnb.set) }
    
    ## foreground
    #ibG = c(InfIG(sdf)$MG, InfIG(sdf)$UG, InfII(sdf)$UG)
    ibG <- rbind(
            rnb.set@M[annot$Design=="I" & annot$Color=="Grn",], 
            rnb.set@U[annot$Design=="I" & annot$Color=="Grn",]
            )

    if(background == "internal") ibG <- rbind(ibG, rnb.set@M[annot$Design=="II",])

    ibG[is.na(ibG)] <- 0
    #ibG[ibG == 0] <- 1 # set signal to 1 if 0
    
   
    ##################### RED CHANNEL ESTIMATION
    oobR <- rbind(
            rnb.set@M0[annot$Design=="I" & annot$Color=="Grn",], 
            rnb.set@U0[annot$Design=="I" & annot$Color=="Grn",])
    oobR[is.na(oobR)] <-0
    #oobR[oobR == 0] <- 1
    
    if (sum(oobR > 0, na.rm=TRUE) < 100) { return(rnb.set) }
    
    ibR <- rbind(
            rnb.set@M[annot$Design=="I" & annot$Color=="Red",], 
            rnb.set@U[annot$Design=="I" & annot$Color=="Red",]
            )
    
    if(background == "internal") ibR <- rbind(ibR, rnb.set@U[annot$Design=="II",])
    
    ibR[is.na(ibR)] <- 0
    #ibR[ibR == 0] <- 1 # set signal to 1 if 0
    
    if(background=="internal"){
        
        #fitG = backgroundCorrectionNoobFit(ibG, oobG)
        oobG[oobG == 0] <- 1
        fitG <- sapply(seq(ncol(ibG)), function(i) backgroundCorrectionNoobFit(ibG[,i],oobG[,i]))
        
        #fitR = backgroundCorrectionNoobFit(ibR, oobR)
        oobR[oobR == 0] <- 1
        fitR <- sapply(seq(ncol(ibR)), function(i) backgroundCorrectionNoobFit(ibR[,i],oobR[,i]))
        
    }else if(background=="bleed-through"){
        
        ibR_other<-rbind(oobG, rnb.set@M[annot$Design=="II",])
        ibR_other[is.na(ibR_other)] <- 0
        
        ibG_other<-rbind(oobR, rnb.set@U[annot$Design=="II",])
        ibG_other[is.na(ibG_other)] <- 0
        
        oobR[oobR == 0] <- 1
        bg.GpredictR <- lapply(seq(ncol(ibG)), function(i) {f<-train.model.lm(ibG[,i], oobR[,i]); return(f)})
        ##ibR.other.channel
        #ibR_other[ibR_other==0] <- 1
        pp.bg.ibR<-sapply(seq(ncol(ibR_other)), function(i) bg.GpredictR[[i]](ibR_other[,i]))
        ibR_other<-NULL; rm(ibR_other)
        ##sset$IG
        pp.bg.oobR<-sapply(seq(ncol(ibG)), function(i) bg.GpredictR[[i]](ibG[,i]))
        
        ibRR<-rbind(ibR, rnb.set@U[annot$Design=="II",])
        ibRR[ibRR==0] <- 1
        alphaR<-sapply(seq(ncol(ibRR)), function(i) pmax(MASS::huber(ibRR[,i] - pp.bg.ibR[,i]$mu)$mu, offset-5))
        ibRR<-NULL; rm(ibRR)
        
        oobG[oobG == 0] <- 1
        bg.RpredictG <- lapply(seq(ncol(ibR)), function(i) {f<-train.model.lm(ibR[,i], oobG[,i]); return(f)})
        ##ibG.other.channel

        #ibG_other[ibG_other==0] <- 1
        pp.bg.ibG<-sapply(seq(ncol(ibG_other)), function(i) bg.RpredictG[[i]](ibG_other[,i]))
        ibG_other<-NULL; rm(ibG_other)
        ##sset$IR
        pp.bg.oobG<-sapply(seq(ncol(ibR)), function(i) bg.RpredictG[[i]](ibR[,i]))
        
        ibGG<-rbind(ibG, rnb.set@M[annot$Design=="II",])
        ibGG[ibGG==0] <- 1
        alphaG<-sapply(seq(ncol(ibGG)), function(i) pmax(MASS::huber(ibGG[,i] - pp.bg.ibG[,i]$mu)$mu, offset-5))
        ibGG<-NULL; rm(ibGG)

    }
    ####### CORRECTION
    tII.len<-sum(annot$Design=="II")
    mg.ind<-c(which(annot$Color=="Grn"), which(annot$Design=="II"))
    ug.ind<-which(annot$Color=="Grn" )
    
    mr.ind<-which(annot$Color=="Red" )
    ur.ind<-c(which(annot$Color=="Red"), which(annot$Design=="II"))
    
    slot_names<-c("M", "U", "M0", "U0")
    ind.lists<-list("Grn"=list(), "Red"=list())
    
    ind.lists[["Grn"]][["M"]]<-mg.ind
    ind.lists[["Grn"]][["U"]]<-ug.ind
    ind.lists[["Grn"]][["M0"]]<-mr.ind
    ind.lists[["Grn"]][["U0"]]<-mr.ind
    
    ind.lists[["Red"]][["M"]]<-mr.ind
    ind.lists[["Red"]][["U"]]<-ur.ind
    ind.lists[["Red"]][["M0"]]<-ug.ind
    ind.lists[["Red"]][["U0"]]<-ug.ind
    
    if(background=="internal"){
        
        fit_param<-list("Grn"=fitG, "Red"=fitR)
        
    }else if(background=="bleed-through"){
        
        ir<-expand.grid(seq(nrow(pp.bg.ibG)), seq(ncol(pp.bg.ibG)))
        
        pp.bg.ibG_M<-pp.bg.ibG
        for(ind in 1:nrow(ir)) pp.bg.ibG_M[[ir[ind,1],ir[ind,2]]]<-pp.bg.ibG_M[[ir[ind,1],ir[ind,2]]][c(seq(ug.ind),(2*length(ug.ind)+1):(2*length(ug.ind)+tII.len))]
        pp.bg.ibG_U<-pp.bg.ibG
        for(ind in 1:nrow(ir)) pp.bg.ibG_U[[ir[ind,1],ir[ind,2]]]<-pp.bg.ibG_U[[ir[ind,1],ir[ind,2]]][c((length(ug.ind)+1):(2*length(ug.ind)))]
        pp.bg.ibR_M<-pp.bg.ibR
        for(ind in 1:nrow(ir)) pp.bg.ibR_M[[ir[ind,1],ir[ind,2]]]<-pp.bg.ibR_M[[ir[ind,1],ir[ind,2]]][c(seq(mr.ind))]
        pp.bg.ibR_U<-pp.bg.ibR
        for(ind in 1:nrow(ir)) pp.bg.ibR_U[[ir[ind,1],ir[ind,2]]]<-pp.bg.ibR_U[[ir[ind,1],ir[ind,2]]][c((length(mr.ind)+1):(2*length(mr.ind)+tII.len))]

        pp.bg.oobG_M<-pp.bg.oobG
        for(ind in 1:nrow(ir)) pp.bg.oobG_M[[ir[ind,1],ir[ind,2]]]<-pp.bg.oobG_M[[ir[ind,1],ir[ind,2]]][seq(mr.ind)]
        pp.bg.oobG_U<-pp.bg.oobG
        for(ind in 1:nrow(ir)) pp.bg.oobG_U[[ir[ind,1],ir[ind,2]]]<-pp.bg.oobG_U[[ir[ind,1],ir[ind,2]]][(length(mr.ind)+1):(2*length(mr.ind))]
        pp.bg.oobR_M<-pp.bg.oobR
        for(ind in 1:nrow(ir)) pp.bg.oobR_M[[ir[ind,1],ir[ind,2]]]<-pp.bg.oobR_M[[ir[ind,1],ir[ind,2]]][seq(ug.ind)]
        pp.bg.oobR_U<-pp.bg.oobR
        for(ind in 1:nrow(ir)) pp.bg.oobR_U[[ir[ind,1],ir[ind,2]]]<-pp.bg.oobR_U[[ir[ind,1],ir[ind,2]]][(length(ug.ind)+1):(2*length(ug.ind))]
        
        fit_param<-list(
                "Grn"=list(
                        "M"=list(
                            "IB"=pp.bg.ibG_M,
                            "OB"=pp.bg.oobG_M
                        ),
                        "U"=list(
                             "IB"=pp.bg.ibG_U,
                             "OB"=pp.bg.oobG_U
                        ),
                        "alpha"=alphaG
                        ), 
                "Red"=list(
                        "M"=list(
                            "IB"=pp.bg.ibR_M,
                            "OB"=pp.bg.oobR_M
                        ),
                        "U"=list(
                            "IB"=pp.bg.ibR_U,
                            "OB"=pp.bg.oobR_U
                        ),
                        "alpha"=alphaR
                        )
        )
    }
    
    for(sl in slot_names){
        if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
            mat_temp<-convert.to.ff.matrix.tmp(matrix(NA_real_, nrow=dim(rnb.set@M)[1], ncol=dim(rnb.set@M)[2]))
        }else{
            mat_temp<-matrix(NA_real_, nrow=nrow(rnb.set@M), ncol=ncol(rnb.set@M))
        }
        for(col in c("Grn", "Red")){
            for(i in seq(ncol(rnb.set@M))){
                if(background=="internal"){
                    mat_temp[ind.lists[[col]][[sl]], i]<-normExpSignal(
                            fit_param[[col]][,i]$mu, 
                            fit_param[[col]][,i]$sigma, 
                            fit_param[[col]][,i]$alpha, 
                            slot(rnb.set,sl)[ind.lists[[col]][[sl]],i]) + offset
                }else{
                    token_mu<-gsub("0","",sl)
                    token_band<-c("IB", "OB")[1L+grepl("0", sl)]
                    input<-slot(rnb.set,sl)[ind.lists[[col]][[sl]],i]
                    input[is.na(input) | input == 0]<-1L
                    mat_temp[ind.lists[[col]][[sl]], i]<-backgroundCorrCh1(
                            fit_param[[col]][[token_mu]][[token_band]][,i],
                            fit_param[[col]]$alpha[i], 
                            input) + offset
                }
            }
        }
        if(!is.null(rnb.set@status) && rnb.set@status$disk.dump && isTRUE(rnb.set@status$discard.ff.matrices)){
            delete(slot(rnb.set,sl))
        }
        
        slot(rnb.set,sl)<-mat_temp
        
        if(!is.null(rnb.set@status) && rnb.set@status$disk.dump){
            rm(mat_temp); rnb.cleanMem()
        }
    }
    
    rnb.set<-update.meth(rnb.set)
    
    return(rnb.set)
}

backgroundCorrectionNoobFit <- function(ib, bg) {
    e = MASS::huber(bg)
    mu = e$mu
    sigma = e$s
    alpha = pmax(MASS::huber(ib)$mu-mu, 10)
    list(mu = mu, sigma = sigma, alpha = alpha)
}

# the following is adapted from Limma
## normal-exponential deconvolution (conditional expectation of
## xs|xf; WEHI code)
normExpSignal <- function(mu, sigma, alpha, x)  {
    
    sigma2 <- sigma * sigma
    
    if (alpha <= 0)
        stop("alpha must be positive")
    if (sigma <= 0)
        stop("sigma must be positive")
    
    mu.sf <- x - mu - sigma2/alpha
    signal <- mu.sf + sigma2 * exp(
            dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
                    pnorm(
                            0, mean = mu.sf, sd = sigma,
                            lower.tail = FALSE, log.p = TRUE))
    
    o <- !is.na(signal)
    if (any(signal[o] < 0)) {
        warning("Limit of numerical accuracy reached with very
                        low intensity or very high background:\nsetting adjusted intensities
                        to small value")
        signal[o] <- pmax(signal[o], 1e-06)
    }
    signal
}

#####

train.model.lm <- function(input, output) {
    fitdata <- data.frame(IB=as.vector(input)+1, LOB=log(as.vector(output)+1))
    m <- lm(LOB~IB, data=fitdata)
    #m <- MASS::rlm(LOB~IB, data=fitdata)
    function(d) {
        force(m)
        pp <- predict(m, newdata=data.frame(IB=as.vector(d)), interval='prediction', level=0.8)
        # use upper bound for mu since true signal is often much higher than noise
        list(mu=exp(pp[,'fit']), sigma=(exp(pp[,'upr'])-exp(pp[,'lwr']))/10.13)
        # list(mu=exp(pp[,'upr']), sigma=log(exp(pp[,'upr'])-exp(pp[,'lwr'])))
    }
}

backgroundCorrCh1 <- function(pp, alpha, x) {
    mu.bg <- pp$mu
    sigma <- pp$sigma
    sigma2 <- sigma * sigma
    if (alpha <= 0)
        stop("alpha must be positive")
    if (any(na.omit(sigma) <= 0))
        stop("sigma must be positive")
    mu.sf <- x - mu.bg - sigma2/alpha
    signal <- mu.sf + sigma2 * exp(
            dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
                    pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log.p = TRUE))
    o <- !is.na(signal)
    if (any(signal[o] < 0)) {
        warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
        signal[o] <- pmax(signal[o], 1e-06)
    }
    signal.min <- min(signal, na.rm = TRUE)
    signal <- signal - signal.min
    signal
}

get.platform.tokens<-function(platform){
    
    channels<-c("Cy5", "Cy3")
    token<-setNames(rep(NA_character_, 2), channels)
    
    dict<-list()
    
    if(platform=="probes27"){
        dict$norm_token["Cy3"] <- 'Norm.G'
        dict$norm_token["Cy5"] <- 'Norm.R'
        dict$bg_token<-"Negative"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    } else if(platform=="probes450"){
        dict$norm_token["Cy3"] <- 'NORM_(C|G)'
        dict$norm_token["Cy5"] <- 'NORM_(A|T)'
        dict$bg_token<-"NEGATIVE"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    } else if(platform=="probesEPIC"){
        dict$norm_token["Cy3"] <- 'NORM_(C|G)'
        dict$norm_token["Cy5"] <- 'NORM_(A|T)'
        dict$bg_token<-"NEGATIVE"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    }else if(platform=="probesMMBC"){
        dict$norm_token["Cy3"] <- 'NORM_(C|G)'
        dict$norm_token["Cy5"] <- 'NORM_(A|T)'
        dict$bg_token<-"NEGATIVE"
        dict$id_col<-"ID"
        dict$trg_col<-"Target"
    }
    
    return(dict)
}

