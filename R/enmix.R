########################################################################################################################
## enmix.R
## created: 2016-10-02
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions for performing exponential-normal (EN) background subtraction on 450k, EPIC and EPICv2 datasets.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.probe.types
#'
#' Gets the indices of the probes of specific type in the given dataset.
#'
#' @param dataset     Infinium 450k, EPIC or EPICv2 dataset as an instance of \code{\linkS4class{RnBeadRawSet}}.
#' @param include.snp Flag indicating if the SNP probes, if present in the dataset, are to be considered.
#' @return Probe indices, split by types as a \code{list} with three elements: \code{"IGrn"}, \code{"IRed"},
#'         \code{"II"}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.probe.types <- function(dataset, include.snp = TRUE) {
	ptypes <- rnb.get.annotation(dataset@target, dataset@assembly)
	pindices <- tapply(dataset@sites[, 3], dataset@sites[, 2], unname, simplify = FALSE)
	i.chrom <- as.integer(names(pindices))
	probe.types <- list()
	for (i in 1:length(pindices)) {
		probe.types[[i]] <- as.data.frame(mcols(ptypes[[i.chrom[i]]])[pindices[[i]], c("Design", "Color")])
	}
	probe.types <- do.call(rbind, probe.types)
	rm(pindices, i.chrom)
	i.probes <- list()
	for (probe.design in levels(probe.types[, 1])) {
		i <- probe.types[, 1] == probe.design
		if (probe.design == "I") {
			for (probe.color in levels(probe.types[, 2])[-1]) {
				i.probes[[paste0(probe.design, probe.color)]] <-
					which(i & probe.types[, 2] == probe.color)
			}
		} else { # probe.design == "II"
			i.probes[[probe.design]] <- which(i)
		}
	}
	if (!include.snp) {
		is.snp <- which(unlist(lapply(ptypes, function(x) { grepl("^rs", names(x)) }), use.names = FALSE))
		i.probes <- lapply(i.probes, setdiff, y = is.snp)
	}
	i.probes
}

########################################################################################################################

#' rnb.enmix.oob
#' 
#' Runs the Exponential-truncated-normal (ENmix) background subtraction using out-of-band (OOB) intensities to estimate
#' the background (noise) distribution parameters.
#' 
#' @param dataset      Infinium 450k, EPIC or EPICv2 dataset as an instance of \code{\linkS4class{RnBeadRawSet}}.
#' @param update.betas Flag indicating if the beta value matrix is to be updated after correcting the M and U intensity
#'                     values.
#' @return The modified dataset.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.enmix.oob <- function(dataset, update.betas = TRUE) {
	ptypes <- rnb.probe.types(dataset)
	if (parallel.isEnabled()) {
		## TODO: Update fexport after integrating this in RnBeads
		fexport <- c("rnb.enmix.em", "rnb.enmix.cm", "rnb.enmix.adjust", "rnb.enmix.oob.sample")
		results <- foreach(i = 1:nrow(dataset@pheno), .packages	= c("RnBeads", "MASS"), .export = fexport) %dopar%
			rnb.enmix.oob.sample(dataset, i, ptypes = ptypes)
		for (i in 1:length(results)) {
			dataset@M[, i] <- results[[i]][[1]]
			dataset@U[, i] <- results[[i]][[2]]
		}
	} else {
		for (i in 1:nrow(dataset@pheno)) {
			result <- rnb.enmix.oob.sample(dataset, i, ptypes)
			dataset@M[, i] <- result[[1]]
			dataset@U[, i] <- result[[2]]
		}
	}
	if (update.betas) {
		if (dataset@status$disk.dump) {
			dataset@meth.sites <- convert.to.ff.matrix.tmp(beta.value(dataset@M[,], dataset@U[,]))
		} else {
			dataset@meth.sites <- beta.value(dataset@M[,], dataset@U[,])
		}
	}
	dataset
}

########################################################################################################################

#' rnb.enmix.oob.sample
#'
#' Processes a single sample from the given dataset using ENmix.
#'
#' @param dataset  Infinium 450k, EPIC or EPICv2 dataset as an instance of \code{\linkS4class{RnBeadRawSet}}.
#' @param i.sample Index of the sample to be processed.
#' @param ptypes   Probe indices, split by types, as returned by \code{\link{rnb.probe.types}}.
#' @return Corrected intensity values for the i-th sample as a \code{list} of two \code{vectors}: one of the M, and
#'         one for the U values, respectively.
#'
#' @author Yassen Assenov
#' @noRd
rnb.enmix.oob.sample <- function(dataset, i.sample, ptypes = rnb.probe.types(dataset, FALSE)) {
	
	## Construct the resulting vectors
	Mhat <- dataset@M[, i.sample]
	Uhat <- dataset@U[, i.sample]
	
	## Compute estimates mu_B and s_B for the green channel
	oob <- c(dataset@M0[ptypes$IRed, i.sample], dataset@U0[ptypes$IRed, i.sample])
	oob[is.na(oob)] <- 10^(-6) # FIXME: these intensities should be set to zero upon loading
	oob[oob <= 0] <- 10^(-6)
	estimates <- unlist(MASS::huber(oob), use.names = FALSE)
	
	## Fit the green channel intensities and adjust to the fit
	Mhat[ptypes$IGrn] <- rnb.enmix.adjust(Mhat[ptypes$IGrn], estimates[1], estimates[2])
	Uhat[ptypes$IGrn] <- rnb.enmix.adjust(Uhat[ptypes$IGrn], estimates[1], estimates[2])
	Mhat[ptypes$II] <- rnb.enmix.adjust(Mhat[ptypes$II], estimates[1], estimates[2])
	
	## Compute estimates mu_B and s_B for the red channel
	oob <- c(dataset@M0[ptypes$IGrn, i.sample], dataset@U0[ptypes$IGrn, i.sample])
	oob[is.na(oob)] <- 10^(-6) # FIXME: these intensities should be set to zero upon loading
	oob[oob <= 0] <- 10^(-6)
	estimates <- unlist(huber(oob), use.names = FALSE)
	
	## Fit the red channel intensities and adjust to the fit
	Mhat[ptypes$IRed] <- rnb.enmix.adjust(Mhat[ptypes$IRed], estimates[1], estimates[2])
	Uhat[ptypes$IRed] <- rnb.enmix.adjust(Uhat[ptypes$IRed], estimates[1], estimates[2])
	Uhat[ptypes$II] <- rnb.enmix.adjust(Uhat[ptypes$II], estimates[1], estimates[2])
	
	return(list(Mhat, Uhat))
}

########################################################################################################################

#' rnb.enmix.adjust
#' 
#' Adjusts the intensity values of the given sample using Exponential-truncated-normal (EN) mixture.
#' 
#' @param x      Vector of intensity values to be adjusted.
#' @param muB    Estimated mu parameter for the background (noise).
#' @param sigmaB Estimated sigma parameter for the background (noise).
#' @return Corrected vector of intensity values.
#' 
#' @author Zongli Xu and Liang Niu, modified by Yassen Assenov
#' @noRd
rnb.enmix.adjust <- function(x, muB, sigmaB) {
	x[is.na(x)] <- 1
	x[x < 1] <- 1
	
	xx <- x[x >= muB] - muB
	estimates <- rnb.enmix.em(xx)
	lambda <- estimates[1]
	muS <- estimates[2]
	sigmaS <- ifelse(estimates[3] <= sigmaB, 0.1, sqrt(estimates[3]^2 - sigmaB^2))
	p <- (length(x) - length(xx) + estimates[4] * length(xx)) / length(x)
	
	result <- rnb.enmix.cm(lambda, muS, sigmaS, p, muB, sigmaB, x)
	result[result <= 0.01] <- 0.01
	return(result)
}

########################################################################################################################

#' rnb.enmix.em
#'
#' Runs the Expectation Maximization (EM) algorithm to fit a model: p * exp(lambda) + (1 - p) * N+(muS, sigmaS^2).
#'
#' @param x Set of read values, in the form of a \code{double} \code{vector}, against which the model is to be fit.
#' @return \code{vector} of the fitted parameters: lambda, mu, sigma, and p.
#'
#' @author Zongli Xu and Liang Niu, modified by Yassen Assenov
#' @noRd
rnb.enmix.em <- function(x) {
	epsilon = c(1e-04, 0.001, 0.001, 0.001)
	p <- 0.5
	lambda <- max(density(x)$y)
	mu <- mean(range(x))
	sigma2 <- (diff(range(x)) / 6)^2
	n <- as.double(length(x))
	repeat {
		z <- p / (p + (1 - p) * exp(dnorm(x, mu, sqrt(sigma2), log = TRUE) - dexp(x, lambda, log = TRUE)))
		sz <- sum(z)
		p.new <- sz / n
		lambda.new <- sz / sum(x * z)
		z <- 1 - z; sz <- n - sz
		mu.new <- sum(x * z) / sz
		sigma2.new <- sum(((x - mu.new)^2) * z) / sz
		finished <- abs(c(lambda.new - lambda, mu.new - mu, sigma2.new - sigma2, p.new - p))
		finished <- all(finished < epsilon)
		p <- p.new
		lambda <- lambda.new
		mu <- mu.new
		sigma2 <- sigma2.new
		if (finished) { break }
	}
	return(c(lambda, mu, sqrt(sigma2), p))
}

########################################################################################################################

#' rnb.enmix.cm
#' 
#' Computes corrected signal intensities given the observed intensities and EN model parameters.
#' 
#' @param lambda The lambda parameter.
#' @param muS    The mu parameter for the signal.
#' @param sigmaS The sigma paraemter for the signal.
#' @param p      The p parameter.
#' @param muB    The mu parameter for the background (noise).
#' @param sigmaB The sigma parameter for the background (noise).
#' @param x      Observed intensities given as a \code{vector} of type \code{double}.
#' 
#' @author Zongli Xu and Liang Niu, modified by Yassen Assenov
#' @noRd
rnb.enmix.cm <- function(lambda, muS, sigmaS, p, muB, sigmaB, x) {
	a <- muB + lambda * sigmaB^2
	C <- (sigmaS^2 * (x - muB) + sigmaB^2 * muS) / (sigmaS^2 + sigmaB^2)
	D <- sigmaS * sigmaB / sqrt(sigmaS^2 + sigmaB^2)
	E <- pnorm(x, C, D) - pnorm(0, C, D)
	G <- pnorm(x, a, sigmaB) - pnorm(0, a, sigmaB)
	T1 <- p * lambda * exp(lambda^2 * sigmaB^2/2 - lambda * (x - muB)) * G
	T2 <- (1 - p) * dnorm(x, muS + muB, sqrt(sigmaS^2 + sigmaB^2)) / (1 - pnorm(0, muS, sigmaS))
	D1 <- T1 * (x - a + sigmaB * (dnorm((x - a) / sigmaB) - dnorm(a / sigmaB)) / G)
	D2 <- T2 * (C * E + D * (exp(-C^2/(2 * D^2)) - exp(-(x - C)^2 / (2 * D^2))) / sqrt(2 * pi))
	(D1 + D2) / (T1 + T2 * E)
}
