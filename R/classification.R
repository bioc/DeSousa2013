##get centroids of PAM classifier
getCentroids <- function (fit, data, threshold) {

    genenames <- data$genenames[fit$gene.subset]
    x <- data$x[fit$gene.subset, fit$sample.subset]
    clabs <- colnames(fit$centroids)
    scen <- pamr.predict(fit, data$x, threshold = threshold, 
        type = "cent")
    dif <- (scen - fit$centroid.overall)/fit$sd
    if (!is.null(fit$y)) {
        nc <- length(unique(fit$y))
    }
    if (is.null(fit$y)) {
        nc <- ncol(fit$proby)
    }
    o <- drop(abs(dif) %*% rep(1, nc)) > 0
    d <- dif[o, ]
    nd <- sum(o)
    genenames <- genenames[o]
    xx <- x[o, ]
    oo <- order(apply(abs(d), 1, max))
    d <- d[oo, ]
    genenames <- genenames[oo]
    return(d)
}
##build classifier using 87 samples
buildClassifier <- function(sigMat, clus.f, nfold=10, nboot=100) {
##	sdat <- diffGenes <- diffGenes.f <- clus.f <- sdat.f <- silh <- NULL
##	data("silh", package="DeSousa2013", envir = environment())
##	data("diffGenes", package="DeSousa2013", envir = environment())
##	data("diffGenes.f", package="DeSousa2013", envir = environment())
##	data("dat", package="DeSousa2013", envir = environment())
##	sigMat <- sdat.f[diffGenes.f, names(clus.f)]
	##1. PAM analysis
	dat <- list()
	dat$x <- sigMat[, names(clus.f)]
	dat$y <- rep(3, length(clus.f))
	dat$y[clus.f==2] <- 2
	dat$y[clus.f==1] <- 1
	dat$y <- factor(dat$y)
	dat$geneid <- rownames(dat$x)
	dat$genenames <- rownames(dat$x)
	pam.rslt <- pamr.train(data=dat)

	pam.cv.rslt.l <- list()
	for(i in 1:nboot) {
		pam.cv.rslt <- pamr.cv(fit=pam.rslt, data=dat, nfold=nfold)
		pam.cv.rslt.l[[i]] <- pam.cv.rslt
	}	
	err <- t(sapply(1:length(pam.cv.rslt.l), function(x) 
		pam.cv.rslt.l[[x]]$error, simplify=TRUE))
	ngenes <- c(sapply(1:(length(pam.rslt$threshold)-1), function(x) 
		nrow(pamr.listgenes(pam.rslt, dat, pam.rslt$threshold[x]))
	), 0)
	colnames(err) <- ngenes

	thresh <- pam.rslt$threshold[18]

	##2. select signatures
	signature <- (pamr.listgenes(pam.rslt, dat, thresh))[, "id"]
	##2a. get centroids
	cents <- getCentroids(pam.rslt, dat, thresh)
	cents <- cents[signature, ]	
	
##	save(signature, pam.rslt, thresh, err, cents, file="classifier.RData")
	return(list(signature=signature, pam.rslt=pam.rslt, thresh=thresh, 
		err=err, cents=cents))
}


##9 use the classifier to classify all 90 AMC samples
pamClassify <- function(datsel, signature, pam.rslt, thresh, postRth=1) {
##	sdat <- uniGenes <- diffGenes.f <- signature <- 
##		pam.rslt <- thresh <- err <- cents <- NULL 
##	data("dat", package="DeSousa2013", envir = environment())
##	data("uniGenes", package="DeSousa2013", envir = environment())
##	data("diffGenes.f", package="DeSousa2013", envir = environment())
##	data("classifier", package="DeSousa2013", envir = environment())
	
##	sdat <- sdat[names(uniGenes), ]
##	rownames(sdat) <- uniGenes	
	##use classifier to predict again including the 3 excluded samples
##	sdat <- sdat[diffGenes.f, ]

	pred <- pamr.predict(pam.rslt, datsel, thresh, type="posterior")
	maxr <- apply(pred, 1, max)
	postR <- maxr/(1-maxr)

	sel.samples.f <- names(postR)[which(postR >= postRth)]
	
	clu.pred <- apply(pred, 1, which.max)
	clu.pred <- sort(clu.pred)	
	clu.pred.reord<-NULL
	for(cl in 1:3) {
		temp <- names(sort(postR[names(clu.pred[clu.pred==cl])]))
		clu.pred.reord <- c(clu.pred.reord, temp)
	}
	clu.pred <- clu.pred[clu.pred.reord]

	nam.ord <- names(clu.pred)

	sdat.sig <- datsel[signature, nam.ord]
	gclu.f <- hclust(as.dist(1-cor(t(sdat.sig))), method="complete")	
	
##	save(sdat.sig, pred, clu.pred, nam.ord, gclu.f, file="predAMC.RData")
	return(list(sdat.sig=sdat.sig, pred=pred, clu.pred=clu.pred, 
		nam.ord=nam.ord, gclu.f=gclu.f))
}
