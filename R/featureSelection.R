################################################################################
#############            Feature Selection/classification          #############
################################################################################
##filtering samples by silhouette information
filterSamples <- function(sdat, uniGenes, clus) {
##	uniGenes <- sdat <- clus <- NULL
##	data("uniGenes", package="DeSousa2013", envir = environment())
##	data("dat", package="DeSousa2013", envir = environment())
##	data("conClust", package="DeSousa2013", envir = environment())

	sdat.ord <- sdat[, names(clus)]
	silh <- silhouette(clus, as.dist(1-cor(sdat.ord)))
	excl <- names(clus)[which(silh[, "sil_width"] < 0)]
	sdat.f <- sdat.ord[, setdiff(colnames(sdat.ord), excl)]
	clus.f <- clus[which(names(clus) %in% colnames(sdat.f))]
	clus.f <- sort(clus.f)

	sdat.f <- sdat.f[names(uniGenes), ]
	rownames(sdat.f) <- uniGenes
##	save(sdat.f, clus.f, silh, file="silh.RData")
	return(list(sdat.f=sdat.f, clus.f=clus.f, silh=silh))
}
##Search for differential genes by SAM
findDiffGenes <- function(sdat.f, clus.f, pvalth=0.01) {
##	uniGenes <- clus <- clus.f <- sdat.f <- silh <- NULL
##	data("uniGenes", package="DeSousa2013", envir = environment())
##	data("conClust", package="DeSousa2013", envir = environment())
##	data("silh", package="DeSousa2013", envir = environment())
	diffGenes <- NULL
	for(cl in 1:3) {
		temp <- clus.f
		temp[clus.f==cl] <- 0
		temp[clus.f!=cl] <- 1
		res.sam <- sam(sdat.f, temp, gene.names=rownames(sdat.f), rand=123456)
		temp.diffGenes <- names(res.sam@q.value)[res.sam@q.value <= pvalth]
		diffGenes <- c(diffGenes, temp.diffGenes)
	}
	diffGenes <- unique(diffGenes)
##	sigdat <- sdat.f[diffGenes, ]
##	save(diffGenes, file="diffGenes.RData")
	return(diffGenes)
}
##select most predictive genes
filterDiffGenes <- function(sdat.f, clus.f, diffGenes, aucth=0.9) {
##	uniGenes <- clus.f <- sdat.f <- silh <- diffGenes <- NULL
##	data("uniGenes", package="DeSousa2013", envir = environment())
##	data("silh", package="DeSousa2013", envir = environment())
##	data("diffGenes", package="DeSousa2013", envir = environment())
	sigdat <- sdat.f[diffGenes, ]
	sigGeneIds <- rownames(sigdat)
	compareTwoClus <- function(sigdat, clus.f, clu=1) {
		
		temp <- clus.f
		temp[clus.f==clu] <- 0
		temp[clus.f!=clu] <- 1

		temp2 <- temp
		temp2[temp==0] <- 1
		temp2[temp==1] <- 0
		
		sigdat <- sigdat[, names(temp)]

		preds1 <- apply(sigdat, 1, prediction, labels = temp)
		preds2 <- apply(sigdat, 1, prediction, labels = temp2)

		aucs1 <- unlist(lapply(preds1, function(x) {
			(performance(x, measure = "auc"))@y.values[[1]]
		}))
		aucs2 <- unlist(lapply(preds2, function(x) {
			(performance(x, measure = "auc"))@y.values[[1]]
		}))
		aucs <- rowMax(cbind(aucs1, aucs2))
		
		diffGenes.f <- rownames(sigdat)[which(aucs >= aucth)]	
		return(diffGenes.f)
	}	
	sigs1 <- compareTwoClus(sigdat, clus.f,1)
	sigs2 <- compareTwoClus(sigdat, clus.f,2)
	sigs3 <- compareTwoClus(sigdat, clus.f,3)
	diffGenes.f <- unique(c(sigs1, sigs2, sigs3))
	return(diffGenes.f)
##	save(diffGenes.f, file="diffGenes.f.RData")
}

