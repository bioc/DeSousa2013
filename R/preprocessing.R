################################################################################
#########################          Preprocessing        ########################
################################################################################
##1 Microarray data preprocessing
##	This program employs fRMA and barcode to preprocess the AMC colon cancer 
##	dataset. The arrays came from two batches: CRC+normals, precursors
geneExpPre <- function(celpath, AMC_sample_head) {
##	AMC_CRC_clinical <- AMC_sample_head <- NULL
##	data("AMC", package="DeSousa2013", envir = environment())
	AllIDs <- as.character(AMC_sample_head[, 'array'])
	############################################################################
	##(1) Use fRMA to preprocess CRC+norm and adenoma samples separately
	############################################################################
	##batch 1: CRC samples + normal samples
	CRCIDs <- AllIDs[grep('GSM', AllIDs)]
	normIDs <- AllIDs[grep('COL', AllIDs)]
	CRCfiles <- file.path(celpath, paste(CRCIDs, 'CEL', sep='.'))
	normfiles <- file.path(celpath, paste(normIDs, 'CEL', sep='.'))
	
	batch1 <- ReadAffy(filenames=c(CRCfiles, normfiles))
	eset1 <- frma(batch1, summarize="random_effect")
	ge.batch1 <- exprs(eset1)
	colnames(ge.batch1)<-as.character(AMC_sample_head[match(colnames(ge.batch1), 
		paste(AMC_sample_head[, "array"], "CEL", sep=".")), "sample"])	
		
	##batch 2: adenoma samples
	adenomaIDs <- AllIDs[grep('EXT', AllIDs)]	
	adenomafiles <- file.path(celpath, paste(adenomaIDs, 'CEL', sep='.'))
	batch2 <- ReadAffy(filenames=adenomafiles)
	eset2 <- frma(batch2, summarize="random_effect")
	ge.batch2 <- exprs(eset2)
	colnames(ge.batch2)<-as.character(AMC_sample_head[match(colnames(ge.batch2), 
		paste(AMC_sample_head[, "array"], "CEL", sep=".")), "sample"])	
			
	##check batch effects
##	pc <- prcomp(t(cbind(ge.batch1, ge.batch2)), scale=TRUE)
##	cols <- rep('black', nrow(pc$x))
##	cols[grep('ext', rownames(pc$x))] <- 'red'
##	plot3d(pc$x[, 1], pc$x[, 2], pc$x[, 3], col=cols, size=2, type='s')
	
	##! indeed there is batch effect
	############################################################################
	##(2) Use ComBat to correct for batch effect
	############################################################################	

	batches <- c(rep(1, ncol(ge.batch1)), rep(2, ncol(ge.batch2)))
	ge.corr <- ComBat(dat=cbind(ge.batch1, ge.batch2), batch=batches, mod=NULL)
	
	##check batch effects again
##	pc <- prcomp(t(ge.corr), scale=TRUE)
##	cols <- rep('black', nrow(pc$x))
##	cols[grep('ext', rownames(pc$x))] <- 'red'
##	plot3d(pc$x[, 1], pc$x[, 2], pc$x[, 3], col=cols, size=2, type='s')	
	
	##! two batches are now well mixed
	
	############################################################################
	##(3) Use fRMA to detect expressed genes
	############################################################################			
	apc1 <- barcode(eset1, cutoff=(-log10(0.05)), output="binary")
	apc2 <- barcode(eset2, cutoff=(-log10(0.05)), output="binary")	
	##find expressed genes in both batches
	selPbs <- rownames(apc1)[which(rowSums(cbind(apc1, apc2))>1)]
	
	ge.all <- ge.corr[, as.character(AMC_sample_head[match(CRCIDs, 
		AMC_sample_head[, 1]),2])]					##CRC samples, all genes
	
	return(list(ge.all=ge.all, selPbs=selPbs))
##	save(ge.all, selPbs, file=savefile)	
}
##2 compute Gap statistic
##	This program checks the stability of getting 3 clusters as optimal with 
##	different numbers of top variable genes
compGapStats <- function(ge.CRC, ntops=c(2, 4, 8, 12, 16, 20)*1000, K.max=6, 
	nboot=100) {
	
##	ge.all <- selPbs <- NULL
##	data("ge.CRC", package="DeSousa2013", envir = environment())
##	ge.CRC <- ge.all[selPbs, ]
	MAD <- apply(ge.CRC, 1, mad)
	ords <- names(sort(MAD, decreasing=TRUE))

	fun <- function(x, k) {
		y <- hclust(as.dist(1-cor(t(x))), method="average")
		clus <- cutree(y, k=k);
		return(list(cluster=clus))
	}
	gaps <- list()	
	for(i in 1:length(ntops)) {
		ntop <- ntops[i]
		sdat <- ge.CRC[ords[1:ntop], ]
		sdat = sweep(sdat,1, apply(sdat,1,median))	
		gaps[[i]] <- clusGap(t(sdat), FUNcluster=fun, K.max=K.max, B=nboot)		
	}
	gapsmat <- matrix(0, K.max, length(ntops))
	gapsSE <- matrix(0, K.max, length(ntops))
	for(i in 1:length(ntops)) {
		gapsmat[, i] <- gaps[[i]]$Tab[, 3]
		gapsSE[, i] <- gaps[[i]]$Tab[, 4]
	}
	
	colnames(gapsmat) <- colnames(gapsSE) <- ntops
	rownames(gapsmat) <- rownames(gapsSE) <- 1:K.max
	##! 14 out of 20 favor three clusters	
	return(list(gapsmat=gapsmat, gapsSE=gapsSE))
##	save(gapsmat, gapsSE, file=savefile)
}

##3	Select genes of top variability
##! The GAP statistic confirms that 3 is indeed the optimal number of clusters
##	regardless of the number of variable genes selected. From this step, we 
##	use genes with MAD > 0.5 to do consensus clustering and GAP statistic 
selTopVarGenes <- function(ge.CRC, MADth=0.5) {
##	load("ge.CRC.expr.RData")
##	ge.all <- selPbs <- NULL
##	data("ge.CRC", package="DeSousa2013", envir = environment())
##	ge.CRC <- ge.all[selPbs, ]
	MAD <- apply(ge.CRC, 1, mad)
	sdat <- ge.CRC[which(MAD > MADth), ]
	sdat = sweep(sdat,1, apply(sdat,1,median))
##	save(sdat, file=savefile)
	return(sdat)
}

##	Transform probesets to unique genes
pbs2unigenes <- function(ge.CRC, sdat) {
##	sdat <- ge.all <- selPbs <- NULL
##	data("ge.CRC", package="DeSousa2013", envir = environment())
##	data("dat", package="DeSousa2013", envir = environment())
##	ge.CRC <- ge.all[selPbs, ]
	x <- hgu133plus2SYMBOL
	mapped_probes <- mappedkeys(x)
	xx <- AnnotationDbi::as.list(x[mapped_probes])
	
	pbs <- rownames(sdat)
	lsyms <- xx[pbs]
	syms <- unlist(lsyms)
	uni.syms <- unique(syms)
	ids <- NULL
	for(i in 1:length(uni.syms)) {
		inds <- which(syms==uni.syms[i])
		tempids <- names(syms)[inds]
		if(length(tempids)>1) {
			ids <- c(ids, tempids[which.max(rowMeans(ge.CRC[tempids, ]))])
		} else 
			ids <- c(ids, tempids)
	}	
	uniGenes <- uni.syms
	names(uniGenes) <- ids
##	save(uniGenes, file=savefile)
	return(uniGenes)
}
