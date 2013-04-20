CRCPipeLine <- function(
	##1.a. for microarray preprocessing
	celpath=".", 
	AMC_sample_head,
	AMC_CRC_clinical,
	preprocess=FALSE, 
	##1.b. for gap stats
	gap.ntops = c(2, 4, 8, 12, 16, 20)*1000, 
	gap.K.max = 6, 
	gap.nboot = 100, 
	##1.c. for selecting top variable genes
	MADth=0.5, 
	##2. for consensus clustering
	conClust.maxK=12, 
	conClust.reps=1000, 
	##3.b. for searching for differential genes
	diffG.pvalth=0.01, 
	##3.c. select most predictive genes
	diffG.aucth=0.9, 
	savepath="."
	) {
	if(preprocess) {
		##1.a. preprocessing microarrays
		ge.pre <- geneExpPre(celpath, AMC_sample_head)
		ge.all <- ge.pre$ge.all
		selPbs <- ge.pre$selPbs
	} else {
		ge.all <- selPbs <- NULL
		data("ge.CRC", package="DeSousa2013", envir = environment())	
	}
	ge.CRC <- ge.all[selPbs, ]
	##1.b. comput gap statistics
	gaps <- compGapStats(ge.CRC, ntops=gap.ntops, K.max=gap.K.max, nboot=gap.nboot)
	gapsmat <- gaps$gapsmat
	gapsSE <- gaps$gapsSE
	##1.c. select genes of top variability
	sdat <- selTopVarGenes(ge.CRC, MADth=MADth)
	##1.d. transform probesets to unique genes
	uniGenes <- pbs2unigenes(ge.CRC, sdat)
	
	##2. consensus clustering
	clus <- conClust(sdat, maxK=conClust.maxK, reps=conClust.reps)
	
	##3.a. filtering samples by silhouette information
	samp.f <- filterSamples(sdat, uniGenes, clus)
	sdat.f <- samp.f$sdat.f
	clus.f <- samp.f$clus.f
	silh <- samp.f$silh
	##3.b. Search for differential genes by SAM
	diffGenes <- findDiffGenes(sdat.f, clus.f, pvalth=diffG.pvalth)
	##3.c. select most predictive genes
	diffGenes.f <- filterDiffGenes(sdat.f, clus.f, diffGenes, aucth=diffG.aucth)
	##3.d. build classifier using filtered samples and genes
	sigMat <- sdat.f[diffGenes.f, names(clus.f)]
	classifier <- buildClassifier(sigMat, clus.f)
	signature <- classifier$signature
	pam.rslt <- classifier$pam.rslt
	thresh <- classifier$thresh
	err <- classifier$err
	##3.e classification
	datsel <- sdat[names(uniGenes), ]
	rownames(datsel) <- uniGenes	
	datsel <- datsel[diffGenes.f, ]
	pamcl <- pamClassify(datsel, signature, pam.rslt, thresh, postRth=1)
	sdat.sig <- pamcl$sdat.sig
	pred <- pamcl$pred
	clu.pred <- pamcl$clu.pred
	nam.ord <- pamcl$nam.ord
	gclu.f <- pamcl$gclu.f
	
	##4. prognosis
	prog <- progAMC(AMC_CRC_clinical, AMC_sample_head, clu.pred)
	surv <- prog$surv
	survstats <- prog$survstats
	data4surv <- prog$data4surv

	##reproduce figures
	setwd(savepath)
	##
	pdf(file=file.path(savepath, "gap.pdf"), width=8, height=8)
	figGAP(gapsmat, gapsSE)
	graphics.off()
	##
	pdf(file=file.path(savepath, "silh.pdf"), width=8, height=8)
	figSilh(silh)
	graphics.off()
	##
	pdf(file=file.path(savepath, "PAMcv.pdf"), width=10, height=7)
	figPAMCV(err)
	graphics.off()
	##
	pdf(file=file.path(savepath, "classification.pdf"), width=10, height=10)
	figClassify(AMC_CRC_clinical, pred, clu.pred, sdat.sig, gclu.f, nam.ord)
	graphics.off()
	##
	pdf(file=file.path(savepath, "KM.pdf"), width=7, height=7)
	figKM(surv, survstats)
	graphics.off()
}






















