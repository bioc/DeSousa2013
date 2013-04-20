COL_GE_HIGH <- "darkorange"
COL_GE_LOW <- "#1483E4"
COL_GE_ZERO <- "white"
##cluster label color
COL_SUB_CLASSICAL <- "#98B9D5"
COL_SUB_MSI <- "#B2DF8A"
COL_SUB_SERRATED <- "#00008B"
COL_SUB_CLASSICAL_PRE <- "brown"
COL_SUB_SERRATED_PRE <- "purple"
##mutation mark
COL_MUT <- colors()[99]
##COL_MUT_BK <- "aquamarine"
COL_MUT_BK <- "white"
COL_MUT_NA <- gray(0.5)
##heatmap setting
HEATMAP_COL <- c(colorRampPalette(c(COL_GE_LOW, COL_GE_ZERO))(7)[1:6], 
	COL_GE_ZERO, 
	colorRampPalette(c(COL_GE_ZERO, COL_GE_HIGH))(7)[2:7])
HEATMAP_BKS <- c(-10, -2, -1.5, -0.8, -0.5, -0.2, -0.01, 0.01, 0.2, 0.5, 
	0.8, 1.5, 2, 10)
HEATMAP_SAMP <- c(-2.1, c(-2, -1.5, -0.8, -0.5, -0.2, -0.01, 0.01, 0.2, 0.5, 
	0.8, 1.5, 2)+0.01)
##for cell line classification
HEATMAP_COL3 <- c(colorRampPalette(c(COL_GE_LOW, "white"))(120)[1:99], 
	"white", 
	colorRampPalette(c("white", COL_GE_HIGH))(120)[22:120])
HEATMAP_BKS3 <- c(-20, -99:(-1)*0.02, 1:99*0.02, 20)
HEATMAP_SAMP3 <- c(-20, c(-9:(-1)*0.2,0, 1:9*0.2)+0.01, 20)
##association heatmap 
ASSO_HEATMAP_COL <- c(colorRampPalette(c("#8CFEEA", COL_GE_ZERO))(8)[1:7], 
	COL_GE_ZERO, 
	colorRampPalette(c(COL_GE_ZERO, "#FE8CA0"))(8)[2:8])
ASSO_HEATMAP_BKS <- c(0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.49, 0.51, 
	0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1)
ASSO_HEATMAP_SAMP <- c(0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.49, 0.51, 
	0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1)
##GAP curves
figGAP <- function(gapsmat, gapsSE) {
##	gapsmat <- gapsSE <- NULL
##	data("gaps", package="DeSousa2013", envir = environment())
	par(mfrow=c(2, 3))
	for(i in 1:ncol(gapsmat)) {
		plotCI(2:6, gapsmat[2:6, i], uiw = gapsSE[2:6, i], 
			liw = gapsSE[2:6, i], xlab="No. of clusters", ylab="GAP")
		lines(2:6, gapsmat[2:6, i])
		title(main=paste("top ", colnames(gapsmat)[i], sep=""))
	}
}

##Silhouette information
figSilh <- function(silh) {
##	clus.f <- sdat.f <- silh <- NULL
##	data("silh", package="DeSousa2013", envir = environment())
	plot(silh, main="Silhouette Information", col=c(COL_SUB_CLASSICAL, 
		COL_SUB_MSI, COL_SUB_SERRATED))
}

##PAM cross validation
figPAMCV <- function(err) {
##	pam.rslt <- signature <- thresh <- err <- cents <- NULL
##	data("classifier", package="DeSousa2013", envir = environment())
	boxplot(err, xlab="No. of genes", ylab="error rate")
}

##Classification of AMC-AJCCII-90 set
figClassify <- function(AMC_CRC_clinical, pred, clu.pred, sdat.sig, gclu.f, nam.ord) {
##	AMC_CRC_clinical <- AMC_sample_head <- pam.rslt <- signature <- 
##		thresh <- err <- cents <- clu.pred <- sdat.sig <- pred <- 
##		gclu.f <- nam.ord <- NULL
##	data("AMC", package="DeSousa2013", envir = environment())
##	data("classifier", package="DeSousa2013", envir = environment())
##	data("predAMC", package="DeSousa2013", envir = environment())

	layout(matrix(c(1:4, 5, 5, 5, 0), 8, 1), 
		heights=c(1.8, 20, rep(0.8, 5), 1))
	##1.1 cluster labels
	par(mar=c(0.1, 8, 2.1, 0.5))
	image(matrix(clu.pred, length(clu.pred),1), col=c(COL_SUB_CLASSICAL, 
		COL_SUB_MSI, COL_SUB_SERRATED), axes=FALSE)	
	axis(side=3, at=c(sum(clu.pred==1)/2, sum(clu.pred==1)+sum(clu.pred==2)/2, 
		sum(clu.pred==1)+sum(clu.pred==2)+sum(clu.pred==3)/2)/length(clu.pred), 
		labels=c("CCS1", "CCS2", "CCS3"), line=-1, tick=FALSE, cex.axis=1.8)
	box(lwd=0.5, col=gray(0.5))
	##1.2 heatmap
	par(mar=c(0.5, 8, 0.1, 0.5))
	image(t(sdat.sig[rev(gclu.f$order), ]), col=HEATMAP_COL3, 
		breaks=HEATMAP_BKS3, axes=FALSE)	
	mtext("classifier", side=2, cex=1.2, line=2)
	box(lwd=0.5, col=gray(0.5))
	##1.3 predicted posterior probs
	par(mar=c(0.1, 8, 0.1, 0.5))
	barplotmat <- t(pred[nam.ord, ])
	colnames(barplotmat) <- NULL
	barplot(barplotmat, beside=FALSE, axes=FALSE, xaxs="i", yaxs="i", xlab="", 
		col=c(COL_SUB_CLASSICAL, COL_SUB_MSI, COL_SUB_SERRATED), border=NA)
	mtext("Prob.", side = 2, line =  0.2, las=2, cex=0.8)
	box(lwd=0.5, col=gray(0.5))
	##1.4 relapse info
	par(mar=c(0.1, 8, 0.1, 0.5))
	Recur.mat <- rep(1, ncol(sdat.sig))
	names(Recur.mat) <- colnames(sdat.sig)
	Recur.mat[intersect(names(Recur.mat), 
		rownames(AMC_CRC_clinical)[AMC_CRC_clinical[, "Met"]=="yes"])] <- 2
	image(matrix(Recur.mat, length(Recur.mat), 1), 
		col=c(COL_MUT_BK, COL_MUT), axes=FALSE)
	mtext("Relapse", side = 2, line = 0.2, las=2, cex=0.8)
	grid(nx=length(clu.pred), ny=0, col=gray(0.8), lty="solid", lwd=0.5)
	box(lwd=0.5, col=gray(0.5))
	##color legend for the heatmap
	par(mar=c(1.5, 45, 2, 0.5))
	image(t(matrix(HEATMAP_SAMP3, nrow=1)), col=HEATMAP_COL3, 
		breaks=HEATMAP_BKS3, axes=FALSE)	
	axis(side=1, at=c(0.02, 0.98), labels=c('low', 'high'), tick=FALSE, 
		cex=1, cex.axis=1, line=-1)
	box(lwd=0.5, col=gray(0.5))
}
##KM plot for AMC subtype prognosis
figKM <- function(surv, survstats) {

##	event <- as.character(AMC_CRC_clinical[, "Met"])
##	time <- AMC_CRC_clinical[, "timeMETRec"]
##	names(event) <- names(time) <- rownames(AMC_CRC_clinical)
##	x <- rep("CCS3", length(clu.pred))
##	x[clu.pred==1] <- "CCS1"
##	x[clu.pred==2] <- "CCS2"
##	status <- ifelse(event=="yes", 1, 0)
##	data4surv <- data.frame(time=time[names(clu.pred)], 
##		status=status[names(clu.pred)], x=x)
##	surv <- survfit(Surv(time, status) ~ x, data = data4surv)	
##	survstats <- survdiff(Surv(time, status) ~ x, data = data4surv)
##	survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) - 1)

	plot(surv, col=c(COL_SUB_CLASSICAL, COL_SUB_MSI, COL_SUB_SERRATED), 
		xlab="Follow up (months)", ylab="DFS (prob.)", lwd=3, 
		ylim=c(0.2, 1), axes=FALSE)
	axis(side=1, at=1:6*20*30, labels=1:6*20)
	axis(side=2, at=1:5*0.2, labels=1:5*0.2)
	box()
	text(x=800, y=0.25, labels=paste("p-value: ", 
		format(survstats$p.value, scientific=TRUE, digits=3), sep=""), cex=1)
	legend(legend = c("CCS1", "CCS2","CCS3"), lty=c(1,1,1),lwd=c(3,3,3), 
		col=c(COL_SUB_CLASSICAL, COL_SUB_MSI, COL_SUB_SERRATED), bty='n', cex=1,
		x=2800, y=0.4)		
}













