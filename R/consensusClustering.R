##4. consensus clustering
conClust <- function(sdat, maxK=12, reps=1000, savepath=".") {
##	sdat <- NULL
##	data("dat", package="DeSousa2013", envir = environment())
	res <- ConsensusClusterPlus(sdat, maxK=maxK, reps=reps, pItem=0.98, 
		innerLinkage="average", finalLinkage="average", pFeature=1, 
		title=savepath, clusterAlg="hc", distance="pearson", plot="pdf")
	areaK <- res[["areaK"]]
	res[["areaK"]] <- NULL
	icl <- calcICL(res, title=savepath, plot="pdf")
	
	con <- icl[[2]][icl[[2]][, "k"]==3, ]
	conMat <- cbind(con[con[, 2]==1, 4], con[con[, 2]==2, 4], con[con[, 2]==3, 4])
	rownames(conMat) <- con[con[, 2]==1, 3]
	clus <- apply(conMat, 1, which.max)
##	lclus <- sapply(1:3, function(x) rownames(conMat)[which(clus==x)])
##	save(clus, file="conClust.RData")
	return(clus)
}
