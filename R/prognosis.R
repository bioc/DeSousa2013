################################################################################
########################      Clinical data analysis     #######################
################################################################################
##Survival analysis
progAMC <- function(AMC_CRC_clinical, AMC_sample_head, clu.pred) {
##	AMC_CRC_clinical <- AMC_sample_head <- sdat.sig <- pred <- 
##		clu.pred <- nam.ord <- gclu.f <- NULL 
##	data("AMC", package="DeSousa2013", envir = environment())
##	data("predAMC", package="DeSousa2013", envir = environment())
	
	event <- as.character(AMC_CRC_clinical[, "Met"])
	time <- AMC_CRC_clinical[, "timeMETRec"]
	names(event) <- names(time) <- rownames(AMC_CRC_clinical)
	
	x <- rep(3, length(clu.pred))
	x[clu.pred==1] <- 1
	x[clu.pred==2] <- 2
	status <- ifelse(event=="yes", 1, 0)
	data4surv <- data.frame(time=time[names(clu.pred)], status=status[names(clu.pred)], 
		x=x)
	surv <- survfit(Surv(time, status) ~ x, data = data4surv)
	survstats <- survdiff(Surv(time, status) ~ x, data = data4surv)
	survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) - 1)
	
	print(survstats)
##	save(surv, survstats, data4surv, file="survival.RData")
	return(list(surv=surv, survstats=survstats, data4surv=data4surv))
}
