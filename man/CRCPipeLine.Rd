\name{CRCPipeLine}
\alias{CRCPipeLine}
\alias{DeSousa2013}
\title{
Pipeline function for De Sousa 2013 
}                                                                        
\description{
This function reproduces the main pipeline to analyze the AMC-AJCCII-90 
microarray data set in De Sousa et al. accepted by Nature Medicine in 2013.
}
\usage{
CRCPipeLine(celpath=".", AMC_sample_head, AMC_CRC_clinical, preprocess=FALSE, 
gap.ntops = c(2, 4, 8, 12, 16, 20)*1000, gap.K.max = 6, gap.nboot = 100, 
MADth=0.5, conClust.maxK=12, conClust.reps=1000, diffG.pvalth=0.01, 
diffG.aucth=0.9, savepath=".")
}
\arguments{
\item{celpath}{
the path to a directory where .CEL files of the AMC-AJCCII-90 set are placed.
}
\item{AMC_sample_head}{
a data frame include clinical information and mapping information between 
microarray ids and sample ids (details in \code{\link[DeSousa2013:AMC_sample_head]{AMC_sample_head}}).
}
\item{AMC_CRC_clinical}{
a data frame including  include clinical information (details in \code{\link[DeSousa2013:AMC_CRC_clinical]{AMC_CRC_clinical}}).
}
\item{preprocess}{
a logical value specifying whether preprocessing of microarrays should be 
performed. If FALSE (default), preprocessed expression data \code{ge.CRC} will 
be loaded automatically.
}
\item{gap.ntops}{
an integer vector of top variable genes, measured by MAD (median absolute 
deviation).
}
\item{gap.K.max}{
an integer value specifying the maximal number of clusters to compute GAP 
statistics.
}
\item{gap.nboot}{
an integer value specifying the number of bootstraps, which is an argument \code{B} 
of function \code{\link[cluster:clusGap]{clusGap}}. 
}
\item{MADth}{
an numeric value specifying the cutoff of MAD (median absolute deviation).
}
\item{conClust.maxK}{
an integer value specifying the maximal number of clusters to perform consensus 
clustering.
}
\item{conClust.reps}{
an integer value specifying resampling times.
}
\item{diffG.pvalth}{
a numeric value specifying the fdr cutoff to select differential genes between 
subtypes by SAM.
}
\item{diffG.aucth}{
a numeric value specifying the AUC cutoff.
}
\item{savepath}{
the path to a directory where figures will be saved.
}
}
\details{
The function wraps up the preprocessing, feature selection, classification and 
subtype clinical characterization steps altogether in a signal function. It 
reproduces the main results and figures of De Sousa 2013, Nature Medicine. 
}
\value{
This function will save figure of gap statistics, Silhouette information, PAM 
cross-validation error rates, classification heatmap as well as Kaplan Meier 
plot to \code{savepath}.
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[DeSousa2013:geneExpPre]{geneExpPre}}, \code{\link[DeSousa2013:compGapStats]{compGapStats}}, \code{\link[DeSousa2013:selTopVarGenes]{selTopVarGenes}}, 
\code{\link[DeSousa2013:pbs2unigenes]{pbs2unigenes}}, \code{\link[DeSousa2013:conClust]{conClust}}, \code{\link[DeSousa2013:filterSamples]{filterSamples}}, 
\code{\link[DeSousa2013:findDiffGenes]{findDiffGenes}}, \code{\link[DeSousa2013:filterDiffGenes]{filterDiffGenes}}, \code{\link[DeSousa2013:buildClassifier]{buildClassifier}}, 
\code{\link[DeSousa2013:pamClassify]{pamClassify}}, \code{\link[DeSousa2013:progAMC]{progAMC}}
}
\examples{
\dontrun{
##This function may take a long time (hours) to finish.
##Please use R-2.15 and corresponding frma package to run the function. 
##Please contact the author if there is any confusion
data(AMC, package="DeSousa2013")
CRCPipeLine(celpath=".", AMC_sample_head, AMC_CRC_clinical, savepath=".")
}
}












