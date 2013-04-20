\name{data-AMC}
%AMC
\alias{AMC_CRC_clinical}
\alias{AMC_sample_head}
%ge.CRC
\alias{ge.all}
\alias{selPbs}
%gaps
\alias{gapsmat}
\alias{gapsSE}
%dat
\alias{sdat}
%uniGenes
\alias{uniGenes}
%conClust
\alias{clus}
%silh
\alias{clus.f}
\alias{sdat.f}
\alias{silh}
%diffGenes
\alias{diffGenes}
%diffGenes.f
\alias{diffGenes.f}
%classifier
\alias{pam.rslt}
\alias{signature}
\alias{thresh}
\alias{err}
\alias{cents}
%predAMC
\alias{clu.pred}
\alias{sdat.sig}
\alias{pred}
\alias{gclu.f}
\alias{nam.ord}
%survival
\alias{surv}
\alias{survstats}
\alias{data4surv}
\title{
Data used in the package DeSousa2013
}                                                                        
\description{
See `details' for the description of each data files included in this package. 
}
\usage{
##see example for details
}
\details{
AMC: \code{AMC_CRC_clinical} and \code{AMC_sample_head} include clinical information and 
mapping information between microarray ids and sample ids, respectively. 

ge.CRC: \code{ge.all} is a numeric matrix including expression levels of all probesets 
of the 90 colon cancer samples. \code{selPbs} is a character vector of probeset ids 
that are present in any sample.

gaps: \code{gapsmat} is a numeric matrix of GAP statistics. \code{gapsSE} is a numeric matrix 
of standard errors of means of the GAP statistics. 

dat: \code{sdat} is a matrix of median centered expression values of top variable 
probesets. 

uniGenes: \code{uniGenes} is a character vector of unique gene symbols named by 
probesets. 

conClust: \code{clus} is a numeric vector of cluster labels named by colon cancer 
samples. 

silh: \code{sdat.f} is a matrix of median centered expression values of top variable 
genes for selected samples. \code{clus.f} is a numeric vector of cluster labels for 
selected cancer samples. \code{silh} is an object of class \code{\link[cluster:silhouette]{silhouette}}. 

diffGenes: \code{diffGenes} is a character vector of differential genes. 

diffGenes.f: \code{diffGenes.f} is a character vector of most predictive genes. 

classifier: \code{signature} is a character vector of signature genes. \code{pam.rslt} is 
an list of training results returned by \code{\link[pamr:pamr.train]{pamr.train}}. 
\code{thresh} is the selected shrinkage threshold. \code{err} is a matrix of cross validation 
error rates for different shrinkage thresholds. \code{cents} is a numeric matrix of 
PAM centroids for three subtypes. 

predAMC: \code{sdat.sig} is a numeric matrix of expression values of signature genes 
for the 90 samples. \code{pred} is a numeric matrix of posterior probabilities for 
samples to be classified to subtypes. \code{clu.pred} is a numeric vector of 
classification labels named by colon cancer samples. \code{nam.ord} is a character 
vector of samples ordered by their classification probabilies for visualization. 
\code{gclu.f} is the result of hierarchical clustering on the expression of signature 
genes for visualization. 

survival: \code{surv} is the result of \code{\link[survival:survfit]{survfit}} containing survival 
curves of the AMC data set. \code{survstats} is the result of \code{\link[survival:survdiff]{survdiff}}. 
\code{data4surv} is the data to perform survival analysis. 
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
data(AMC, package="DeSousa2013")
}










