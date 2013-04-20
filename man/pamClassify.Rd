\name{pamClassify}
\alias{pamClassify}
\alias{figClassify}
\title{
Classifying the AMC samples
}                                                                        
\description{
These functions classify and visualize the 90 AMC colon cancer samples.
}
\usage{
pamClassify(datsel, signature, pam.rslt, thresh, postRth=1)
figClassify(AMC_CRC_clinical, pred, clu.pred, sdat.sig, gclu.f, nam.ord)
}
\arguments{
\item{datsel}{
a numeric value specifying the cutoff of posterior odds in favor of optimal 
subtype to the other two. 
}
\item{signature}{
a character vector of signature genes. 
}
\item{pam.rslt}{
an list of training results returned by \code{\link[pamr:pamr.train]{pamr.train}}. 
}
\item{thresh}{
the selected shrinkage threshold. 
}
\item{postRth}{
a numeric value specifying the cutoff of posterior odds in favor of optimal 
subtype to the other two. 
}
\item{AMC_CRC_clinical}{
a data frame including  include clinical information (details in \code{\link[DeSousa2013:AMC_CRC_clinical]{AMC_CRC_clinical}}).
}
\item{pred}{
a numeric matrix of posterior probabilities for samples to be classified to 
subtypes.
}
\item{clu.pred}{
a numeric vector of classification labels named by colon cancer samples.
}
\item{sdat.sig}{
a numeric matrix of expression values of signature genes for the 90 samples.
}
\item{gclu.f}{
the result of hierarchical clustering on the expression of signature genes for 
visualization.
}
\item{nam.ord}{
a character vector of samples ordered by their classification probabilies for 
visualization.
}
}
\details{
The function \code{pamClassify} classifies the 90 AMC colon cancer samples 
using the gene expression based classifier built by PAM.
The function \code{figClassify} generates classification results with a heatmap 
of median centred expression values across samples, a track indicating 
classification probabilities, as well as a track indicating relapse of patients.
}
\value{
This function \code{pamClassify} will return a list including \code{sdat.sig} 
(a numeric matrix of expression values of signature genes for the 90 samples), 
\code{pred} (a numeric matrix of posterior probabilities for samples to be 
classified to subtypes), \code{nam.ord} (a character vector of samples ordered 
by their classification probabilies for visualization), \code{gclu.f} (result 
of hierarchical clustering on the expression of signature genes for 
visualization).
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
\code{\link[DeSousa2013:buildClassifier]{buildClassifier}}
}
\examples{
data(AMC, package="DeSousa2013")
data(dat, package="DeSousa2013")
data(uniGenes, package="DeSousa2013")
data(diffGenes.f, package="DeSousa2013")
data(classifier, package="DeSousa2013")
datsel <- sdat[names(uniGenes), ]
rownames(datsel) <- uniGenes	
datsel <- datsel[diffGenes.f, ]
pamcl <- pamClassify(datsel, signature, pam.rslt, thresh, postRth=1)
figClassify(AMC_CRC_clinical, pamcl$pred, pamcl$clu.pred, pamcl$sdat.sig, 
pamcl$gclu.f, pamcl$nam.ord)
}












