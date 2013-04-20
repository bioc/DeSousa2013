\name{buildClassifier}
\alias{buildClassifier}
\alias{figPAMCV}
\title{
Build a gene expression based classifier
}                                                                        
\description{
This function employs PAM to build a gene expression based classifier.
}
\usage{
buildClassifier(sigMat, clus.f, nfold=10, nboot=100)
figPAMCV(err)
}
\arguments{
\item{sigMat}{
a matrix of median centered expression values of selected genes (by function 
\code{\link[DeSousa2013:filterDiffGenes]{filterDiffGenes}}) for selected cancer samples (by function \code{\link[DeSousa2013:filterSamples]{filterSamples}}).
}
\item{clus.f}{
a numeric vector of cluster labels for selected cancer samples.
}
\item{nfold}{
an integer value specifying the fold of cross validation.
}
\item{nboot}{
an integer value specifying the number of bootstraps.
}
\item{err}{
a matrix of cross validation error rates for different shrinkage thresholds.
}
}
\details{
The expression data of the retained most predictive genes were trained 
by PAM to build a robust classifier. To select the optimal threshold for 
centroid shrinkage, we performed 10-fold cross-validation over a range of 
shirnkage thresholds for 1000 iterations, and selected the one yielding a 
good performance (error rate < 2\%) with the least number of genes. 
}
\value{
This function will return a list including \code{signature} (a character vector 
of signature genes), \code{pam.rslt} (an list of training results returned by 
\code{\link[pamr:pamr.train]{pamr.train}}), \code{thresh} (the selected shrinkage 
threshold), \code{err} (a matrix of cross validation error rates for different 
shrinkage thresholds), \code{cents} (a numeric matrix of PAM centroids for three 
subtypes).
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

Tibshirani, Robert and Hastie, Trevor and Narasimhan, Balasubramanian and Chu, 
Gilbert (2002). Diagnosis of multiple cancer types by shrunken centroids of 
gene expression. PNAS, 99(10), 6567-6572.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[DeSousa2013:filterSamples]{filterSamples}}, \code{\link[DeSousa2013:findDiffGenes]{findDiffGenes}}, \code{\link[DeSousa2013:filterDiffGenes]{filterDiffGenes}}
}
\examples{
data(diffGenes.f, package="DeSousa2013")
data(silh, package="DeSousa2013")
sigMat <- sdat.f[diffGenes.f, names(clus.f)]
##set a small number of bootstraps, only for testing the function
classifier <- buildClassifier(sigMat, clus.f, nboot=10)
data(classifier, package="DeSousa2013")
figPAMCV(err)
}












