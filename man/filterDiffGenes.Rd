\name{filterDiffGenes}
\alias{filterDiffGenes}
\title{
Select the most predictive genes from differential genes
}                                                                        
\description{
This function evaluates AUC (area under ROC curve) to select the most predictive 
genes from differential genes.
}
\usage{
filterDiffGenes(sdat.f, clus.f, diffGenes, aucth=0.9)
}
\arguments{
\item{sdat.f}{
a matrix of median centered expression values of top variable genes for selected 
samples.
}
\item{clus.f}{
a numeric vector of cluster labels for selected cancer samples.
}
\item{diffGenes}{
a character vector of differential genes.
}
\item{aucth}{
a numeric value specifying the AUC cutoff.
}
}
\details{
After obtaining differential genes between subtypes, we calculated AUC 
(area under ROC curve, using package \code{ROCR}) to assess each gene's ability 
to separate one subtype from the others. 
}
\value{
This function will return \code{diffGenes.f}, which is a character vector of 
most predictive genes.
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

Sing, T., Sander, O., Beerenwinkel, N., Lengauer, T. (2005). ROCR: visualizing 
classifier performance in R. Bioinformatics 21(20):3940-3941.  
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[DeSousa2013:findDiffGenes]{findDiffGenes}}
}
\examples{
data(uniGenes)
data(silh)
data(diffGenes)
##select randomly part of the whole differential genes only for testing the function
diffGenes <- diffGenes[sample(1:length(diffGenes), 100)]
diffGenes.f <- filterDiffGenes(sdat.f, clus.f, diffGenes, aucth=0.9)
}












