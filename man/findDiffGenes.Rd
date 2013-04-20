\name{findDiffGenes}
\alias{findDiffGenes}
\title{
Search for differential genes between subtypes 
}                                                                        
\description{
This function employs SAM to search for differential genes between subtypes
}
\usage{
findDiffGenes(sdat.f, clus.f, pvalth=0.01)
}
\arguments{
\item{sdat.f}{
a matrix of median centered expression values of top variable genes for selected 
samples.
}
\item{clus.f}{
a numeric vector of cluster labels for selected cancer samples.
}
\item{pvalth}{
a numeric value specifying the fdr cutoff to select differential genes between 
subtypes by SAM.
}
}
\details{
In this function, we find differential genes between each two of three subtypes 
and take the unique genes of combined differential genes altogether.
}
\value{
This function will return \code{diffGenes}, which is a character vector of 
differential genes.
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

Tusher, Virginia Goss and Tibshirani, Robert and Chu, Gilbert (2001). 
Significance analysis of microarrays applied to the ionizing radiation response.
PNAS, 98(9), 5116-5121.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[DeSousa2013:filterDiffGenes]{filterDiffGenes}}
}
\examples{
data(uniGenes, package="DeSousa2013")
data(conClust, package="DeSousa2013")
data(silh, package="DeSousa2013")
##select randomly part of the whole data set only for testing the function
diffGenes <- findDiffGenes(sdat.f[sample(1:nrow(sdat.f), 500), ], clus.f, pvalth=0.01)
}












