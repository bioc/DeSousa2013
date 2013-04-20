\name{conClust}
\alias{conClust}
\title{
Consensus clustering
}                                                                        
\description{
This function performs consensus clustering for the AMC-AJCCII-90 data set.
}
\usage{
conClust(sdat, maxK=12, reps=1000, savepath = ".")
}
\arguments{
\item{sdat}{
a matrix of median centered expression values of top variable probesets.
}
\item{maxK}{
an integer value specifying the maximal number of clusters to perform consensus 
clustering.
}
\item{reps}{
an integer value specifying resampling times.
}
\item{savepath}{
the path to a directory where figures will be saved.
}
}
\details{
Using the most variable probesets (MAD>0.5), we performed hierarchical 
clustering with agglomerative average linkage to cluster these samples. 
Consensus clustering was employed, with 1000 iterations and 0.98 subsampling 
ratio, to assess the clustering stability. 
}
\value{
This function will return \code{clus}, which is a numeric vector of cluster 
labels named by colon cancer samples.
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

Monti, Stefano and Tamayo, Pablo and Mesirov, Jill and Golub, Todd (2003). 
Consensus clustering: a resampling-based method for class discovery and 
visualization of gene expression microarray data. Machine learning, 52(1), 
91-118.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\examples{
library(ConsensusClusterPlus)
data(dat, package="DeSousa2013")
##set a small number of bootstraps only for testing the function
clus <- conClust(sdat, maxK=6, reps=10, savepath=".")
}












