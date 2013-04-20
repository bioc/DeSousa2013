\name{compGapStats}
\alias{compGapStats}
\alias{figGAP}
\title{
Computing Gap statistics to identify the optimal number of subtypes
}                                                                        
\description{
Compute Gap statistics to identify the optimal number of subtypes
}
\usage{
compGapStats(ge.CRC, ntops=c(2, 4, 8, 12, 16, 20)*1000, K.max=6, nboot=100)
figGAP(gapsmat, gapsSE)
}
\arguments{
\item{ge.CRC}{
a numeric matrix of expression data of genes expressed in at least one sample. 
}
\item{ntops}{
an integer vector of top variable genes, measured by MAD (median absolute 
deviation).
}
\item{K.max}{
an integer value specifying the maximal number of clusters to compute GAP 
statistics.
}
\item{nboot}{
an integer value specifying the number of bootstraps, which is an argument \code{B} 
of function \code{\link[cluster:clusGap]{clusGap}}. 
}
\item{gapsmat}{
a numeric matrix of GAP statistics.
}
\item{gapsSE}{
standard errors of means of the GAP statistics.
}
}
\details{
GAP statistic is a popular method to estimate the number of clusters in a set 
of data by comparing the change in observed and expected within-cluster 
dispersion. To identify the optimal number of clusters, GAP statistic can be 
computed for k=1 to \code{K.max} with \code{nboot} bootstraps for \code{ntops} 
top variable genes in the AMC data set. 

The function \code{figGAP} is designed to visualize GAP curves.
}
\value{
This function will return a list including \code{gapsmat} (a numeric matrix of 
GAP statistics) and \code{gapsSE} (standard errors of means of the GAP 
statistics).
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

Tibshirani, Robert andWalther, Guenther and Hastie, Trevor (2001). Estimating
the number of clusters in a data set via the gap statistic. Journal of the 
Royal Statistical Society: Series B (Statistical Methodology), 63(2), 411-423.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\examples{
data(ge.CRC, package="DeSousa2013")
ge.CRC <- ge.all[selPbs, ]
gaps <- compGapStats(ge.CRC, ntops=c(2, 4)*1000, K.max=6, nboot=10)
figGAP(gaps$gapsmat, gaps$gapsSE)
}












