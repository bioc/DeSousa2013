\name{filterSamples}
\alias{filterSamples}
\alias{figSilh}
\title{
Filter colon cancer samples by Silhouette width
}                                                                        
\description{
This function computes Silhouette widths for the 90 colon cancer samples.
}
\usage{
filterSamples(sdat, uniGenes, clus)
figSilh(silh)
}
\arguments{
\item{sdat}{
a matrix of median centered expression values of top variable probesets.
}
\item{uniGenes}{
a character vector of unique gene symbols named by probesets.
}
\item{clus}{
a numeric vector of cluster labels named by colon cancer samples.
}
\item{silh}{
an object of class \code{\link[cluster:silhouette]{silhouette}}.
}
}
\details{
Silhouette width was computed to identify the most representative samples 
within each cluster. Samples with positive silhouette width were retained 
to build the PAM classifier.
}
\value{
This function will return a list including \code{sdat.f} (a matrix of median 
centered expression values of top variable genes for selected samples), 
\code{clus.f} (a numeric vector of cluster labels for selected cancer samples) 
and \code{silh} (an object of class \code{\link[cluster:silhouette]{silhouette}}).
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

Rousseeuw, Peter J (1987). Silhouettes: a graphical aid to the interpretation
and validation of cluster analysis Journal of computational and applied 
mathematics, 20, 53-65.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
\code{\link[DeSousa2013:selTopVarGenes]{selTopVarGenes}}
}
\examples{
data(uniGenes, package="DeSousa2013")
data(dat, package="DeSousa2013")
data(conClust, package="DeSousa2013")
samp.f <- filterSamples(sdat, uniGenes, clus)
figSilh(samp.f$silh)
}












