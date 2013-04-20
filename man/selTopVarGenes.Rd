\name{selTopVarGenes}
\alias{selTopVarGenes}
\title{
Select top variable probesets
}                                                                        
\description{
This function selects top variable probesets across the 90 cancer samples by 
median absolute deviations. 
}
\usage{
selTopVarGenes(ge.CRC, MADth=0.5)
}
\arguments{
\item{ge.CRC}{
a numeric matrix of expression data of genes expressed in at least one sample. 
}
\item{MADth}{
an numeric value specifying the cutoff of MAD (median absolute deviation).
}
}
\value{
This function will return \code{sdat}, which is a matrix of median centered 
expression values of top variable probesets.
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
data(ge.CRC, package="DeSousa2013")
ge.CRC <- ge.all[selPbs, ]
sdat <- selTopVarGenes(ge.CRC, MADth=0.5)
}












