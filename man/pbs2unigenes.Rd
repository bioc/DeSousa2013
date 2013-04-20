\name{pbs2unigenes}
\alias{pbs2unigenes}
\title{
Generate a mapping file between probesets and unique gene symbols
}                                                                        
\description{
This function takes the expression data of top variable (MAD) genes to generate 
a mapping file between probesets and unique gene symbols.
}
\usage{
pbs2unigenes(ge.CRC, sdat)
}
\arguments{
\item{ge.CRC}{
a numeric matrix of expression data of genes expressed in at least one sample. 
}
\item{sdat}{
a matrix of median centered expression values of top variable probesets. 
}
}
\details{
To facilitate the use of the classifier on other platforms, we collapsed the 
expression levels for probesets to genes. For each gene the probeset with 
highest overall expression was selected. 
}
\value{
This function will save \code{uniGenes} (a character vector of unique gene 
symbols named by probesets) to a file called \code{uniGenes.RData} in the 
current working directory.
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
\code{\link[DeSousa2013:selTopVarGenes]{selTopVarGenes}}
}
\examples{
library(hgu133plus2.db)
data(ge.CRC, package="DeSousa2013")
data(dat, package="DeSousa2013")
ge.CRC <- ge.all[selPbs, ]
uniGenes <- pbs2unigenes(ge.CRC, sdat)
}












