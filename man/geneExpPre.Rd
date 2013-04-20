\name{geneExpPre}
\alias{geneExpPre}
\title{
Preprocessing of the AMC-AJCCII-90 microarray data
}                                                                        
\description{
Preprocessing of the AMC-AJCCII-90 microarray data
}
\usage{
geneExpPre(celpath, AMC_sample_head)
}
\arguments{
\item{celpath}{
the path to a directory where .CEL files of the AMC-AJCCII-90 set are placed.
}
\item{AMC_sample_head}{
a data frame include clinical information and mapping information between 
microarray ids and sample ids (details in \code{\link[DeSousa2013:AMC_sample_head]{AMC_sample_head}}).
}
}
\details{
This function reproduces the preprocessing of the AMC-AJCCII-90 colon cancer 
microarrays, together with 13 adenomas and 6 normal samples. The microarrays of 
the 90 cancer samples and normal samples were generated in one batch, and the 
adenoma samples in a different batch. These two batches were first 
normalized separately by \code{\link[frma:frma]{frma}}, and then corrected for batch 
effect using \code{\link[sva:ComBat]{ComBat}}. The \code{\link[frma:barcode]{barcode}} 
function is used to call present probesets in the two batches, 
respectively. Finally, probesets present in any one of the cancer, adenoma or 
normal samples are selected for the following analysis.
}
\value{
This function will return a list including \code{ge.all} (a numeric matrix 
including expression levels of all probesets of the 90 colon cancer samples) 
and \code{selPbs} (a character vector of probeset ids that are present in any 
sample).
}
\references{
De Sousa E Melo, F. and Wang, X. and Jansen, M. et al. Poor prognosis colon 
cancer is defined by a molecularly distinct subtype and precursor lesion. 
accepted

McCall, Matthew N and Bolstad, Benjamin M and Irizarry, Rafael A (2010). 
Frozen robust multiarray analysis (fRMA). Biostatistics, 11(2), 242-253. 

McCall, Matthew N and Uppal, Karan and Jaffee, Harris A and Zilliox, Michael J 
and Irizarry, Rafael A (2011). The Gene Expression Barcode: leveraging public 
data repositories to begin cataloging the human and murine transcriptomes. 
Nucleic acids research, 39(suppl 1), D1011-D1015. 

Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in 
microarray expression data using Empirical Bayes methods. Biostatistics 
8(1):118-127.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\examples{
\dontrun{
##Please make sure that all colon cancer, adenomas as well as normal samples 
##have been collected from GEO before running this function.  
##This function may take a long time (> half an hour) to finish.
##Please use R-2.15 and corresponding frma package to run the function. 
##Different versions of dependent microarray data preprocessing packages and annotation 
##packages may give slightly different results.
##Please contact the author if there is any confusion
data("AMC")
ge.pre <- geneExpPre(celpath=".", AMC_sample_head)
}
}












