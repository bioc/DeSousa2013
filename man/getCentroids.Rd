\name{getCentroids}
\alias{getCentroids}
\title{
Retrieve centroids of a PAM classifier
}                                                                        
\description{
This function retrieves centroids of a PAM classifier
}
\usage{
getCentroids(fit, data, threshold)
}
\arguments{
\item{fit}{
the fitting result of a PAM classifier to training data.
}
\item{data}{
the training data used for building the PAM classifier
}
\item{threshold}{
the threshold for PAM classification
}
}
\details{
This is an internal function called by \code{\link[DeSousa2013:buildClassifier]{buildClassifier}} 
to retrieve the centroids of PAM classifier built using the AMC-AJCCII-90 
data set.
}
\value{
This function will return a numeric matrix of centroids of the PAM 
classifier. 
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
\code{\link[DeSousa2013:getCentroids]{getCentroids}}
}












