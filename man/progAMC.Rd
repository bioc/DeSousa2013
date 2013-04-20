\name{progAMC}
\alias{progAMC}
\alias{figKM}
\title{
Prognosis of the three subtypes in the AMC data set
}                                                                        
\description{
These functions performs survival analysis for the AMC data set, and generate 
a kaplan-meier plot to diagnose the prognosis of three subtypes.
}
\usage{
progAMC(AMC_CRC_clinical, AMC_sample_head, clu.pred)
figKM(surv, survstats)
}
\arguments{
\item{AMC_CRC_clinical}{
a data frame including  include clinical information (details in \code{\link[DeSousa2013:AMC_CRC_clinical]{AMC_CRC_clinical}}).
}
\item{AMC_sample_head}{
a data frame include mapping information between microarray ids and sample ids 
(details in \code{\link[DeSousa2013:AMC_sample_head]{AMC_sample_head}}).
}
\item{clu.pred}{
a numeric vector of classification labels named by colon cancer samples.
}
\item{surv}{
the result of \code{\link[survival:survfit]{survfit}} containing survival curves 
of the AMC data set.
}
\item{survstats}{
the result of \code{\link[survival:survdiff]{survdiff}} for the AMC data set.
}
}
\details{
The function \code{progAMC} performs survival analysis for the AMC data set and 
compares the prognosis of the three subtypes. 
The function \code{figKM} helps visualize the prognosis in a kaplan-meier plot.
}
\value{
This function \code{progAMC} will save \code{surv} (result of \code{\link[survival:survfit]{survfit}} 
containing survival curves of the AMC data set), \code{survstats} (result of 
\code{\link[survival:survdiff]{survdiff}}), \code{data4surv} (the data to perform 
survival analysis) to a file called \code{survival.RData} in the current 
working directory.
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
library(survival)
data(AMC, source="DeSousa2013")
data(predAMC, source="DeSousa2013")
prog <- progAMC(AMC_CRC_clinical, AMC_sample_head, clu.pred)
surv <- prog$surv
survstats <- prog$survstats
figKM(surv, survstats)
}












