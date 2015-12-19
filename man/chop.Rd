\name{chop}
\alias{chop}
\docType{data}
\title{Survival of CHOP for diffuse large B cell lymphoma}
\description{
Microarray data for DLBCL patients undergoing CHOP treatment.
}
\usage{data(chop)}
\format{
  The format is:
 num [1:181, 1:3835] 
}
\details{
Microarray data of DLBCL of 181 patients treated with a combination chemotherapy with cyclophosphamide, doxorubicin, vincristine
and prednisone (CHOP). The original data have 54675 probe sets or covariates. Due to the nature of high-dimensional data, a preselection procedure was conducted to filter out the genes with lower variations if a sample variance for a gene was smaller than the 10th percentile for that gene.
The first column if the survival times. The second column is an indicator whether an the survival time was observed or right censoring occurred. 0=alive, 1=dead. There are 3833 genes after the filtering process. 
}
\source{
Lenz, et al. (2008). Stromal gene signatures in large-{B}-cell lymphomas. \emph{New England Journal of Medicine},
\bold{359(22)}, 2313--2323
}

\examples{
data(chop)
str(chop)
}
\keyword{datasets}
