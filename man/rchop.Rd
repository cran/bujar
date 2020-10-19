\name{rchop}
\alias{rchop}
\docType{data}
\title{Survival of R-CHOP for diffuse large B cell lymphoma}
\description{
Microarray data for DLBCL patients undergoing R-CHOP treatment.
}
\usage{data(rchop)}
\format{
  The format is:
 num [1:233, 1:3835] 
}
\details{
Microarray data of DLBCL of 233 patients treated with the current gold standard R-CHOP including rituxima immunotherapy in addition to the chemotherapy CHOP.  The original data have 54675 probe sets or covariates. Due to the nature of high-dimensional data, a preselection procedure was conducted to filter out the genes to match those in \code{chop}.
The first column if the survival times. The second column is an indicator whether an the survival time was observed or right censoring occurred. 0=alive, 1=dead. There are 3833 same genes as in \code{chop}. The data set is used to validate the prediction accuracy for models developed using training data \code{chop}. 
}
\source{
Lenz, et al. (2008). Stromal gene signatures in large-{B}-cell lymphomas. \emph{New England Journal of Medicine},
\bold{359(22)}, 2313--2323
}

\examples{
data(rchop)
str(rchop)
}
\keyword{datasets}
