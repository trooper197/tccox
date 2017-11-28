\name{itcfitter}
\alias{itcfitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fits ITC Cox model using coxph() function
}
\description{
  Adds ITC covariates to dataset

  Calls coxph() function to fit ITC model
}
\usage{
itcfitter(dataset, ncov, cov_names, nitc_start, itc_start_endpoint, nitc_stop,
 itc_stop_endpoint, tstartp, tstopp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
     data.frame organized as expected by tc()
}
  \item{ncov}{
     number of baseline covariates (including treatment) to be included in model
}
  \item{cov_names}{
     vector of baseline covariate names (including treatment)
}
  \item{nitc_start}{
     number of ITC starting intervals
}
  \item{itc_start_endpoint}{
     vector containing the endpoints of the ITC starting intervals
}
  \item{nitc_stop}{
     number of ITC stopping intervals
}
  \item{itc_stop_endpoint}{
     vector containing the endpoints of the ITC starting intervals
}
  \item{tstartp}{
     matrix whose columns are the ITC starting covariate values
}
  \item{tstopp}{
     matrix whose columns are the ITC stopping covariate values
}
}
\value{
  \item{treatcoef_tc3 }{fitted coefficient of treatment term}
  \item{treatse_tc3 }{fitted SE of treatment term}
  \item{fit_itc}{fit of ITC model}
  \item{cov_names3}{vector containing the covariate names of the ITC model}
}
\references{
 Troendle, JF, Leifer, E, Zhang Z, Yang, S, and Tewes H (2017) How to Control for Unmeasured Confounding in an Observational Time-To-Event Study With Exposure Incidence Information: the Treatment Choice Cox Model.  Statistics in Medicine 36: 3654-3669.
}
\author{
  James F. Troendle
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
