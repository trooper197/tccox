\name{ptcfitter}
\alias{ptcfitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fits PTC Cox model using coxph() function
}
\description{
  Calls coxph() function to fit PTC model
}
\usage{
ptcfitter(dataset, ncov, cov_names)
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
}
\value{
  \item{treatcoef_tc1 }{fitted coefficient of treatment term}
  \item{treatse_tc1 }{fitted SE of treatment term}
  \item{fit_ptc}{fit of PTC model}
  \item{cov_names1}{vector containing the covariate names of the PTC model}
}
\references{
 Troendle, JF, Leifer, E, Zhang Z, Yang, S, and Tewes H (2017) How to Control for Unmeasured       Confounding in an Observational Time-To-Event Study With Exposure Incidence Information: the Treatment Choice Cox Model.  Statistics in Medicine 36: 3654-3669.
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
