\name{cox}
\alias{cox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calls coxfitter() function
}
\description{
  Calls coxfitter() function to fit Cox model
}
\usage{
cox(dataset1, ncov, cov_names, nmaxint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset1}{
     data.frame organized as expected by tc()
}
  \item{ncov}{
     number of baseline covariates (including treatment) to be included in model
}
  \item{cov_names}{
     vector of baseline covariate names
}
\item{nmaxint}{
     maximum number of TC intervals allowed (ignored)
}
}
\value{
  \item{fit_cox}{fit of Cox model}
  \item{nperson }{number of subjects in dataset1}
  \item{numevents }{number of events in dataset1}
  \item{medianfollowup }{median followup for subjects in dataset1}
}
\references{
 Troendle, JF, Leifer, E, Zhang Z, Yang, S, and Tewes H(2017) How to Control for Unmeasured Confounding in an Observational Time-To-Event Study With Exposure Incidence Information: the Treatment Choice Cox Model.  Statistics in Medicine 36: 3654-3669.
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
