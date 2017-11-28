\name{ptc}
\alias{ptc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prepares and fits a PTC Cox model
}
\description{
  Calls treatinit() to prepare dataset

  Calls addtc() to create TC intervals

  Calls ptcfitter() to fit PTC model
}
\usage{
ptc(dataset, ncov, cov_names, maxfollow, nmaxint, interval_width, min_future_events)
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
\item{maxfollow}{
     maximum followup for any subject in dataset
}
  \item{nmaxint}{
     maximum number of TC intervals allowed
}
  \item{interval_width}{
     width of the TC intervals
}
  \item{min_future_events}{
     minimum number of events expected of future starters(stoppers) of treatment for determining upper bound on starting(stopping) TC intervals
}
}
\value{
  \item{fit_ptc }{fit of PTC model}
  \item{nstartint }{number of TC starting intervals}
  \item{startint }{vector containing the TC starting interval endpoints}
  \item{nstopint }{number of TC stopping intervals}
  \item{stopint }{vector containing the TC stopping interval endpoints}
  \item{cov_names1 }{vector containing the covariate names of the model}
  \item{nperson }{number of subjects in dataset}
  \item{numevents }{number of events in dataset}
  \item{medianfollowup }{median followup for subjects in dataset}
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
