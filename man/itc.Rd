\name{itc}
\alias{itc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prepares and fits an ITC Cox model
}
\description{
  Calls treatinit() to prepare dataset

  Calls itcadd() to create ITC intervals

  Calls addtc() to create TC intervals

  Calls itcfitter() to fit ITC model
}
\usage{
itc(dataset, ncov, cov_names, maxfollow, nmaxint, interval_width, min_exp_events,
 min_future_events, nitc_fixed, n_start_fixed, n_stop_fixed,
 interval_stop_beginning)
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
  \item{min_exp_events}{
     minimum number of events expected of subjects in each cell for determining ITC intervals
}
  \item{min_future_events}{
     minimum number of events expected of future starters(stoppers) of treatment for determining upper bound on starting(stopping) TC intervals
}
  \item{nitc_fixed}{
     indicator (0 or 1) that potential ITC intervals are fixed
}
  \item{n_start_fixed}{
     number of fixed ITC starting intervals (only applicable if nitc_fixed=1)
}
  \item{n_stop_fixed}{
     number of fixed ITC stopping intervals (only applicable if nitc_fixed=1)
}
  \item{interval_stop_beginning}{
     smallest ITC stopping interval endpoint (only applicable if nitc_fixed=1)
}
}
\value{
  \item{fit_itc }{fit of ITC model}
  \item{nitc_start }{number of ITC starting intervals}
  \item{itc_start_endpoint }{vector containing the ITC starting interval endpoints}
  \item{nitc_stop }{number of ITC stopping intervals}
  \item{itc_stop_endpoint }{vector containing the ITC stopping interval endpoints}
  \item{nstartint }{number of TC starting intervals}
  \item{startint }{vector containing the TC starting interval endpoints}
  \item{nstopint }{number of TC stopping intervals}
  \item{stopint }{vector containing the TC stopping interval endpoints}
  \item{cov_names3 }{vector containing the covariate names of the model}
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
