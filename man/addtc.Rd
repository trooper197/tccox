\name{addtc}
\alias{addtc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Determines TC intervals and breaks up dataset into TC intervals
}
\description{
  First the starting and stopping TC intervals are determined according to the expected event criteria; ITC intervals are truncated if necessary or dropped if meaningless. Next the dataset is broken up according to the TC intervals.Finally, the dataset is compressed to only the meaningful intervals.
}
\usage{
addtc(dataset, ncov, maxfollow, start_times, stop_times, min_future_events, numevents,
 nperson, nmaxint, maxobs, interval_width, nitc_start, itc_start_endpoint, nitc_stop,
 itc_stop_endpoint, tti, tts, followup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
     data.frame organized as expected by tc()
}
  \item{ncov}{
     number of baseline covariates (including treatment) to be included in model
}
\item{maxfollow}{
     maximum followup for any subject in dataset
}
  \item{start_times}{
     vector of ordered times when starting of treatment occurs in dataset
}
  \item{stop_times}{
     vector of ordered times when stopping of treatment occurs in dataset
}
  \item{min_future_events}{
     minimum number of events expected of future starters(stoppers) of treatment for determining upper bound on starting(stopping) TC intervals
}
  \item{numevents}{
     number of events in dataset
}
  \item{nperson}{
     number of subjects in dataset
}
  \item{nmaxint}{
     maximum number of TC intervals allowed
}
  \item{maxobs}{
     maximum number of observations (intervals of time) allowed for dataset
}
  \item{interval_width}{
     width of the TC intervals
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
     Vector containing the endpoints of the ITC starting intervals
}
  \item{tti}{
     vector of same length as dataset containing times when starting occurs or 0 if subject does not start
}
  \item{tts}{
     vector of same length as dataset containing times when stopping occurs or 0 if subject does not stop
}
  \item{followup}{
     vector of same length as dataset containing followup times
}
}
\value{
  \item{dataset }{data.frame with dataset broken up according to TC intervals}
  \item{tstartp }{matrix whose columns are the ITC starting covariate values}
  \item{tstopp }{matrix whose columns are the ITC stopping covariate values}
  \item{nstartint }{number of TC starting intervals}
  \item{startint }{vector containing the TC starting interval endpoints}
  \item{nstopint }{number of TC stopping intervals}
  \item{stopint }{vector containing the TC stopping interval endpoints}
  \item{nitc_start }{number of ITC starting intervals}
  \item{itc_start_endpoint }{vector containing the ITC starting interval endpoints}
  \item{nitc_stop }{number of ITC stopping intervals}
  \item{itc_stop_endpoint }{vector containing the ITC stopping interval endpoints}
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
