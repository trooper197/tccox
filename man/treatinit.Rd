\name{treatinit}
\alias{treatinit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Gets dataset into standard form
}
\description{
  Creates sequential IDs

  Calculates starting and stopping times
}
\usage{
treatinit(dataset,nmaxint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
     data.frame organized as expected by tc()
}
\item{nmaxint}{
     maximum number of TC intervals allowed
}
}
\value{
  \item{dataset }{data.frame with dataset broken up according to TC intervals}
  \item{obsnum }{number of intervals of time in dataset}
  \item{nperson }{number of subjects in dataset}
  \item{maxobs }{maximum number of observations (intervals of time) allowed for dataset}
  \item{numevents }{number of events in datsaet}
  \item{start_times }{vector of ordered times when starting of treatment occurs in dataset}
  \item{stop_times }{vector of ordered times when stopping of treatment occurs in dataset}
  \item{followup }{vector of same length as dataset containing followup times}
  \item{tti }{vector of same length as dataset containing times when starting occurs or 0 if subject does not start}
  \item{tti }{vector of same length as dataset containing times when stopping occurs or 0 if subject does not stop}
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
