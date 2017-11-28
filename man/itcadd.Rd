\name{itcadd}
\alias{itcadd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Determines potential ITC intervals according to expected event criteria or as fixed
}
\description{
  Either speciifies the ITC interval endpoints as fixed or else determines candidates based on expected event criteria
}
\usage{
itcadd(dataset, nmaxint, interval_width, min_exp_events, nitc_fixed, n_start_fixed,
 n_stop_fixed, interval_stop_beginning)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
     data.frame organized as expected by tc()
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
  \item{nitc_fixed}{
     indicator that potential ITC intervals are fixed
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
