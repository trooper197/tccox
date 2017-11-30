\name{tc}
\alias{tc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Treatment Choice Cox Model
}
\description{
  Builds time-varying covariates needed and fits Treatment Choice Cox models (Parametric Treatment Choice (PTC), Hybrid Treatment Choice (HTC), or Interval Treatment Choice (ITC)) for observational time-to-event studies.
}
\usage{
tc(type = "PTC", dataset, cov_names, maxfollow=100, nmaxint = 80, interval_width = 0.1,
 min_exp_events = 50, min_future_events = 50, nitc_fixed = 0, n_start_fixed = 10,
 n_stop_fixed = 10, interval_stop_beginning = 1.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
     character indicating the type of TC model to be fit ('PTC' for Parametric, 'HTC' for Hybrid, or 'ITC' for Interval)
}
  \item{dataset}{
     data.frame containing the data to be used to fit the TC model

     dataset should have all baseline covariates, starting with treatment (0 or 1), in the leading positions

     following the baseline covariates should be in order the variables: id, start, stop, status

     id is a unique number for each subject

     start is the beginning of each time interval where treatment is constant

     stop is the endpoint of each time interval where treatment is constant

     status is an indicator (0 or 1) of an event occuring at the corresponding stop time

     dataset should be ordered by start values within each level of id

     for each id the first entry should have treatment=0
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
  \item{fit }{fit of TC model}
  \item{nitc_start }{number of ITC starting intervals}
  \item{itc_start_endpoint }{vector containing the ITC starting interval endpoints}
  \item{nitc_stop }{number of ITC stopping intervals}
  \item{itc_stop_endpoint }{vector containing the ITC stopping interval endpoints}
  \item{nstartint }{number of TC starting intervals}
  \item{startint }{vector containing the TC starting interval endpoints}
  \item{nstopint }{number of TC stopping intervals}
  \item{stopint }{vector containing the TC stopping interval endpoints}
  \item{cov_names }{vector containing the covariate names of the model

                    tstart is the cumulative constant starting term (PTC only)

                    tstart1 is the cumulative linear starting term (PTC only)

                    tstop is the cumulative constant stopping term (PTC only)

                    tstop1 is the cumulative linear stopping term (PTC only)

                    tstart0 is the cumulative constant starting term outside of ITC intervals (HTC only)

                    tstop0 is the cumulative constant stopping term outside of ITC intervals (HTC only)

                    treatstartp.# is the #'th ITC starting term (ITC and HTC only)

                    treatstopp.# is the #'th ITC stopping term (ITC and HTC only)}
  \item{nperson }{number of subjects in dataset}
  \item{numevents }{number of events in datsaet}
  \item{medianfollowup }{median followup for subjects in dataset}
}
\references{
 Troendle, JF, Leifer, E, Zhang Z, Yang, S, and Tewes H (2017) How to Control for Unmeasured Confounding in an Observational Time-To-Event Study With Exposure Incidence Information: the Treatment Choice Cox Model.  Statistics in Medicine 36: 3654-3669.
}
\author{
  James F. Troendle
}
\examples{
# Use simulated data in example.dat to build and fit a PTC model
#
require(stats)
require(survival)
cov_names=names(example.dat)[1:2]
z=tc(type="PTC", dataset=example.dat, cov_names = cov_names, min_exp_events = 10,
 min_future_events = 50)
z[[1]]
#
# Use simulated data in example.dat to build and fit an HTC model
#
require(stats)
require(survival)
cov_names=names(example.dat)[1:2]
z=tc(type="HTC", dataset=example.dat, cov_names = cov_names, min_exp_events = 10,
 min_future_events = 50)
z[[1]]
#
# Use simulated data in example.dat to build and fit an ITC model
#
require(stats)
require(survival)
cov_names=names(example.dat)[1:2]
z=tc(type="ITC", dataset=example.dat, cov_names = cov_names, min_exp_events = 10,
 min_future_events = 50)
z[[1]]
}
\keyword{survival}

