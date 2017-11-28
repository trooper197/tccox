\name{coxfitter}
\alias{coxfitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fits Cox model using coxph() function
}
\description{
  Calls coxph() function to fit Cox model
}
\usage{
coxfitter(dataset1, ncov, cov_names)
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
     vector of baseline covariate names (including treatment)
}
}
\value{
  \item{treatcoef_cox1 }{fitted coefficient of treatment term}
  \item{treatse_cox1 }{fitted SE of treatment term}
  \item{fit_cox}{fit of Cox model}
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
