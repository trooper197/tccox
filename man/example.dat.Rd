\name{example.dat}
\docType{data}
\alias{example.dat}
\title{Simulated data from Observational Study}
\description{
  data.frame containing simulated data from an observational study

     treatment (0 or 1)

     baseage (baseline age in years)

     id (integer identifying subjects)

     start (beginning of each time interval where treatment is constant)

     stop (endpoint of each time interval where treatment is constant)

     status is an indicator (0 or 1) of an event occuring at the corresponding stop time

     ordered by start values within each level of id with first enty for each id having treatment=0
}
\usage{example.dat}
\format{data.frame}
\references{
 Troendle, JF, Leifer, E, Zhang Z, Yang, S, and Tewes H(2017) How to Control for Unmeasured Confounding in an Observational Time-To-Event Study With Exposure Incidence Information: the Treatment Choice Cox Model.  Statistics in Medicine 36: 3654-3669.
}
\keyword{datasets}
