% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotfunctions.R
\name{SurvCF.KMPlot}
\alias{SurvCF.KMPlot}
\title{SurvCF.KMPlot}
\usage{
SurvCF.KMPlot(x)
}
\arguments{
\item{x}{An SurvCF object created with SurvCF()}
}
\description{
KM plot of a SurvCF object using GGally and some default settings.
}
\examples{

sim.df <- simStudy()
x <- SurvExt(Surv(os.t, os.e) ~ I(x.trt==1),
   Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
   AdminCensTime = t.censor,
   data = sim.df)
SurvExt.KMPlot(x)
}
\keyword{RPSFT,}
\keyword{Survival}
