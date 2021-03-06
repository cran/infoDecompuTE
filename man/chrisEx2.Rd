% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{chrisEx2}
\alias{chrisEx2}
\title{Randomised Block design consisted of 8 blocks and 2 plots.}
\format{
A data frame with 18 rows and 5 variables: \describe{
\item{Blocks}{Block factor containing 6 levels} \item{Plots}{Plot factor
containing 3 levels} \item{A}{Treatment factor A containing 2 levels}
\item{B}{Treatment factor B containing 2 levels} \item{C}{Treatment factor
C containing 3 levels} }
}
\usage{
chrisEx2
}
\description{
Randomised Block design consisted of 8 blocks and 2 plots.
}
\examples{
 data(chrisEx2)

summaryAovOnePhase(chrisEx2, "Blocks/Plots", "A*B*C")
}
\keyword{datasets}
