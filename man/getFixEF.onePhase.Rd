\name{getFixEF.onePhase}
\alias{getFixEF.onePhase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the Fixed Components' coefficients and Efficiency Factors of Single-Phase Experiments.
}
\description{
Construct a matrix contain the coefficients of the fixed Components and the average efficiency factors of single-phase experiments.
}
\usage{
getFixEF.onePhase(effFactors, trt.Coef, T, Rep, table.legend)
}
\arguments{
  \item{effFactors}{
a list numeric vector from \code{trtProkMat} function.
}
  \item{trt.Coef}{
a numeric vector from  \code{getTrtCoef} function.
}
  \item{T}{
a list of matrix from  \code{makeTreatProjectors} function.
}
  \item{Rep}{
a numeric vector from  \code{getReplicationList} function.
}
  \item{table.legend}{
a logical allows the users to use the legend for the variance components of the ANOVA table for a large design. Default is \code{FALSE}, which uses the original names.  
}
}
\value{
This function returns a matrix.
}
\author{
Kevin Chang
}
