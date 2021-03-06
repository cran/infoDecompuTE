% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invInfMat.R
\name{invInfMat}
\alias{invInfMat}
\title{Invert the Information Matrix}
\usage{
invInfMat(C, N, T)
}
\arguments{
\item{C}{a matrix of block projector for a single stratum.}

\item{N}{a matrix representation the smallest unit of block or treatment
effects generated by \code{\link{makeOverDesMat}}.}

\item{T}{a list of contrast matrices from \code{\link{makeContrMat}}.}
}
\value{
This function returns a matrix.
}
\description{
Using the eigenvalue decomposition method to invert the information matrix.
}
\examples{


m <- matrix(rnorm(10), 10, 10)

invInfMat(m, identityMat(10), identityMat(10))
     


}
\references{
Nelder JA (1965b). "The Analysis of Randomized Experiments with
Orthogonal Block Structure. II. Treatment Structure and the General
Analysis of Variance." \emph{Proceedings of the Royal Society of London.
Series A, Mathematical and Physical Sciences}, 283(1393), 163-178.
}
\author{
Kevin Chang
}
