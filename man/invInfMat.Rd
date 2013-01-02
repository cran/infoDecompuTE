\name{invInfMat}
\alias{invInfMat}
\title{
Invert the Information Matrix 
}
\description{
Using the eigenvalue decomposition method to invert the infromation matrix.
}
\usage{
invInfMat(C, N, T)
}
\arguments{
  \item{C}{
a matrix of block projector for a single stratum.
}
  \item{N}{
a matrix of treatment incidence matrix from \code{getIncidenceMatrix} function.
}
   \item{T}{
a list of matrix from  \code{makeTreatProjectors} function.
}
}
\value{
This function returns a matrix.
}
\references{
Nelder JA (1965b). "The Analysis of Randomized Experiments with Orthogonal Block Structure. II. Treatment Structure and the General Analysis of Variance." \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences}, 283(1393), 163-178.
}
\author{
Kevin Chang
}
