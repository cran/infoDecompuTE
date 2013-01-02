\name{trtProkMat}
\alias{trtProkMat}
\title{
Compute the Efficiency Factors by Treatment Information Decomposition  
}
\description{
Pre- and post-multiply NTginvATN by block projection matrices of a single stratum. It will then compute the eigenvalues and the averege efficiency factors. 
}
\usage{
trtProkMat(z, T, N, Rep)
}
\arguments{
  \item{z}{
a matrix of block projector for a single stratum.
}
  \item{T}{
a list of treatment contrast matrix from \code{makeTreatProjectors} function.
}
  \item{N}{
a matrix of treatment incidence matrix from \code{getIncidenceMatrix} function.
}
   \item{Rep}{
a numeric vector of replication ratio.
}
}
\value{
This function will return a list of numeric vectors.
}
\author{
Kevin Chang
}
