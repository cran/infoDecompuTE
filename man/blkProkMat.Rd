\name{blkProkMat}
\alias{blkProkMat}
\title{
Construct PNTginvATNP Matrix by Block Information Decomposition 
}
\description{
Pre- and post-multiply NTginvATN by block projection matrices of a single stratum.
}
\usage{
blkProkMat(z, T, N)
}
%- maybe also 'usage' for other objects documented here.
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
}
\value{
This function will return a list of matrices.
}
\author{
Kevin Chang
}

