\name{projMat}
\alias{projMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Projection Matrix
}
\description{
Compute the projection matrix.
}
\usage{
projMat(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a matrix.
}
}
\value{
This function returns a square matrix.
}

\author{
Kevin Chang
}

\examples{
m = matrix(1, nrow = 10, ncol = 3)
projMat(m) 
}
