\name{mK}
\alias{mK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Averging Matrix 
}
\description{
Construct a n-by-n averaging matrix.}
\usage{
mK(n)
}
\arguments{
  \item{n}{
a numeric describes the dimension of the averaging matrix.
}
}
\value{
This function returns a n-by-n square matrix with all elements equal 1/n. 
}
\references{
John J, Williams E (1987). \emph{Cyclic and computer generated Designs}. Second edition. Chapman
& Hall.
}
\author{
Kevin Chang
}

\examples{
mK(10)
}
