\name{adjustMissingLevels}
\alias{adjustMissingLevels}
\title{
Adjust the Missing Levels
}
\description{
Adjust for appropriate syntax describing the effects matching the structural formula.
}
\usage{
adjustMissingLevels(design.df, str.for)
}
\arguments{
\item{design.df}{
a data frame containing the experimental design. Requires every column be a \code{\link{factor}}.
}
\item{str.for}{
a single string of characters containing the structural formula using the Wilkinson-Rogers' syntax.
}
}

\value{
A list containing a data frame with the experimental design and a single string of characters containing the structural formula.
}
\author{
Kevin Chang
}
\examples{
design.df = data.frame( Blk = factor(1:16),
                      	Ani = factor(c(	1,1,2,2,
                                      	1,1,2,2,
                                      	1,1,2,2,
                                      	1,1,2,2)),
                      	Trt = factor(c(	1,2,3,4,
                                      	1,2,3,4,
                                      	1,2,3,4,
                                     		1,2,3,4)))
 
adjustMissingLevels(design.df, str.for = "Ani/Trt") 
}
