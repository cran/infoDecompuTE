\name{makeBlkDesMatrix}
\alias{makeBlkDesMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Block Design Matrix
}
\description{
Construct a binary matrix representing the block design. The rows are corrosponding to the observations and the columns are corrosponding to the blocks.
}
\usage{
makeBlkDesMatrix(design.df, blkTerm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
  \item{blkTerm}{
a list of block terms generated from the terms function.
}
}
\value{
This function returns a list of the binary matrices.  
}

\author{
Kevin Chang
}


\seealso{
 \code{\link{terms}}
}
\examples{
  design1 <- local({ 
    Ani = as.factor(LETTERS[c(1,2,3,4,
                              5,6,7,8)])
    Trt = as.factor(letters[c(1,1,1,1,
                              2,2,2,2)])
    data.frame(Ani, Trt)
  })

    blk.str = "Ani"
    
		rT = terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE) 

    blkTerm = attr(rT,"term.labels")
		Z = makeBlkDesMatrix(design1, rev(attr(rT,"term.labels")))

}
