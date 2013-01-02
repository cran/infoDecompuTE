\name{makeBlockProjectors}
\alias{makeBlockProjectors}
\title{
Make Block Projectors 
}
\description{
Construct a list of projection matrices of every stratum. 
}
\usage{
makeBlockProjectors(BlkDesignMatrixList, initial = diag(nrow(BlkDesignMatrixList[[1]])) - mK(nrow(BlkDesignMatrixList[[1]])))
}
\arguments{
  \item{BlkDesignMatrixList}{
a list of block design matrix from the \code{makeBlkDesMatrix}. }
  \item{initial}{
a matrix contains the first projection matrices for preforming the decomposition on. This arguement is essential for the two-phase experiment.
}
}
\value{
This function returns a list of projection matrices of every stratum. 
}
\author{
Kevin Chang
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
		Pb = makeBlockProjectors(Z)
}
