\name{getVMat.onePhase}
\alias{getVMat.onePhase}
\title{
Get Variance Matrix for Single-phase Experiments
}
\description{
Construct the matrix for each variance components for the single-phase experiments.
}
\usage{
getVMat.onePhase(Z.Phase1, design.df, var.comp = NA)
}
\arguments{
  \item{Z.Phase1}{
a list of block design matrix from \code{makeBlkDesMatrix} function.
}
 \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
   \item{var.comp}{
a vector of characters containing the variance components of interest this allows the user to specify the variance components to be shown on the ANOVA table. This also allows the user to specify artificial stratum to facilitate decomposition. Default is \code{NA}, which uses every random factor as the variance components with the first phase's variance components in \code{random.terms1} appear before the second phase's variance components in \code{random.terms2}.
}
}
\value{
This function returns a list of matrices.
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

    V = getVMat.onePhase(Z, design1)
}
