\name{getVMat.onePhase}
\alias{getVMat.onePhase}
\alias{getVMat.twoPhase}

\title{
Get the Variance Matrices for Sngle-Phase or Two-Phase experiment
}
\description{
Construct the matrix for each variance components for the single-phase or two-phase experiment.
}
\usage{
getVMat.onePhase(Z.Phase1, design.df, var.comp = NA)
getVMat.twoPhase(Z.Phase1, Z.Phase2, design.df, var.comp = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z.Phase1}{
a list of block design matrix from \code{makeBlkDesMat} function from Phase 1 block structure.
}
  \item{Z.Phase2}{
a list of block design matrix from \code{makeBlkDesMat} function from Phase 2 block structure.
}
   \item{design.df}{
a data frame containing the experimental design. Requires every column be a \code{\link{factor}}.
}
   \item{var.comp}{
a vector of characters containing the variance components of interest this allows the user to specify the variance components to be shown on the ANOVA table. This also allows the user to specify artificial stratum to facilitate decomposition. Default is \code{NA}, which uses every random factor as the variance components with the first phase's variance components in appear before the second phase's variance components.
}
}
\value{
A list of matrices.
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
		Z = makeBlkDesMat(design1, rev(attr(rT,"term.labels")))

    V = getVMat.onePhase(Z, design1)
    
    design2 <- local({ 
  Run = as.factor(rep(1:4, each = 4))
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8,
                            3,4,1,2,
                            7,8,5,6)])
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 4)])
  Trt = as.factor(letters[c(1,2,1,2,
                            2,1,2,1,
                            1,2,1,2,
                            2,1,2,1)])
  data.frame(Run, Ani, Tag, Trt)
})

    blk.str1 = "Ani"
    blk.str2 = "Run"
   
	rT1 = terms(as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE) 
	#random terms phase 1
	rT2 = terms(as.formula(paste("~", blk.str2, sep = "")), keep.order = TRUE) 
	#random terms phase 2

	blkTerm1 = attr(rT1,"term.labels")
	blkTerm2 = attr(rT2,"term.labels")

	Z1 = makeBlkDesMat(design2, rev(blkTerm1))
	Z2 = makeBlkDesMat(design2, rev(blkTerm2))

	V = getVMat.twoPhase(Z1, Z2, design2, var.comp = NA)
}
