\name{getVCs.onePhase}
\alias{getVCs.onePhase}
\title{
Get the Variance Components of Single-Phase Experiments
}
\description{
Computes the coefficients of the variance components for the expected mean squares for single-phase experiments. The function accepts a data frame of the experimetnal design with the structural formulae of the block and treatment factors. Two tables containing the variance components of the random effects and fixed effects are returned.  
}
\usage{
getVCs.onePhase(design.df, random.terms, fixed.terms, var.comp = NA,
                trt.contr = NA, table.legend = FALSE)
}

\arguments{
  \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
  \item{random.terms}{
a single string of characters containing the structural formula for the block factors using the Wilkinson-Rogers' syntax. 
}
  \item{fixed.terms}{
a single string of characters containing the structural formula for the treatment factors using the Wilkinson-Rogers' syntax. 
}
  \item{var.comp}{
a vector of characters containing the variance components of interest this allows the user to specify the variance components to be shown on the ANOVA table. This also allows the user to specify artificial stratum to facilitate decomposition. Default is \code{NA}, which uses every random factor as the variance components from \code{random.terms}.
}
  \item{trt.contr}{
a list of treatment contrast matrices, this allows the user to specify the contrasts for each treatment factor. Note that if this arguement is used, it is necessary to specify the contrasts for every treatment factor with the same order as \code{fixed.terms}. Default is \code{NA}, which uses the C matrix described by John and Williams (1987).
}
  \item{table.legend}{
a boolean allows the users to use the legend for the variance components of the ANOVA table for a large design. Default is \code{FALSE}, which uses the original names.  
}
}

\value{
 The values returned depends on the value of the \code{table.legend} arguement. If \code{table.legend = FALSE}, this function will return a list of two data frames. The first data frame contains the random effects and the second data frame contains the fixed effects. If the \code{table.legend} arguement is \code{TRUE}, then it will return a list containing two lists. The first list consists of a data frame of random effects and a character string for the legend. The second list consists of a data frame of fixed effects and a character string for the legend.
}

\references{
John J, Williams E (1987). \emph{Cyclic and computer generated Designs}. Second edition. Chapman
& Hall.

Nelder JA (1965b). "The Analysis of Randomized Experiments with Orthogonal Block Structure. II. Treatment Structure and the General Analysis of Variance." \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences}, 283(1393), 163-178.

Wilkinson GN, Rogers CE (1973). "Symbolic Description of Factorial Models for Analysis of Variance." \emph{Applied Statistics}, 22(3), 392-399.
}

\seealso{
\code{\link{terms}} for more information on the structural formula.
}

\author{
Kevin Chang 
}

\examples{
design <- local({ 
  Run = as.factor(rep(1:4, each = 4))
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8,
                            3,4,1,2,
                            7,8,5,6)])
  Sam = as.factor(as.numeric(duplicated(Ani)) + 1)
  Tag = as.factor(c(114,115,116,117)[rep(1:4, 4)])
  Trt = as.factor(letters[c(1,2,1,2,
                            2,1,2,1,
                            1,2,1,2,
                            2,1,2,1)])
  data.frame(Run, Ani, Sam, Tag, Trt)
})

getVCs.onePhase(design, random.terms = "Ani/Sam", fixed.terms = "Trt") 
                                
} 

\keyword{design}