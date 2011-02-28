\name{getVCs.twoPhase}
\alias{getVCs.twoPhase}
\title{
Get the Variance Components of the Two Phase Experiments
}
\description{
Computes the coefficients of the variance components for the expected mean squares for the two phase experiments. Inputting a data.frame of a design with the relationships of the block and treatment factors, it returns two tables containing the variance components of the random effects and fixed effects.  
}
\usage{
getVCs.twoPhase(design.df, random.terms1, random.terms2, fixed.terms, var.comp = NA, trt.contr = NA, table.legend = FALSE)
}
\arguments{
 \item{design.df}{
data.frame contain the design of the experiment 
}
 \item{random.terms1}{
a single string of characters contains the structural formula for the block factors of the first phase using the Wilkinson-Rogers' syntax . 
}
  \item{random.terms2}{
a single string of characters contains the structural formula for the block factors of the second phase using the Wilkinson-Rogers' syntax . 
}
  \item{fixed.terms}{
a single string of characters contains the structural formula for the treatment factors using the Wilkinson-Rogers' syntax . 
}
  \item{var.comp}{
 a vector of characters contains the varaince componenets of interestm this allows the user to specify the variance componenets to be shown on the ANOVA table. Default is NA, which uses every random factors as the variance components with the 1st phase's variance components appear before the 2nd phase's variance components. 
}
  \item{trt.contr}{
a list of treatment contrast matrices, this allows the user to specify the contrasts for each treatment factor. Note that if the user uses this parameter, they need to specify the contrasts for every treatment factor. Default is NA, which uses the C matrix proposed by John and Williams.
}
  \item{table.legend}{
 a logical, allow the users to using the lengend for the variance compoenents of the ANOVA table for a large design. Default is FALSE, which uses the orignal names.    
}
}
\value{
 The values returned depends wether legend arguement has specified. If it has not, this function will return a list of two data frames. The first data frame contains the random effects and the second data frame contains the fixed effects. If the legends arguement has sepcified, then it will return a list containing two lists. The first list consists of a data frame of random effects and a character string for the legend. The second list consists of a data frame of fixed effects and a character string for the legend.
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

getVCs.twoPhase(design, random.terms1 = "Ani/Sam", random.terms2 = "Run", 
fixed.terms = "Trt + Tag")                                    

}
