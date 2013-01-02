\name{getTrtCoef}
\alias{getTrtCoef}
\title{
Get the Treatment Coefficient
}
\description{
Compute the overall treatment coefficients.
}
\usage{
getTrtCoef(design.df, trtCols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
  \item{trtCols}{
a list of treatment terms generated from the terms function.
}
}
\value{
This function returns a the numeric vector.
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

  trt.str = "Trt"
  
  	fT = terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

	trtTerm = attr(fT,"term.labels")
		effectsMatrix = attr(fT,"factor") 
		
		trt.Coef = getTrtCoef(design1, trtTerm)

}
