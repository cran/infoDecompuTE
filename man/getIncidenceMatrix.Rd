\name{getIncidenceMatrix}
\alias{getIncidenceMatrix}
\title{
Get the Treatment Incidence Matrix
}
\description{
Construct the treatment incidence matrix base on the experimental design.
}
\usage{
getIncidenceMatrix(design.df, trtCols)
}
\arguments{
  \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
  \item{trtCols}{
a list of treatment terms generated from the terms function.
}
}
\value{
This function returns a matrix describing the design of overall treatment structure.
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
    
	N =  getIncidenceMatrix(design1, trtTerm)
       

}
