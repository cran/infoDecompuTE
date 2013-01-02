\name{makeTreatProjectors}
\alias{makeTreatProjectors}
\title{
Make the Treatment Projectors}
\description{
Construct a list of treatment contrast matrices.
}
\usage{
makeTreatProjectors(design.df, trtCols, effectsMatrix, trt.contr, contr.matrix)
}
\arguments{
  \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
  \item{trtCols}{
a list of treatment terms generated from the terms function.
}
  \item{effectsMatrix}{
a matrix contains the 
}
  \item{trt.contr}{
a list of treatment contrast vectors, this allows the user to specify the contrasts for each treatment factor. Note that if this arguement is used, it is necessary to specify the contrasts for every treatment factor with the same order as \code{fixed.terms}. Default is \code{NA}, which uses the C matrix described by John and Williams (1987).
}
\item{contr.matrix}{
a logical allows the users to use the contrast matrices or contrast vectors. If \code{contr.matrix = FALSE}, the outputs will use the contrast vectors each contributes one DF.  If the \code{contr.matrix} arguement is \code{TRUE} the output will conbine the contrast vector for each treatment factors. The default depends on whether \code{trt.contr} has defined, if \code{trt.contr} has defined, then \code{contr.matrix} becomes \code{FALSE}.  
}
}
\value{
This function returns a list of treatment contrast matrices.
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

		T =  makeTreatProjectors(design1, trtTerm, effectsMatrix, trt.contr = NA , contr.matrix = TRUE)
		
		#Fit with the pre-defined treatment contrasts
    Trt1 = rep(c(1,-1), each = 4)
    
    Trt2 =  rep(c(1,-1), time = 4)
    
    Trt3 = Trt1*Trt2
  
		T =  makeTreatProjectors(design1, trtTerm, effectsMatrix, 
      trt.contr =list(Trt = list(Trt1 = Trt1, Trt2 = Trt2, Trt3 = Trt3)) , contr.matrix = FALSE)

}
