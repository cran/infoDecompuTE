\name{getVcCoefMS.onePhase}
\alias{getVcCoefMS.onePhase}
\title{
Get Variance Components' Coefficients and Mean Squares for Single-Phase Experiments
}
\description{
Construct a matrix contain the coefficients of the variance componenets of single-phase experiments. The mean squares is calculated is the \code{response} arguement is used.
}
\usage{
getVcCoefMS.onePhase(PNTginvATNP, design.df, v.mat, response, table.legend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PNTginvATNP}{
a list of matrix from the \code{blkProkMat} function.
}
 \item{design.df}{
a data frame containing the design of the experiment. Every column must be factors. 
}
  \item{v.mat}{
a list of matrix from the \code{getVMat.onePhase} function.
}
  \item{response}{
a numeric vector contains the responses from the experiment. 
}
 \item{table.legend}{
a logical allows the users to use the legend for the variance components of the ANOVA table for a large design. Default is \code{FALSE}, which uses the original names.  
}
}
\value{
This function returns a matrix. 
}
\author{
Kevin Chang
}
