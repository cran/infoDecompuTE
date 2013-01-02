\name{toLatexTable}
\alias{toLatexTable}
\title{
Generate the Latex scripts for the Latex table
}
\description{
Print the Latex scripts on the screen for the users to output the table from the Latex output.
}
\usage{
toLatexTable(ANOVA, EF, fixed.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ANOVA}{
a list of matrix from the \code{getVcCoefMS.onePhase} or \code{getVcCoefMS.twoPhase} function.
}
 \item{EF}{
a list of matrix from the \code{getFixEF.onePhase} or \code{getFixEF.onePhase} function.
}
  \item{fixed.names}{
a vector of character allows the users to modify symbols for the fixed effects.
}
}
\author{
Kevin Chang
}
