\name{print.lem}
\alias{print.lem}

\title{Linear Effects Models - printing basic info}
\description{Prints the Linear Effects Model parameters}
\usage{
	\method{print}{lem}(x, \dots)
}

\arguments{
  \item{x}{Lem object (result of calling \code{lem}).}
 \item{...}{Stands for nothing, just used for agreement with the generic.}
}
\details{  
 This function will print the adjacency matrix of the graph, the beta coefficients, the score, and the inference and parameter estimation methods used for learning.
 
 
}
\value{
  none  
}

\author{Ewa Szczurek}


\seealso{ \code{\link{lem}}}

\keyword{graphs}
\keyword{models}
