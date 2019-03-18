\name{simulateLEM}
\alias{simulateLEM}

\title{Simulate a Linear Effects Model}
\description{Generates random graph, experiments matrix, beta contributions and data (measurements) }
\usage{
	simulateLEM(n, b = 1, sigma = 0.01, rep.no = 3, Y.no=1, up.to=NULL, collapse=FALSE)
}

\arguments{
	  \item{n}{Number of genes perturbed in experiments (the genes in the pathway modeled by LEM). }
  \item{b}{The value for the hyper-prior parameter b of the distribution of the beta parameters. We assume the betas distribution is a zero mean
isotropic Gaussian with precision b, By default, b=1.}
\item{sigma}{The standard deviation for the Gaussian noise around the simulated measurements. By default, sigma=0.01.}
\item{rep.no}{The number of times each experiment is repeated.}
\item{Y.no}{The number of measurement vectors to simulate. }
\item{up.to}{The upper bound on the number of combinations of perturbations are to be simulated. Can either be set to 2 or 3. }
\item{collapse}{Whether the cliques in the graph should be collapsed.}

  
 
}
\details{  
	 The rows correspond to experiments and columns to perturbed genes. 
	
	For up.to=NULL, the method will simulate experiments required for identifiability. For $n<=5$, the experiments have to at least perturb each single and each pair of nodes. In this case, the generated experiment matrix E will have rep.no *( n + (n choose 2)) rows. For $n> 5$ nodes, experiments perturbing all single, all pairs and all triplets of genes are required for identifiability (yielding rep.no *( n + (n choose 2) + (n choose 3)) rows). 
	
	If the user speifies a specific number up.to, then all single, all pairs, ..., all up.to-combinations will be simulated.
	
For Y.no = 1, one vector of beta coefficients will be sampled, and one vector of measurements (one vector entry for each experiment) will be generated according to the model. For Y.no >1, a matrix of coefficients will be sampled, and one matrix of measurements will be returned, with one column per one measurement vector. This corresponds to assuming that the pathway is always the same, but if a number of phenotypes is measured, each of them will be different and will have different gene contributions.

The function returns a list.
 
}
\value{
 \item{Y }{The simulated measurement vector (matrix for Y.no>1).}
 \item{betas}{The simulated gene contributions vector (matrix for Y.no>1)}
 \item{E}{The experiment matrix.}
 \item{G}{The simulated graph as an adjacency matrix.}  
}

\author{Ewa Szczurek}


\keyword{graphs}
\keyword{models}
