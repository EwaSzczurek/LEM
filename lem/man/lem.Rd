\name{lem}
\alias{lem}
\title{Linear Effects Models - main function}
\description{The main function to perform model estimation from data}
\usage{
	lem(Y, E, Gs=NULL, inference="greedy", parameter.estimation="linear.reg", verbose=FALSE, greedy.init.no = 10, collapse.res = TRUE)
}

\arguments{
  \item{Y}{Effect measurements: can be a vector (for one measured phenotype) or a matrix (the different phenotypes then stand for columns). Contains continuous values. Each element of the vector or each row of the matrix corresponds to a perturbation experiment under which the measurement was taken.}
  \item{E}{Binary experiment matrix (called the perturbation matrix). Each row corresponds to an experiment, each column to a pathway gene. The entry in row e and column g is equal 1 if experiment e targets gene g. }
  \item{Gs}{Optional list of graphs (as matrices) to search the best from. Default = \code{NULL}.}
  \item{inference}{\code{search} to use exhaustive enumeration over provided set of graphs or over all possible, \code{greedy} for greedy hillclimbing}
   \item{parameter.estimation}{ \code{linear.reg} (default) for linear regression and evaluation using BIC, \code{bayes.linear.reg} for bayesian linear regression.}
  \item{verbose}{Allows printout of progression statements. Default: FALSE} 
  \item{greedy.init.no}{Number of initial graphs to initialize the greedy search with. The set of greedy.init.no graphs is sampled at random. Always to this initial set an empty graph (separate nodes with no edges) is added, since it proved to often be a good sarting point in simulations. }
  \item{collapse.res}{ States whether the resulting graph should have cliques collapsed into single nodes. Defaults to TRUE.}
}
\details{  
  For the method to provably work, we need at least one perturbation experiment per each single gene, and at least one experiment per each pair. Check your experiments matrix E! The effect measurements can be any real vectors. 
  Bayesian linear regression (for parameter.estimation = "bayes.linear.reg") takes much longer to compute than linear regression (parameter.estimation = "linear.reg").)   
}
\value{
Object of class "lem", with attributes:
  \item{graph}{inferred directed pathway gene graph (represented as an adjacency matrix)}   
  \item{score}{scores of model(s). In the case when bayesian regression was used for inference, it is the marginal log likelihood. In the case when linear regression was used, it is the Bayesian Information Criterion (BIC).}   
  \item{beta}{estimate of the vector of contributions of the genes in the pathway graph to the phenotype}
  \item{sigma}{standard deviation of the phenotype}
  \item{b}{ For parameter.estimation = "bayes.linear.reg" it is the prior parameter b. By default, for parameter.estimation = "linear.reg", b is NULL. }  
  \item{all.scores}{vector of scores of all models considered in the inference of the pathway graph.}
  \item{inference}{Method used for infering the graph. Either "greedy" for greedy search or "search" for exhaustive search over user-defined set of graphs or over all possible graphs of given size.}
  \item{parameter.estimation}{Method used for estimating model parameters and scoring the model. Either "linear.reg" for linear regression or "bayes.linear.reg" for bayesian linear regression.}
  
}
\references{
Szczurek, E & Beerenwinkel, N., Linear Effects Models from Combinatorial Perturbation Data. Bioinformatics, 2016
}

\author{Ewa Szczurek}
\examples{
	% ### simulate the data and the model
	set.seed(111)
  sim.set <- simulateLEM(n=5, b = 1, sigma = 0.01, rep.no = 5, Y.no=1) 
  lem.res <- lem( sim.set$Y, sim.set$E, inference = "greedy", 
  parameter.estimation = "bayes.linear.reg",verbose=FALSE)
  print("The true simulated graph:")
  print(sim.set$G)
  print("The learned model")
 print(lem.res)
 
 }
\keyword{graphs}
\keyword{models}
