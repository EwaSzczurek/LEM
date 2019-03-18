## coll.model - collapsed graph
score1 <- function( coll.model, Y, parameter.estimation="linear.reg", verbose){
	
	if (	is.vector(Y) ){
		Y = matrix(Y, ncol = 1)
	} 
	### Y is a matrix with columns corresponding to different response vectors. 
	### We learn different betas and sigmas for each column of Y, 
	### the likelihood is a product of the likelihoods over the columns 

	E=coll.model$E
	G=coll.model$G
	
		
	
	models <- apply(Y, 2, inference, G = G, E = E, parameter.estimation=parameter.estimation, verbose  =verbose )	
	
	
	LL = NULL
	beta = NULL
	c.vector = NULL
	b = NULL

	for (model in models){
		LL = c(LL, model$LL)
		beta = cbind(beta, model$beta)
		c.vector = c(c.vector, model$c)
		b  = c( b, model$b )	
	}
	
	### We sum the log likelihoods across the response vectors
	LL.all = sum(LL)
	
	if (parameter.estimation == "linear.reg"){
		### Deal with different number of parameter values using BIC
		n = nrow(E) # number of experiments
		k = ncol(E) # number of parameters to estimate
		score = -log(n)*k + 2*LL.all
	}else{
		score=LL.all
	}
	if (verbose){
		print("scored:")
		print(G)
		print(score)
	}
	list( score=score, beta=beta, sigma=sqrt(1/c.vector), b=b )

}


score <- function(coll.models, Y, parameter.estimation,  verbose=TRUE) {
	   
	### Compute the score for each model
	results <- sapply(coll.models, score1, Y,  parameter.estimation, verbose) 
	
  	s <- unlist(results["score",])
  	beta  <- results["beta",]
	score = results["score", ]
	sigma     <- results["sigma",] 	  
	b = results["b",]
	  	
 	# winning model
  	best <- which.max(s)       
  	winner <- coll.models[[best]]$G
  	diag(winner) <- 1  
  	beta.winner <- beta[[best]]
  	rownames(beta.winner)=colnames(winner)
  	b.winner=b[[best]]
  	sigma.winner = sigma[[best]]
  	score.winner <- score[[best]]
  	
  	gR = winner

	res <- list(graph = gR, beta = beta.winner,  score = score.winner, sigma = sigma.winner, b = b.winner,  all.scores = s)
  	
  	return(res)  
}
 