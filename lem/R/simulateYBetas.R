

simulateYBetas <- function( G, E, b=1, sigma, Y.no, distr = "rnorm"){
	
	n = ncol(G)
	m = nrow(E)
	
	if (distr=="rnorm"){
	
		betas =	sapply(1:Y.no, function(i){	be = rnorm(n, 0, sd = sqrt(1/b) ); names(be) = colnames(G); be } )
	}else{
		betas =	sapply(1:Y.no, function(i){	be = rgamma(n, shape = 10, rate = 9); signs = sample( c(-1, 1), n, replace=T );  be  = be*signs; names(be) = colnames(G); be } )
	}
	betas[betas==0]<- 1e-05 # add a pseudocount to zero beta entries - zero entries are not admissable
	
	epsilon = matrix( rnorm( m*Y.no, 0, sigma), ncol=Y.no)
	X = ( (E%*%G)>0)
	
	Y = X%*%betas + epsilon
	
	return( list( Y = Y, betas = betas, E = E, G=G))
}