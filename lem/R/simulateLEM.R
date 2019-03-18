

simulateLEM <- function( n, b = 1, sigma = 0.01, rep.no = 3, Y.no=1,  up.to=NULL,  collapse=FALSE){
	if (is.null(up.to)){
		
		if (n>5){
			up.to=3
		}else
			up.to=2
	}
	E = simulateEXP( n, rep.no ,  up.to)
	G = sampleRndGraph(E, collapse=collapse)$G
	m = nrow(E)
	betas =	sapply(1:Y.no, function(i){	be = (rnorm(n, 0, sd = sqrt(1/b))); names(be) = colnames(G); be } )
	betas[betas==0]<- 1e-05 # add a pseudocount to zero beta entries - zero entries are not admissable
	
	epsilon = matrix( rnorm( m*Y.no, 0, sigma), ncol=Y.no)
	X = ( (E%*%G)>0)
	
	Y = X%*%betas + epsilon
	
	return( list( Y = Y, betas = betas, E = E, G=G))
}