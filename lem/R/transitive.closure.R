## This code is copied from package nem (Author: Florian Markowetz), and slightly modified
transitive.closure <- function(g){

    if (!(class(g)%in%c("matrix"))) stop("Input must be an adjacency matrix")
    g <- as(g, "matrix")
    

	n <- ncol(g)
	matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
		while(i < pow) {
			z <- z %*% x
			y <- y+z
			i <- i+1
		}
		return(y)
	}
	
	h <- matExpIterativ(g,n)
	h <- (h>0)*1   
	dimnames(h) <- dimnames(g)
	diag(h) <- rep(1,n)
       
    return(h)
}
