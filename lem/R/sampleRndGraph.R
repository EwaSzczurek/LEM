sampleRndGraph = function(E, collapse=TRUE){	
	
	allone=TRUE
	while (allone){ # we do not want a graph which will collapse into one node after transitive closure
		G = sampleRndNetwork(colnames(E))
		diag(G)=1
		allone=all(G==1)
	}
	if (collapse){
		GE	<-  collapse.cycles(G=G, E = E )
	}else{
		GE <- list(G  = G, E = E)
	}
	
	GE

}
