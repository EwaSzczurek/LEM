collapse.cycles<- function( G, E=NULL, beta=NULL ){
	cycle_found = TRUE
	diag(G) = 1

	while( cycle_found ){
		cycle_found=FALSE
		if ( ncol(G) > 1 ){
			for ( i in  (1: (ncol(G)-1) )){
				for ( j in ( (i+1) : ncol(G) )){
					
					if ( (G[i,j]==1) & (G[j,i]==1) ){
						
						#collapse
						node_i_name = colnames(G)[i]
						node_j_name = colnames(G)[j]
						new_node_name = paste( c( node_i_name, "_", node_j_name), collapse="")
						
						in.edges = as.numeric( ( G[, node_i_name] + G[, node_j_name]) > 0 )
						out.edges = G[i, ]
						# add new collapsed node
						G = cbind(in.edges, G)
						G = rbind(c(1, out.edges), G)
						colnames(G) = rownames(G) = c(new_node_name, colnames(G)[2:ncol(G)]) 

						## add the node to E 
						if (! is.null(E)){
							in.experiments <- as.numeric( ( E[, node_i_name] + E[, node_j_name]) > 0 )
							E = cbind(in.experiments, E)
							colnames(E) <- c(new_node_name, colnames(E)[2:ncol(E)])
						}
						
						## make the total contribution
						if (!is.null(beta)){
							beta = rbind(sum(beta[i,], beta[j,]), beta )
							rownames(beta ) = c(new_node_name, rownames( beta)[2:nrow(beta)] )
						}

						# remove the nodes i and j
						G = G[ -( which(rownames(G) == node_i_name)), , drop=F ]
						G = G[ -( which(rownames(G) == node_j_name)), , drop=F ]
						G = G[ , -( which(colnames(G) == node_j_name)), drop=F ]
						G = G[ , -( which(colnames(G) == node_i_name)), drop=F ]

						if (! is.null(E)){
							E = E[,-( which(colnames(E) == node_i_name)), drop=F ]
							E = E[,-( which(colnames(E) == node_j_name)), drop=F ]
						}	
						
						if (!is.null(beta)){
							beta = beta[-( which(rownames(beta) == node_i_name)), ,drop=F]
							beta = beta[-( which(rownames(beta) == node_j_name)), ,drop=F]
						}
	
						cycle_found = TRUE
						break();
					} 					
				}
				if (cycle_found) break()
			}
		}
	}
	
	
	res = list( G = G, E = E, beta = beta)

	res
}

uncollapse.cycles <- function( G, E=NULL, beta=NULL ){
	
	for ( cn in  colnames(G) ){
		
		spl = strsplit( cn, "_")[[1]]
		l.spl = length(spl)
		
		if (l.spl >1)	{
			
			# mark the index of the node in terms of rows and columns
			r.index = which( rownames(G) == cn ) 
			c.index = which( colnames(G) == cn ) 
			
			# mark the in and out edges of the collapsed node
			in.edges =  G[, cn]
			out.edges = G[cn, ]

			# remove the collapsed node from G
			G = G[ -( r.index ), , drop=F ]
			G = G[ , -( c.index ), drop=F ]
			in.edges = in.edges[-r.index]
			out.edges = out.edges[-c.index]

			#add the uncollapsed node to G
			spl.m = matrix(1, nrow = l.spl, ncol = l.spl)
			out.m=matrix( rep(out.edges, l.spl ), nrow=l.spl, byrow = TRUE)
			out.m = cbind(spl.m, out.m)
			colnames(out.m) = c( spl, colnames(G) )
			rownames(out.m) = spl
			
			in.m = matrix( rep(in.edges, l.spl ), ncol=l.spl, byrow = FALSE)			
			in.m = cbind(in.m, G )
			colnames(in.m)=c(spl, colnames(G))
			rownames(in.m)=rownames(G)
			G = rbind(out.m, in.m)


			## add the uncolapsed nodes to E 
			if (! is.null(E) ){
				E.col = E[, cn]
				E.spl = matrix( rep(E.col, l.spl), ncol = l.spl, byrow = FALSE )
				colnames(E.spl ) = spl
				E = cbind(E.spl, E)
				# remove the collapsed node
				E = E[,-( which(colnames(E) == cn) ), drop=F ]
	
			}
						
			## make the split contribution
			if (!is.null(beta)){
				beta.spl = matrix( rep( beta[ cn,]/l.spl, l.spl), nrow = l.spl, byrow = TRUE )
				rownames(beta.spl)=spl
				colnames(beta.spl)=colnames(beta)
				beta = rbind(beta.spl, beta )
				beta = beta[-( which(rownames(beta) == cn)), ,drop=F]
			}
		}
	}
	

	res = list( G = G, E = E, beta = beta)

	res
}

makeGEs <- function(Gs, E, verbose=FALSE, collapse=TRUE){
	### Remove the model which will be collapsed into 1 node
	alln1 <- sapply( Gs, function(G){ !all(G==1) })
	Gs <- Gs[alln1]
	if (collapse){
		Gs	<- lapply(Gs, collapse.cycles, E = E )
	}else{
		Gs = lapply( Gs, function(G){ list(G=G, E=E) })
	}
	Gs
}

enumerateGEs<-function(E, verbose=FALSE, collapse=TRUE){
	Gs<-enumerate.models(ncol(E), colnames(E), trans.close=TRUE, verbose=verbose)
	makeGEs(Gs=Gs, E=E, verbose=verbose, collapse=collapse)
}

