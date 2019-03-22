simulateEXP <- function( n, rep.no, up.to=3){
	
	if (!up.to%in%c(2,3)){
		cat("up.to can either be 2 or 3\n")
		stop()
	}
	
	# set the number of experiments
	m = n+ choose(n, 2)
	if (up.to >2){
		m = m+ choose(n, 3)
	}
	m = m*rep.no
	
	E = matrix(0, ncol = n, nrow= m)
	colnames(E) = giveNames(n)
	e.double = rep.no*n + 1
	
	for( i in seq(1, n-1)){
		for (j in seq(i+1, n)){
			for (r in seq(1, rep.no)){
				E[e.double, i ] = 1
				E[e.double, j ] = 1
				e.double = e.double+1
			}
			
		}
	}
	
	e.single = 1
	for( i in seq(1, n) ){
		for ( r.s in seq( 1, rep.no ) ){
			E[e.single, i] = 1
			e.single = e.single+1
		}
	}

	if (up.to > 2){	
		e.triple = (n + choose(n, 2))*rep.no+ 1
	
		for( i in seq(1, n-2)){
			for (j in seq(i+1, n-1)){
				for( l in seq(j+1, n)){
					for (r in seq(1, rep.no)){
						E[e.triple, i ] = 1
						E[e.triple, j ] = 1
						E[e.triple, l ] = 1 
						e.triple = e.triple+1
					}
				}	
			}
		}
	}
			
	E
}

