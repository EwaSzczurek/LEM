get.insertions = function(Phi){ ## All Phi transitively closed in the end.
    idx = which(Phi == 0)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible new edges
            Phinew = Phi
            Phinew[idx[i]] = 1
            Phinew = transitive.closure(Phinew) 
            models[[i]] <- Phinew
        }
    } 
    models       
}

delete.edge.transclose <- function( Phinew,  rowno , colno){
	Phinew[rowno, colno] = 0
	### update that the parents of gene rowno also remove edges to gene colno
	# get the parents of rowno
	parents = which( Phinew[,rowno] == 1 )
	# remove edges parents -- colno
	Phinew[parents, colno] = 0
	# still some parents could feed into colno via independent routes. For them we need to restore their edges to colno by transitive closure
	Phinew= transitive.closure(Phinew) 
	Phinew	
}

get.deletions = function(Phi){
    Phi = Phi - diag(ncol(Phi))
    idx = which(Phi == 1)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible edge deletions
            Phinew = Phi

            ### retreive the column at which we made the switch to 0
            colno = ceiling(idx[i]/nrow(Phi))
            ### retreive the row at which we made the switch to 0. This means that the parent rowno no longer has an edge to the gene colno
			rowno = idx[i] - (colno-1)*nrow(Phi)
			
			if (colno!=rowno){ # we do not make an attempt to remove the diagonal
				Phinew = delete.edge.transclose( Phinew,  rowno , colno)
				diag(Phinew) = 1
				if (!all(Phi==Phinew))
					models[[i]] <- Phinew
			}

        }
    } 
    models       
}

get.reversions = function(Phi){
    idx = which(Phi + t(Phi) == 1, arr.ind=TRUE)
    models = list()
    if(NROW(idx) > 0){
        for(i in 1:NROW(idx)){ # test all possible edge reversions
            Phinew = Phi
            Phinew[idx[i,1],idx[i,2]] = 0
            Phinew[idx[i,2],idx[i,1]] = 1
            diag(Phinew) = 1
    	        models[[i]] <- Phinew
        }
    } 
    models       
}

lem.greedy.sub <- function(Phi, Y, E,  parameter.estimation="linear.reg", verbose=TRUE){ 
    n <- ncol(E) 

    
    best.model = lem(Y, E, Gs=list(Phi), inference="search", parameter.estimation= parameter.estimation, verbose=verbose, collapse.res =FALSE)
    
    finished <- FALSE
    i = 1
    scores =best.model$all.scores
	sco0= best.model$score
    if (verbose){
    		cat("Greedy hillclimber for",n,"genes...\n\n")
    		print(Phi)
    		print(best.model)
    }
    while(!finished){
        models <- list()
#       propose new edges     
        models = get.insertions(Phi)
        models.del = get.deletions(Phi)
        models = c(models, models.del)
        models <- unique(models)
        if(verbose)
            cat(length(models), " local models to test ...\n")

        finished <- TRUE # this means no improving edge could be inserted
      
        if(length(models) > 0){		
            sconew <- lem(Y, E, Gs=models, inference="search", parameter.estimation= parameter.estimation, verbose=FALSE, collapse.res =FALSE) ## this returns one best model and a record of all acquired scores 
                
			if(sconew$score > sco0){ # this means an improving edge can be inserted or deleted
            	finished <- FALSE
            	
                	if(verbose)
                    	cat("step",i,"--> Edge added or removed\n")
                	sco0 <- sconew$score
                	winner = which.max(sconew$all.scores)
                	Phi <- models[[winner]]  ##store the best graph without collapsing
                	if(verbose)
                   	print(Phi)
                	best.model=sconew # store the best model 
                	scores=c(scores, sco0)         
            }
		}       
    		i = i+1
    }
    if (verbose){
    		cat("lem greedy made",i,"iterations\n\n")
    }
    res <- best.model
    res$all.scores=scores #saving the path of traversed scores
    if(verbose)
        cat("score of the best model = ", res$score,"\n\n")
    return(res)
}

lem.greedy <- function(Y, E, initial, parameter.estimation="linear.reg", greedy.init.no=10, verbose=TRUE){ 
	n <- ncol(E)
    Phi.d <- matrix(0, nrow=n, ncol=n)
    colnames(Phi.d)=colnames(E)
    rownames(Phi.d)=colnames(E)
    diag(Phi.d) <- 1    
            
    if (!is.null(initial)){
    		Phis <- initial
    	}else{
    		if (greedy.init.no > 0)
		    Phis = lapply( 1: greedy.init.no, function( i ){sampleRndGraph(E, collapse=FALSE)$G} )
    		else{
    			Phis = list()
			greedy.init.no  = 0
		}
    }
    Phis[[greedy.init.no+1]]=Phi.d
    Phis <- unique(Phis)
    evaluates = lapply( Phis, lem.greedy.sub,  Y, E,  parameter.estimation, verbose )
    scores = unlist( lapply(evaluates, function(el){el$score}) )
    best.model.ind = which.max(scores)
    best.model = evaluates[[best.model.ind]]
    return(best.model)
 }