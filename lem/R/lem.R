## Main function 
## Arguments : 
## Y : the phenotype (response variable)
## E :  the experiments matrix
## Gs :  initial list of graphs (represented as adjacency matrices)
## inference :  the method with which the graph space is searched through
## parameter.estimation :  the method with which parameters are estimated
## verbose : whether the function shall print out messages
## collapse.res : whether the resulting graph should have cliques collapsed into a single node.

lem <- function(Y, E, Gs=NULL, inference="greedy", parameter.estimation="linear.reg", verbose=FALSE, greedy.init.no = 10, collapse.res = TRUE){

#------------------------------
# Sanity checks               
if (!(inference %in% c("greedy","search") )) 
    stop("\nlem> argument 'inference' is not valid\n")
if (!(parameter.estimation %in% c("linear.reg","bayes.linear.reg") )) 
    stop("\nlem> argument 'parameter.estimation' is not valid\n")

#------------------------------
# GREEDY                    
if(inference == "greedy"){
    result <- lem.greedy(Y, E, initial=Gs, parameter.estimation=parameter.estimation, greedy.init.no=greedy.init.no, verbose=verbose)
    
}

#------------------------------
# Score all given models/search exhaustively over all possible 
else if (inference == "search"){    
	if (is.null(Gs)){
		if (parameter.estimation == "bayes.linear.reg"){
			GEs <- enumerateGEs(E, verbose=verbose, collapse=FALSE)
		}else{ ## linear.reg
			GEs <- enumerateGEs(E, verbose=verbose, collapse=TRUE)
    		}
	}else{
		if (parameter.estimation == "bayes.linear.reg"){
	        GEs = makeGEs(Gs, E, verbose= verbose, collapse=FALSE)	
        }else{
	        GEs = makeGEs(Gs, E, verbose= verbose, collapse= TRUE)	
        	}
    }
	if (length(GEs)>0){
		result <- score(GEs, Y, parameter.estimation, verbose)
    }else{
        	result <-  list(graph = NULL, beta = NULL, score=-Inf, sigma = NULL, b = NULL, all.scores = c(-Inf))
    }		
}

result$parameter.estimation <- parameter.estimation
result$inference <- inference

if (collapse.res){
	##the final result needs to be collapsed
	res.col = collapse.cycles(G=as( result$graph, "matrix"), beta=result$beta )
	result$beta = res.col$beta
}
class(result)="lem"
#------------------------------
# OUTPUT                       
return(result)

}
