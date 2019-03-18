# this code is copied ans slightly modified from package nem (author: Florian Markowetz)
sampleRndNetwork = function(Sgenes, scaleFree=TRUE, gamma=2.5, maxOutDegree=length(Sgenes), maxInDegree=length(Sgenes), trans.close=TRUE, DAG=FALSE){		
	n = length(Sgenes)
	S = diag(n) # network of S-genes	
	maxOutDegree = min(n,maxOutDegree)	
	degprob = (0:maxOutDegree)^(-gamma)
	degprob[1] = 1
	degprob = degprob/sum(degprob)	
	for(i in 1:n){ # connect network randomly
		if(scaleFree)
			outdeg = sample(0:maxOutDegree,1,prob=degprob)# power law for out-degree => scale-free network
		else
			outdeg = sample(0:maxOutDegree,1)		
		if(outdeg > 0){
			if(!DAG)
				idx0 = which(S[i,] == 0)
			else
				idx0 = which(S[i,] == 0 & 1:n < i) # sample on lower triangle matrix
			if(length(idx0) > 0){
				idx = which(colSums(S[,idx0,drop=FALSE]) <= maxInDegree)
				if(length(idx) > 0){			
					idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)
					S[i,idx] = 1					
				}			
			}	
		}
	}						
	if(trans.close)
		S = transitive.closure(S)
	# we assume the diagonal is always 1 (the perturbed node is perturbed)				
	diag(S) = 1			
	colnames(S) = Sgenes
	rownames(S) = Sgenes
	S
}
