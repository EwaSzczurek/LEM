
converge <- function(pars, I, XtX, XtY , X, Y, E, n, m){
  b=as.numeric( pars[1] )
  c=as.numeric( pars[2] )
  l.old = 0
  l = -Inf
  d = Inf
  iter = 0
  ### estimation of computing mu and V and hyperparameters b and c via iteration until convergence
  while( (d >= 0.0001) & (iter<= 1000)){
  	iter = iter+1
  	V.1 = b*I + c*XtX	
  	out <- tryCatch(
        {
            suppressWarnings( solve(V.1)  )
        },
        error=function(cond) {
            return(NULL)
        },
        warning=function(cond) {
            return(NULL)
        },
        finally={  
        }
    )
    if (! is.null(out) & (! all(is.na(out)))){
    	V = out
  		mu = ( c*V%*%XtY )
	  	mu.scalar = as.numeric( t(mu)%*%mu )
  		mu.sq.err  = as.numeric( sum(( Y - X%*%mu)^2) )
  		E.mu = as.numeric( c/2*mu.sq.err + b/2*mu.scalar )
  		l.old = l 
  		l = n/2*log(b) + m/2*log(c)  - E.mu + 1/2*log( det( V )) - m/2 * log(2*pi)
  		b = as.numeric( n/mu.scalar )
  		c = as.numeric( m/mu.sq.err )
  		d = abs(l - l.old)
  		if (l == -Inf)
  			d = 0	
	}else{
		l = -Inf
		d = 0
	}
  }
   unconverged = is.na(l)
  
   if (!unconverged){
   	unconverged <- any(is.na(mu))
   }
  
   if (unconverged)
   	l = -Inf
  
  #adding pseudocounts 	    
   if (!unconverged){
  	 if (any(mu == 0)){ 
	  	 mu[mu==0] <- 0.00001
  		 mu.scalar = as.numeric( t(mu)%*%mu )
  		 mu.sq.err  = as.numeric( sum(( Y - X%*%mu)^2) )
  		 b = as.numeric( n/mu.scalar )
  		 c = as.numeric( m/mu.sq.err )
  		 E.mu = as.numeric( c/2*mu.sq.err + b/2*mu.scalar )
  		 l = n/2*log(b) + m/2*log(c)  - E.mu + 1/2*log( det( V )) - m/2 * log(2*pi)
  	 }
   }
  
  l=as.numeric(l)
  c( b=b, c=c, l=l, mu=as.numeric( mu )) 
}

inference.bayes.linear.reg <- function(G, Y, E, verbose=FALSE) {
  X = ( ( E%*%G ) > 0 )
  X[X]<-1
  n = ncol(G)
  m = nrow(E)
  
  I = matrix( 0, nrow = n, ncol = n)
  diag( I ) = 1
  XtX = t(X)%*%X
  XtY = t(X)%*%Y
  
  b = c( 0.1, 0.5, 1, 2, 4, 1/var(Y))
  c = c( 0.1, 0.5, 1, 2, 4, 1000)
  pars=expand.grid(b,c)
  conv=apply(pars,1, converge, I, XtX,XtY,X, Y, E, n, m)
  liks =(conv["l",])
  if (!all(is.na(liks))){	
	  sel = which.max(liks)
	  res=converge( as.numeric( pars[sel,] ), I, XtX, XtY, X,Y, E, n, m)
  	  res=list(LL=res["l"], beta=res[4:length(res)], b=res["b"], c=res["c"])
  }else{
  	  res=list(LL=-Inf, beta=NULL, b=NULL, c=NULL)
  }
  
  res
}

inference.linear.reg <- function(G, Y, E, verbose=FALSE) {
	
	X = ( ( E%*%G ) > 0 )
	X[X]<-1
	
	# the graph can be collapsed to one node only due to being a clique. in such a case, fit a model only with the intercept 
	if (ncol(G)==1){
		model=stats::lm( formula = Y ~ 1)
	}else{
	# fit a simple linear model to Y. We know X has enough rows and has such values that $\beta$ is identifiable from Y
		model=stats::lm( formula = Y ~ 0+., data=data.frame(X))
	}
	ll = as.numeric(stats::logLik(model))
	beta = model$coefficients

    res=list(LL=ll, beta=beta, b=NULL, c=1/var(Y))
  	res
}

inference <- function(G, Y, E, parameter.estimation="linear.reg", verbose=FALSE) {
	
	if (parameter.estimation=="linear.reg")
		res = inference.linear.reg(G, Y, E, verbose) 
	else
		res=inference.bayes.linear.reg(G, Y, E, verbose)
   
  	res
}
