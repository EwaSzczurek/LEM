print.lem <- function(x,...) {

  # general
  cat("Object of class ",class(x),"\n")
  cat("\n")
  
  # slots
  cat("$graph: adjacency matrix for",ncol(x$graph),"genes\n") 
  print( as(x$graph, "matrix") ) 
  cat("Score: ", x$score, "\n")
  cat("Beta coefficients: ",x$beta,"\n")
  cat("Inference method:", x$inference,"\n")
  cat("Parameter estimation method:", x$parameter.estimation, "\n")  
}
