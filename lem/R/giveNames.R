add0s <- function( ch, maxlen){
	lendiff = maxlen - nchar(ch)
	if (lendiff > 0){
		toapp = paste(rep(0, lendiff), collapse="")
		ch = paste(toapp, ch, sep="")
	}

	ch
}

## this function generates such names which are lexicographically sorted by their number
giveNames<-function( n ){
	name <- as.character(1:n)
		## now we will append 0s in front
		maxlen = max(sapply(name, nchar))
		name = sapply(name, add0s, maxlen = maxlen)
		name = paste("g", name, sep="")
		name
}
