
# You can add this code to the UDFs.R file

# Specify new weight matrix for nth-order neighbours
update.W <- function(weights, W){
	N <- nrow(W)
	W.list <- vector("list", length(weights))
	W.list[[1]] <- W * weights[1]
	for(k in 2:length(weights)){
		W.list[[k]] <- matrix(0, N, N)
		for(i in 1:N){
			neighbours.i <- which(W.list[[k-1]][i,] > 0) # (k-1)th order neighbours of area i
			if(length(neighbours.i) > 0){
				for(j in neighbours.i){
					neighbours.j <- W.list[[1]][j,]		# 1st order neighbours of area j
					neighbours.j[c(i, neighbours.i)] <- 0
					neighbours.j <- which(neighbours.j> 0)
					W.list[[k]][i, neighbours.j] <- weights[k]
				}
			}
		}
	}
	# Remove duplicates from layers
	for(k in 2:length(W.list)){
		for(kk in seq(1, (k-1), by = 1)){
			W.list[[k]][which(W.list[[kk]] != 0, arr.ind = TRUE)] <- 0
		}
	}
	W.new <- Reduce('+', W.list)

	return(W.new)
}


# EOF



