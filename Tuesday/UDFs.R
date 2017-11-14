#===========================================================================
# User-Defined Functions
#===========================================================================

# Get weight matrix from adjacency data
adj2W <- function(adj.data){
	neighbours <- adj.data$adj
	num <- adj.data$num
	W <- matrix(0, length(num), length(num))
	for(i in 1:length(num)){
		if(i == 1){
			start <- 1
		}else{
			start <- sum(num[1:(i-1)])+1
		}
		end <- sum(num[1:i])
		if(start <= end){
			W[i,neighbours[start:end]] <- 1
		}
	}
	if(!isSymmetric(W)){
		stop("W is not symmetric. Check adj.data.")
	}
	return(W)
}


#--------------------------------------------------------------------------
# Assessing Model Fit and Adequacy
#--------------------------------------------------------------------------

# Watanabe-Akaike (or Widely Applicable) Information Criterion (WAIC)
get.WAIC <- function(theta, y, L, type){
	# Pointwise evaluation of the Likelihood and log-likelihood
	if(is.list(theta) & length(theta) == 2){
		LogLike <- matrix(NA, nrow(theta[[1]]), length(y))
		for(i in 1:length(y)){
			LogLike[,i] <- L(y[i], theta[[1]][,i], theta[[2]][,i], log = TRUE)
		}
	}else if(is.atomic(theta)){
		LogLike <- matrix(NA, nrow(theta), length(y))
		for(i in 1:length(y)){
			LogLike[,i] <- L(y[i], theta[,i], log = TRUE)
		}
	}
	Like <- exp(LogLike)
	S.log <- apply(LogLike, 2, mean)
	S <- apply(Like, 2, mean)
	WAIC <- data.frame(WAIC = type, value = NA)
	for(t in 1:length(type)){
		switch(as.character(type[t]),
			   "1" = {
			   	WAIC[t,2] <- 2 * sum(log(S) - 2*S.log)
			   },
			   "2" = {
			   	Var <- apply(LogLike, 2, var)
			   	WAIC[t,2] <- 2 * sum(Var - log(S))
			   }
		)
	}
	return(WAIC)
}


# Morans.I
get.Morans.I <- function(y, W, normalise = TRUE){
	# Normalise weights if requested
	if(normalise == TRUE){
		W <- W/apply(W, 1, sum)
	}
	
	# Computes Moran's I statistic, and the expected value and variance
	# under a null hypothesis that of no spatial autocorrelation.
	N <- nrow(W)
	y.diff <- y - mean(y)
	
	# Observed value (test statistic)
	s <- sum(W)
	top <- sum(W * outer(y.diff, y.diff))
	bot <- sum((y.diff)^2)
	I <- N * top / (s * bot)
	
	# Expected value
	Exp.I <- -1 / (N - 1)
	
	# Variance
	s.sq <- s^2
	S1 <- 0.5 * sum((W + t(W))^2)
	S2 <- sum((apply(W, 1, sum) + apply(W, 2, sum))^2)
	S3 <- 1/N * sum(y.diff ^ 4) / ((1/N * sum(y.diff^2)) ^ 2)
	S4 <- (N ^ 2 - 3 * N + 3) * S1 - N * S2 + 3 * s.sq
	S5 <- (N ^ 2 - N) * S1 - 2 * N * S2 + 6 * s.sq
	
	##	Var.I <- (N * S4 - S3 * S5) / ((N - 1) * (N - 2) * (N - 3) * s.sq) - Exp.I^2
	Var.I <- (N * S4 - S3 * S5) / ((N - 1) * (N - 2) * (N - 3) * s.sq) - 1/ (N - 1)^2
	
	# Z scores
	Z.score.I <- (I - Exp.I) / sqrt(Var.I)
	
	# p-value for a two sided test is the probability of observing a
	# test statistic (z value) larger than the absolute value, times 2
	# to account for lower tail.
	p.1 <- pnorm(abs(Z.score.I), lower.tail = FALSE)
	p.2 <- p.1 * 2
	
	# Return values
	return(data.frame(Stat = I, Expected = Exp.I, Variance = Var.I, 
		z = Z.score.I, p.1.sided = p.1, p.2.sided = p.2, row.names = ""))
}


#--------------------------------------------------------------------------
# Mapping
#--------------------------------------------------------------------------

# Shorthand for creating clean map
theme.map <- theme_grey() + theme(
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank(),
	plot.title = element_text(vjust = -0.1, size = rel(0.8)),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	panel.background = element_blank(),
	strip.background = element_blank()
)


# EOF



