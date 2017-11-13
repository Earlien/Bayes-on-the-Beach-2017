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



