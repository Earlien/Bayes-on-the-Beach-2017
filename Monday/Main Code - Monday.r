# Create Models in CARBayes
#
# Author: Earl Duncan
# Created: 10/11/2017
#==========================================================================

# Load packages
library(rgdal)				# For readOGR()
library(ggplot2)			# For fortify(), ggplot()
library(spdep)				# For poly2nb() and nb2WB()
library(CARBayes)			# For S.CARbym(), S.CARleroux(), etc.
library(gridExtra)		    # For grid.arrange()
library(scales)				# (Might use in Tuesday session)

# Set filepaths
fp.wd <- getwd()    # Change as required
setwd(fp.wd)
fp.UDFs <- paste0(fp.wd,"/UDFs.R")
fp.data <- paste0(fp.wd,"/Data.csv")
fp.shp <- paste0(fp.wd,"/Scotland Shapefile.shp")
fp.plots <- paste0(fp.wd, "/Create Plots.R")

# Load user-defined functions
source(fp.UDFs)

# Read in data
data <- read.csv(fp.data)
map <- readOGR(fp.shp)			# Shapefile
map.df <- fortify(map)			# Shapefile as data.frame (for plots)

# Re-order records to match order in shapefile
original.order <- match(map$id, data$Area)
data <- data[original.order,]

# Set constants
y <- data$Obs			# Observed values
E <- data$Exp		 	# Expected values
x <- data$Cov			# Covariate
N <- length(y)			# Number of areas

# Add small perturbation to E if zero (otherwise offset in model will error)
E[which(E == 0)] <- 0.1


#==========================================================================
# Create spatial weights matrix
#==========================================================================

# Created neighbourhood/adjacency matrix
adj.data <- nb2WB(poly2nb(map))

# Create first-order, binary spatial weights matrix
W <- adj2W(adj.data)

# Check neighbourhood/adjacency information is correct
W.test <- W
diag(W.test) <- 2   # Highlight selected area in plot

Area <- 1

Append <- cbind(
    area = rep(Area, N),
    id = 1:N, 
    "neigh" = c(W.test[,Area])
)
dat <- map.df
dat$id <- as.numeric(dat$id) + 1			# Make shapefile dataframe SLA ID numeric
dat <- merge(dat, Append, by = "id")		# Merge values with shapefile dataframe

ggplot(dat, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = factor(neigh)), color = "black", size = 0.2) +
    scale_fill_manual(values = c("white", "orange", "red"), guide = FALSE) +
    theme_bw()

# NB: the Scotland adjacency info has errors.
# Please import corrected version from file
load("Scotland Adjacency Info.RData")
W <- adj2W(adj.data)

# You can check the errors have been corrected be re-running above code.


#==========================================================================
# Create model
#==========================================================================

# Set MCMC parameters
M.burnin <- 10000
M <- 5000
n.thin <- 5

# Choose model
Models <- c("ICAR", "BYM", "Leroux")
Model <- Models[1]

# Estimate parameters using CARBayes
if(Model == "ICAR"){
	MCMC <- S.CARleroux(
		formula = y ~ offset(log(E)) + x, 
		data = data.frame(y, E, x),
		family = "poisson",
		prior.mean.beta = # Add value(s) here,
		prior.var.beta = # Add value(s) here,
		prior.tau2 = # Add value(s) here,
		fix.rho = # Add value(s) here,
		rho = # Add value(s) here,
		W = W,
		burnin = M.burnin,
		n.sample = M.burnin + (M * n.thin),
		thin = n.thin
	)
}else if(Model == "BYM"){
	# Add code here
}else if(Model == "Leroux"){
    # Add code here
}

# Save estimates to disk (saves us re-running model again  tomorrow)
save(MCMC, file = paste0("MCMC Output, ", Model, ".RData"))


#==========================================================================
# Posterior summary plots
#==========================================================================

source(fp.plots)

# Export
DPI <- 200			# Resolution
png(filename = "Posterior Distributions.png", 
	width = 700*DPI/72, height = 400*DPI/72, pointsize = 15, res = DPI)
if(Model == "ICAR"){
	grid.arrange(gg.post[[1]], gg.post[[2]], gg.post[[3]], 
				 gg.post[[4]], gg.post[[5]], ncol = 3)
}else{
	grid.arrange(gg.post[[1]], gg.post[[2]], gg.post[[3]], 
				 gg.post[[4]], gg.post[[5]], gg.post[[6]],ncol = 3)
}
dev.off()



# Add some code to the "Create Plots" script to create spatial maps
# for log-relative risk, spatial random effect, covariate effect, etc.



png(filename = "Risk Surface Deconstruction.png", 
	width = 1200*DPI/72, height = 600*DPI/72, pointsize = 15, res = DPI)
grid.arrange(gg.spat[[1]], gg.spat[[2]], gg.spat[[3]], gg.spat[[4]], nrow = 1)
dev.off()


#==========================================================================
# Model Checking
#==========================================================================

# To do tomorrow

# If you get this far on Monday, see "Additional exercises" in powerpoint
# slides.


# EOF

