# Create Models in WinBUGS 
#
# Author: Earl Duncan
# Created: 10/11/2017
#==========================================================================

# Load packages
library(rgdal)				# For readOGR()
library(ggplot2)			# For fortify(), ggplot()
library(spdep)				# For poly2nb() and nb2WB()
library(R2WinBUGS)		    # For bugs()
library(gridExtra)		    # For grid.arrange()
library(scales)				# For rescale() (for plotting risk surface)

# Set filepaths
fp.wd <- getwd()    # Change as required
fp.wd <- "D:/11 BotB 2017/Workshop/R Files/1 Monday"
setwd(fp.wd)
fp.UDFs <- paste0(fp.wd,"/UDFs.R")
fp.data <- paste0(fp.wd,"/Data.csv")
fp.shp <- paste0(fp.wd,"/Scotland Shapefile.shp")
fp.plots <- paste0(fp.wd, "/Create Plots.R")
fp.plots.2 <- paste0(fp.wd, "/Create Plots (PPC).R")

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
n.thin <- 1

# Choose model
Models <- c("ICAR", "BYM", "Leroux")
Model <- Models[2]

# Initial values for MCMC chain
inits <- function(){list(
    alpha = 1,
    beta = 0,
    S = rnorm(N, 0, 1),
    U = 0.5,
    tau.S = 0.5,
    tau.U = 0.5
)}

# Paramaters to monitor
params.monitor <- c("mu", "alpha", "beta","S","U", "sigma.S", "sigma.U")

# Estimate parameters using WinBUGS
if(Model == "ICAR"){
    # Add code here
}else if(Model == "BYM"){
    MCMC <- bugs(
        data = c(list(y = y, E = E, x = x, N = N), adj.data),
        inits = inits,
        parameters.to.save = params.monitor,
        model.file = paste("BUGS Code for",Model, "Model.bug"),
        n.chains = 1,
        n.iter = M.burnin + (M * n.thin),
        n.burnin = M.burnin,
        n.thin = n.thin,
        debug = TRUE,
        bugs.directory = "c:/Program Files/WinBUGS14/"
    )
}else if(Model == "Leroux"){
    # Add code here
}


# Save estimates to disk (optional)
##	save(MCMC, file = paste0("MCMC Output (WinBUGS), ", Model, ".RData"))


#==========================================================================
# Posterior summary plots
#==========================================================================

# Modify code from CARBayes plots (should be similar)

# EOF

