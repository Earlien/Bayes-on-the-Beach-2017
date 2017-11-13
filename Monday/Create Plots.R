
# ===================================================================
# Recreate/modify parameters
# ===================================================================

pars <- MCMC$samples

# names(pars)		# BEFORE
# ICAR and Leroux: 	"beta"   "phi"    "tau2"   "rho"    "fitted" "Y"  
# BYM: 				"beta"   "psi"    "tau2"   "sigma2" "fitted" "Y" 

# Renaming phi or psi as R
names(pars)[which(names(pars) %in% c("phi", "psi"))] <- "R"

# Renaming sigma2 as sigma.sq.u (BYM only)
names(pars)[which(names(pars) == "sigma2")] <- "sigma.sq.u"

# Renaming tau2 as sigma.sq.s
names(pars)[which(names(pars) == "tau2")] <- "sigma.sq.s"

# Splitting beta and into beta and alpha
##pars$beta <- pars$beta[,2]
pars$alpha <- pars$beta[,1]
pars$beta <- pars$beta[,2]

# Recreate mu
if(exists("x")){
	pars$mu <- pars$alpha + outer(pars$beta, x) + pars$R
}else{
	pars$mu <- pars$alpha + pars$S
}

# Remove rho for ICAR
if(Model == "ICAR"){
	pars$rho <- NULL
}

# names(pars)	# AFTER
# ICAR : 		"beta"     "R"    "sigma.sq.s"              "fitted"     "Y"    "alpha"    "mu"  
# BYM: 			"beta"     "R"    "sigma.sq.s" "sigma.sq.u" "fitted"     "Y"    "alpha"    "mu"  
# Leroux: 		"beta"     "R"    "sigma.sq.s" "rho"        "fitted"     "Y"    "alpha"    "mu"   

# Number of parameters to be plotted
par.num <- length(pars) - 2

# ===================================================================
# Posterior summary plots
# ===================================================================

message("Rendering posterior summary plots ...")
flush.console()

gg.post <- vector("list", par.num)


# Plot for mu
dat <- pars$mu
mean <- apply(dat, 2, mean)
CI.lower <- apply(dat, 2, quantile, probs = 0.025)
CI.upper <- apply(dat, 2, quantile, probs = 0.975)
dat <- data.frame(id = original.order, mean, CI.lower, CI.upper)
gg.post[[1]] <- ggplot(dat, aes(x = id, y = mean)) + 
	geom_pointrange(aes(ymin = CI.lower, ymax = CI.upper), 
					alpha = 0.3, fatten = 0.1) +
	geom_point(size = 1) +
	scale_x_continuous("Area") +
	scale_y_continuous(bquote("Posterior mean and 95% CI for "~mu))

# Plot for alpha
dat <- data.frame(val = as.numeric(pars$alpha))
gg.post[[2]] <- ggplot(dat, aes(x = val)) + 
	geom_density(fill = "black", alpha = 0.5, color = NA) +
	scale_x_continuous("") +
	scale_y_continuous(bquote("Posterior density for "~alpha))

# Plot for beta
dat <- data.frame(val = as.numeric(pars$beta))
gg.post[[3]] <- ggplot(dat, aes(x = val)) +
	geom_density(fill = "black", alpha = 0.5, color = NA) +
	scale_x_continuous("") +
	scale_y_continuous(bquote("Posterior density for "~beta))

# Plot for R
dat <- pars$R
mean <- apply(dat, 2, mean)
CI.lower <- apply(dat, 2, quantile, probs = 0.025)
CI.upper <- apply(dat, 2, quantile, probs = 0.975)
dat <- data.frame(id = original.order, mean, CI.lower, CI.upper)
gg.post[[4]] <- ggplot(dat, aes(x = id, y = mean)) + 
	geom_pointrange(aes(ymin = CI.lower, ymax = CI.upper), 
					alpha = 0.3, fatten = 0.1) +
	geom_point(size = 1) +
	scale_x_continuous("Area") +
	scale_y_continuous(bquote("Posterior mean and 95% CI for "~R))

# Plot for sigma_squared_s
dat <- data.frame(val = as.numeric(pars$sigma.sq.s))
gg.post[[5]] <- ggplot(dat, aes(x = val)) + 
	geom_density(fill = "black", alpha = 0.5, color = NA) + 
	scale_x_continuous("") + 
	scale_y_continuous(bquote("Posterior density for "~sigma[s]^2))

if(Model == "BYM"){
	# Plot for sigma_squared_u
	dat <- data.frame(val = as.numeric(pars$sigma.sq.u))
	gg.post[[6]] <- ggplot(dat, aes(x = val)) + 
		geom_density(fill = "black", alpha = 0.5, color = NA) + 
		scale_x_continuous("") + 
		scale_y_continuous(bquote("Posterior density for "~sigma[u]^2))
}else if(Model == "Leroux"){
	# Plot for rho
	dat <- data.frame(val = as.numeric(pars$rho))
	gg.post[[6]] <- ggplot(dat, aes(x = val)) + 
		geom_density(fill = "black", alpha = 0.5, color = NA) + 
		scale_x_continuous("") + 
		scale_y_continuous(bquote("Posterior density for "~rho))
}





#==========================================================================
# Spatial plots (NEW)
#==========================================================================

message("Rendering posterior spatial plots ...")
flush.console()

gg.spat <- vector("list", 4)


Append <- data.frame(
    id = 1:N,
    R = apply(pars$R, 2, mean),
    Cov = mean(pars$beta) * x,
    log.SIR = apply(pars$mu, 2, mean),
    log.y = apply(log(pars$fitted), 2, mean)
)

# Create dataframe
dat <- map.df
dat$id <- as.numeric(dat$id) + 1
dat <- merge(dat, Append, by = "id")		# Merge values with shapefile dataframe

gg.base <-  ggplot(dat, aes(x = long, y = lat, group = group)) +
    theme_grey() + theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(vjust = -0.1, size = rel(0.8)),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom"
    ) +     # These themes just make the plot look cleaner/nicer
    ##coord_map() +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0.02, 0.02))  # Reduces white space around border

Colours.all <- c(
    "#003380", "#0047b3", "#005ce6", "#3399ff", "#66ccff", 
    "white", 
    "#ff9966", "#f45000", "#ff0000", "#af0303", "#831628"
)   # Diverging colour schemes away from zero

All <- Append[,2:5] # All variables in the diverging colour map
Lim <- max(abs(range(All)))
Values <- seq(-Lim, Lim, length = 11)
mid <- which(abs(Values) < 1e-10)   # Value corresponding to zero (white)
start <- ceiling(length(Colours.all)/2) - mid + 1
Colours <- Colours.all[start:(start + length(Values) - 1)]

# Rescale the values so that the maps have a consistent colour legend
# (this is very important for being able to compare across these plots)
Values.1 <- rescale(Values, from = range(Append$R))
Values.2 <- rescale(Values, from = range(Append$Cov))
Values.3 <- rescale(Values, from = range(Append$log.SIR))

gg.spat[[1]] <- gg.base + geom_polygon(aes(fill = R), color = "grey70", size = 0.1) +
    scale_fill_gradientn("", colours = Colours, values = Values.1,
                         breaks = seq(-10, 10, 1)) +
    ggtitle("R")

gg.spat[[2]] <- gg.base + geom_polygon(aes(fill = Cov), color = "grey70", size = 0.1) +
    scale_fill_gradientn("", colours = Colours, values = Values.2,
                         breaks = seq(-10, 10, 0.5)) +
    ggtitle(bquote("Covariate Effect"))

gg.spat[[3]] <- gg.base + geom_polygon(aes(fill = log.SIR), color = "grey70", size = 0.1) +
    scale_fill_gradientn("", colours = Colours, values = Values.3,
                         breaks = seq(-10, 10, 2)) +
    ggtitle("Log(SIR)")

gg.spat[[4]] <- gg.base + geom_polygon(aes(fill = log.y), color = "black", size = 0.1) +
    scale_fill_gradient("", low = "white", high = "black",
                        breaks = seq(0, 10, 1), limits = c(0, max(dat$log.y))) +
    ggtitle("Log(Observed Values)")

message("... complete")

# EOF
