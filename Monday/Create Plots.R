
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

message("... complete")

# EOF
