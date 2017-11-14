
# ===================================================================
# Posterior predictive performance plots
# ===================================================================

message("Rendering posterior predictive performance plots...")
flush.console()

# ---------------

# y fitted values and residuals

Append <- data.frame(
	id = 1:N, 
	val = c(apply(pars$fitted, 2, mean),
		   apply(pars$fitted, 2, mean) - y
	),
	var = c(sapply(c("Fitted y", "Residuals"), rep, N))
)

dat <- map.df
dat$id <- as.numeric(dat$id) + 1		# Make shapefile dataframe SLA ID numeric
dat <- merge(dat, Append, by = "id")	# Merge values with shapefile dataframe

gg.y <- ggplot(dat, aes(x = long, y = lat, group = group)) +
	geom_polygon(color = "black", aes(fill = val), size = 0.2) +
	scale_fill_continuous("", low = "white", high = "blue") +
	facet_grid(. ~ var) +
	theme.map


# -------------------------------------------------------------------

# Plot of observed data with prediction bands
dat <- pars$fitted
dat <- data.frame(
	id = original.order,
	mean = apply(dat, 2, mean), 
	CI.lower = apply(dat, 2, quantile, probs = 0.025), 
	CI.upper = apply(dat, 2, quantile, probs = 0.975),
	y = y
)
Colours <- ifelse((dat$CI.upper - dat$CI.lower) > 10, 2, 1)
Colours.2 <- ifelse((dat$CI.upper - dat$CI.lower) > 20, 3, 0)
Colours <- apply(cbind(Colours, Colours.2), 1, max)
dat$col <- Colours

gg.pred <- ggplot(dat, aes(x = id, y = mean)) + 
	geom_linerange(aes(ymin = CI.lower, ymax = CI.upper, color = factor(col)), 
		alpha = 0.5, size = 0.51) +
	scale_color_manual(guide = FALSE, values = c("cyan", "#cc9900", "red")) +
	geom_point(data = dat, aes(x = id, y = y), alpha = 0.6, shape = 16) +
	scale_y_continuous("y") +
	ggtitle("Observed values with 95% CI band for fitted values")


message("Plots Complete")

# EOF

