# This script aims to review several software packages for implementing
# propensity score analysis in R.

# Load required libraries.
suppressMessages(library(causaldrf))
suppressMessages(library(cobalt))
suppressMessages(library(CBPS))
suppressMessages(library(dplyr))
suppressMessages(library(ipw))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(mgcv))
suppressMessages(library(Matching))
suppressMessages(library(MatchIt))
suppressMessages(library(twang))

# Set configuration for themes (non-data components of the ggplot).
ggtheme.config <- function(main, size, xt.angle){
  theme <- list(ggtitle(main),
                theme_bw(base_size = size),
                theme(plot.title = element_text(size = size),
                      axis.title.x = element_text(size = size, angle = 0),
                      axis.title.y = element_text(size = size, angle = 90),
                      axis.text.x = element_text(size = size, angle = xt.angle),
                      strip.text.x = element_blank(),
                      legend.justification = c(1, 0), legend.position = c(1, 0)))
  return(theme)
}

# Analysis of the Bangladesh neurodevelopment outcome data.
DATA.PATH <- stringr::str_c("/Users/cathylee/Documents/Harvard/Postdoc/",
                            "EarlySigns/earlyS.data/bgd_data.txt")
data <- read.table(DATA.PATH, header = TRUE)

# Neurodevelopmental outcome: ccs_z.
bgd <- data %>%
  dplyr::select(repro_id, ccs_z, mn_ln_c, momage_c, homescore_z, momIQ_z,
                age_c, momeducd, smokenv, sex) %>%
  na.omit

# Manganese exposure: mn_ln_z, mn_ln_q, mn_ln_b.
bgd$mn_ln_z <- with(bgd, (mn_ln_c - mean(mn_ln_c)) / sd(mn_ln_c))
bgd$mn_ln_q <- cut(bgd$mn_ln_z, quantile(bgd$mn_ln_z, seq(0, 1, length = 10)),
                   include.lowest = TRUE)
bgd$mn_ln_b <- as.logical(bgd$mn_ln_z > 0)

# Covariates: momage_c, homescore_z, momIQ_z, age_c, approxage, momeducd,
#             smokenv, sex.
bgd$momage_z <- with(bgd, (momage_c - mean(momage_c)) / sd(momage_c))
bgd$age_z <- with(bgd, (age_c - mean(age_c)) / sd(age_c))
bgd$sex <- bgd$sex - 1
bgd$momeducd <- factor(bgd$momeducd)
bgd$smokenv <- factor(bgd$smokenv)
bgd$sex <- factor(bgd$sex)

# Set up for response, exposure and covariate vector and matrices.
response <- bgd$ccs_z
cont.exposure <- bgd$mn_ln_z
cat.exposure <- bgd$mn_ln_b
covariates <- with(bgd, cbind(momage_z, homescore_z, momIQ_z, age_z, momeducd,
                              smokenv, sex))

# Set up for treatment formula and grid for plotting.
X.name <- colnames(covariates)
cont.formula <- as.formula(paste("mn_ln_z ~ ", paste(X.name, collapse = "+")))
cat.formula <- as.formula(paste("mn_ln_q ~ ", paste(X.name, collapse = "+")))
bin.formula <- as.formula(paste("mn_ln_b ~ ", paste(X.name, collapse = "+")))
exposure.grid <- quantile(cont.exposure, probs = seq(0, 1, by = 0.01))

# causaldrf: Tools for estimating causal dose response functions.
# Explore several estimators for estimating the average dose response function.
# The additive spline estimator from Bia et al. (2014).
# 1. Fit a treatment model to estimate the GPS.
# 2. Create additive spline bases values for the treatment and GPS.
# 3. Regress outcome on the treatment, GPS, treatment bases and GPS bases.
add_spl_estimate <- causaldrf::add_spl_est(
  Y = ccs_z,
  treat = mn_ln_z,
  treat_formula = mn_ln_z ~ momage_z + homescore_z + momIQ_z + age_z +
    momeducd + smokenv + sex,
  data = bgd,
  grid_val = exposure.grid,
  knot_num = 3,
  treat_mod = "Normal",
  link_function = "identity")

# The generalized additive model estimator.
# 1. Fit a treatment model to estimate the GPS.
# 2. Regress outcome on the treatment and GPS bases from the GPS fit.
gam_estimate <- causaldrf::gam_est(
  Y = ccs_z,
  treat = mn_ln_z,
  treat_formula = mn_ln_z ~ momage_z + homescore_z + momIQ_z + age_z +
    momeducd + smokenv + sex,
  data = bgd,
  grid_val = exposure.grid,
  treat_mod = "Normal",
  link_function = "identity")

# The Hirano-Imbens estimator.
# 1. Regress treatment on a set of covariates to estmate the GPS.
# 2. Regress outcome on the treatment and GPS.
hi_estimate <- causaldrf::hi_est(
  Y = ccs_z,
  treat = mn_ln_z,
  treat_formula = mn_ln_z ~ momage_z + homescore_z + momIQ_z + age_z +
    momeducd + smokenv + sex,
  outcome_formula = ccs_z ~ mn_ln_z + I(mn_ln_z ^ 2) + gps + I(gps ^ 2) +
    mn_ln_z * gps,
  data = bgd,
  grid_val = exposure.grid,
  treat_mod = "Normal",
  link_function = "identity")

# The BART estimator described in Hill (2011).
# Fit a rich outcome model on the treatment and covariates to create a flexible
# response surface.
bart_estimate <- causaldrf::bart_est(
  Y = ccs_z,
  treat = mn_ln_z,
  outcome_formula = ccs_z ~ mn_ln_z + momage_z + homescore_z + momIQ_z + age_z +
    momeducd + smokenv + sex,
  data = bgd,
  grid_val = exposure.grid)

# The estimator described in Flores et al. (2012).
# This is a local linear regression of the outcome on the treatment with a
# weighted kernel. The weighted kernel is weighted by the reciprocal of the GPS
# values.
iw_estimate <- causaldrf::iw_est(
  Y = ccs_z,
  treat = mn_ln_z,
  treat_formula = mn_ln_z ~ momage_z + homescore_z + momIQ_z + age_z +
    momeducd + smokenv + sex,
  data = bgd,
  grid_val = exposure.grid,
  bandw = 2 * bw.SJ(bgd$mn_ln_z),
  treat_mod = "Normal")

# Plot of average dose response function curves.
ylim <- c(-1.5, 1)
par(mar = c(6, 6, 2, 2))
plot(exposure.grid, add_spl_estimate$param, type = "l",
     ylab = "estimate", xlab = "log(Mn concentration)",
     xlim = c(-2, 4), ylim = ylim,
     lwd = 2, cex.axis = 2, cex.lab = 2)
lines(exposure.grid, gam_estimate$param, col = "red", ylim = ylim, lwd = 2)
lines(exposure.grid, hi_estimate$param, col = "blue", ylim = ylim, lwd = 2)
lines(exposure.grid, bart_estimate$param, col = "orange", ylim = ylim, lwd = 2)
lines(exposure.grid, iw_estimate$param, col = "purple", ylim = ylim, lwd = 2)
legend("topleft",
       legend = c("Additive spline", "GAM", "Hirano-Imbens", "BART", "IPWT"),
       col = c("black", "red", "blue", "orange", "purple"),
       lty = 1, lwd = 2, bty = "n", cex = 2, ncol=2)

# Matching: Multivariate and propensity score matching software with automated
# balance optimization. Only works for a binary treatment.

# Estimate our first propensity score model.
ps <- glm(cat.formula, family = binomial, data = bgd)

# Perform one-to-one matching with replacement using our preliminary propensity
# score model where the estimand is the average treatment effect on the treated.
one_to_one <- Matching::Match(Y = bgd$ccs_z, Tr = bgd$mn_ln_b, X = ps$fitted)
Matching::MatchBalance(bin.formula, match.out = one_to_one, nboots = 1000,
                       data = bgd)

# Assessment of overlap/common support: restrict data set to observations that
# overlap and have a common support using the formula from Bia et al. (2014).
overlap_list <- causaldrf::overlap_fun(
  Y = ccs_z,
  treat = mn_ln_z,
  treat_formula = mn_ln_z ~ momage_z + homescore_z + momIQ_z + age_z +
    momeducd + smokenv + sex,
  data = bgd,
  n_class = 3,
  treat_mod = "Normal")
overlapped_data <- overlap_list$overlap_dataset
summary(overlapped_data)

# cobalt: Covariate balance tables and plots.
# The balance statistic used is the Pearson correlation between each covariate
# and the treatment variable. A threshold for balance on correlations can be
# specified using r.threshold; though there are no specific guidelines for this
# value, we choose 0.1 as indicating balance. Because the goal is complete
# independence between treatment and covariates, not simply the absence of a
# linear correlation between treatment and covariates, including interactions
# and squared terms through the use of int = TRUE is recommended.
fitted.CBPS <- CBPS(formula = cont.formula, ATT = 0, data = bgd)
bal.tab(fitted.CBPS, un = TRUE, r.threshold = 0.1, int = TRUE)

# If either fit line is not close to flat (i.e., lying on top of the reference
# line), there may be some remaining dependence between treatment and the
# covariate.
bal.plot(fitted.CBPS, "age_c", un = TRUE)
bal.plot(fitted.CBPS, "age_c")

# For binary/categorical covariates, bal.plot() displays a density plot of the
# treatment variable in each category. If treatment and the covariate are
# independent, the densities for each category should overlap with each other.
bal.plot(fitted.CBPS, "sex", un = TRUE)
bal.plot(fitted.CBPS, "sex")

# Summarize balance in a love plot.
love.plot(bal.tab(CBPS.fit), threshold = 0.1, abs = TRUE,
          var.order = "unadjusted")

# twang: The toolkit for weighting and analysis of nonequivalent groups. The
# propensity score estimation is centered on using boosted logistic regression
# as implemented in gbm.
#fitted.cont.ps <- gbm(cont.formula, data = bgd, distribution = "gaussian",
#                      n.trees = 3000, n.cores = 2)
fitted.ps <- twang::ps(formula = bin.formula, data = bgd, n.trees = 3000,
                       estimand = "ATE", verbose = FALSE)
ps.weights <- twang::get.weights(fitted.ps)

# Perform diagnostic check to obtain an optimal number of "n.trees".
plot(fitted.ps)

# Compute the relative influence of each variable for estimating the
# probability of treatment assignment.
summary(fitted.ps$gbm.obj, n.trees = fitted.ps$desc$ks.mean.ATE$n.trees,
        plot = TRUE)

# Generate useful diagnostic plots from the propensity score object.
plot(fitted.ps, plots = "boxplot")
plot(fitted.ps, plots = "histogram")

data <- with(bgd, data.frame(y = ccs_z, weights = ipw.weights, z = mn_ln_b))
ggplot(data, aes(y, weights = ipw.weights, fill = z)) +
  geom_histogram()

# Repeatedly use the ps() function and compare each treatment to the pooled
# sample of other treatments.
fitted.mnps <- twang::mnps(formula = cat.formula, data = bgd, n.trees = 3000,
                           estimand = "ATE", verbose = FALSE)
mnps.weights <- get.weights(fitted.mnps)

# Perform diagnostic check to obtain an optimal number of "n.trees".
plot(fitted.mnps)

# Generate useful diagnostic plots from the propensity score object.
plot(fitted.mnps, plots = "boxplot")
plot(fitted.mnps, plots = 3)

# ipw: Inverse probability weights (point treatment).
# The weights are computed using the ratio of predicted densities. The numerator
# is the marginal density function of Z and the denominator is the condiitonal
# density function of Z given X.
fitted.ipw <- ipwpoint(exposure = mn_ln_z,
                       family = "gaussian",
                       numerator = ~ 1,
                       denominator = ~ momage_z + homescore_z + momIQ_z +
                         age_z + momeducd + smokenv + sex,
                       data = bgd)
ipw.weights <- fitted.ipw$ipw.weights
ipwplot(ipw.weights)

# CBPS: Covariate balancing propensity score estimation.
# Parametric.
fitted.CBPS <- CBPS(formula = cont.formula, ATT = 0, data = bgd)
CBPS.weights <- fitted.CBPS$weights
summary(fitted.CBPS)
plot(fitted.CBPS, boxplot = TRUE)

# Non-parametric.
fitted.npCBPS <- npCBPS(formula = cont.formula, data = bgd,
                        corprior = 0.1/nrow(bgd))
npCBPS.weights <- fitted.npCBPS$weights
summary(fitted.npCBPS)
plot(fitted.npCBPS, boxplot = TRUE)

# Compare response-exposure functions across GAM models.
#bgd$gps <- ps.weights
#bgd$gps <- mnps.weights
#bgd$gps <- ipw.weights
#bgd$gps <- CBPS.weights
bgd$gps <- npCBPS.weights
gam.standard <- gam(ccs_z ~ s(mn_ln_z) + momage_z + homescore_z + momIQ_z +
                      age_z + momeducd + smokenv + sex, data = bgd)
gam.IPWT <- gam(ccs_z ~ s(mn_ln_z), weights = gps, data = bgd)
gam.regression <- gam(ccs_z ~ s(mn_ln_z) + gps, data = bgd)

plothFuncViaGam <- function(gam.fit, method, title) {
  mn.sequence <- seq(min(bgd$mn_ln_z), max(bgd$mn_ln_z), length = 100)
  if (method == "standard") {
    test.data <- data.frame(mn_ln_z = mn.sequence,
                            momage_z = mean(gam.fit$model$momage_z),
                            homescore_z = mean(gam.fit$model$homescore_z),
                            momIQ_z = mean(gam.fit$model$momIQ_z),
                            age_z = mean(gam.fit$model$age_z),
                            momeducd = 1,
                            smokenv = 1,
                            sex = 1)
  } else if (method == "IPWT") {
    test.data <- data.frame(mn_ln_z = mn.sequence)
  } else if (method == "regression") {
    test.data <- data.frame(mn_ln_z = mn.sequence,
                            gps = mean(gam.fit$model$gps))
  }
  fit <- predict(gam.fit, newdata = test.data, type = "response", se = TRUE)
  predicts <- data.frame(test.data, fit)
  g <- ggplot(predicts, aes_string(x = "mn_ln_z", y = "fit")) +
    geom_smooth(aes(ymin =  predicts$fit - 1.96 * (predicts$se.fit),
                    ymax = predicts$fit + 1.96 * predicts$se.fit),
                fill = "gray80", size = 1, stat = "identity") +
    scale_x_continuous(limits = c(-2, 4)) +
    scale_y_continuous(limits = c(-2, 1)) +
    xlab("log (Mn concentration)") + ylab("estimate") +
    ggtheme.config(title, 25, 0)
  return(g)
}

g1 <- plothFuncViaGam(gam.standard, "standard", "Standard")
g2 <- plothFuncViaGam(gam.IPWT, "IPWT", "GPS-weighting")
g3 <- plothFuncViaGam(gam.regression, "regression", "GPS-regression")
grid.arrange(g1, g2, g3, ncol = 3)
