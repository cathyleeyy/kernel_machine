# This script aims to re-analyse the Bangladesh neurodevelopment outcome via
# bkmr() and gam() by adjusting confounders using generalized propensity score
# regression. The generalized propensity scores are calculated using three
# methods including multivariate linear regression, lm() or heavyLm(), and
# multivariate normal additive model, mvn in gam().

# Load required libraries.
suppressMessages(library(bkmr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(heavy))
suppressMessages(library(mgcv))
suppressMessages(library(mixtools))

# Set constant parameters.
is.runGAM <- FALSE
is.runBKMR <- TRUE
n.outcomes <- 3 # Number of neurodevelopment outcomes.
n.bkmr.iterations <- 100 # Number of iterations for bkmr(), should be 50000.
kTruncateValue <- 0.05 # Value for truncating propensity scores.
kDegreesOfFreedom <- 1e300 # Set it to large to be closer to normal.

# Set configuration for themes (non-data components of the plot).
ggtheme.config <- function(main, size, xt.angle){
  theme <- list(ggtitle(main),
                theme_bw(base_size = size),
                theme(plot.title = element_text(size = size),
                      axis.title.x = element_text(size = size, angle = 0),
                      axis.title.y = element_text(size = size, angle = 90),
                      axis.text.x = element_text(size = size, angle = xt.angle),
                      strip.text.x = element_blank(),
                      legend.justification = c(1, 0),
                      legend.position = c(1, 0)))
  return(theme)
}

# Read in the Bangladesh neurodevelopment outcome data.
DATA.PATH <- stringr::str_c("/Users/cathylee/Documents/Harvard/Postdoc/",
                            "EarlySigns/earlyS.data/bgd_data.txt")
data <- read.table(DATA.PATH, header = TRUE)
data$age_sq <- data$age_c ^ 2
bgd <- data %>% dplyr::select(ccs_z, as_ln_c, mn_ln_c, pb_ln_c, momage_c,
                              homescore_z, momIQ_z, age_c, age_sq, momeducd,
                              smokenv, sex, protein, clinic) %>%
  na.omit

bgd$momeducd <- factor(bgd$momeducd)
bgd$smokenv <- factor(bgd$smokenv)
bgd$sex <- factor(bgd$sex)
bgd$protein <- factor(bgd$protein)
bgd$clinic <- factor(bgd$clinic)

# Standardize exposure and continuous confounder variables if neccessary.
# bgd$as_ln_z <- with(bgd, (as_ln_c - mean(as_ln_c)) / sd(as_ln_c))
# bgd$mn_ln_z <- with(bgd, (mn_ln_c - mean(mn_ln_c)) / sd(mn_ln_c))
# bgd$pb_ln_z <- with(bgd, (pb_ln_c - mean(pb_ln_c)) / sd(pb_ln_c))
# bgd$momage_z <- with(bgd, (momage_c - mean(momage_c)) / sd(momage_c))
# bgd$age_z <- with(bgd, (age_c - mean(age_c)) / sd(age_c))

# Finalize the analysis data. From here onwards, we use unstandardized exposure
# and confounder variables for now.
bgd <- dplyr::select(bgd, ccs_z, as_ln_c, mn_ln_c, pb_ln_c, momage_c,
                     homescore_z, momIQ_z, age_c, age_sq, momeducd, smokenv,
                     sex, protein, clinic)

# Set up exposure matrices and formulaes for models.
exposures <- cbind(bgd$as_ln_c, bgd$mn_ln_c, bgd$pb_ln_c)

cov.part <- paste(stringr::str_c("momage_c + homescore_z + momIQ_z + age_c +",
                                 "age_sq + momeducd + smokenv + sex +",
                                 "protein + clinic"))
num.lm.formula <- as.formula("cbind(as_ln_c, mn_ln_c, pb_ln_c) ~ 1")
den.lm.formula <- as.formula(paste("cbind(as_ln_c, mn_ln_c, pb_ln_c) ~",
                                   cov.part))
num.gam.formula <- list(as_ln_c ~ 1, mn_ln_c ~ 1, pb_ln_c ~ 1)
den.gam.formula <- list(as.formula(paste("as_ln_c ~", cov.part)),
                        as.formula(paste("mn_ln_c ~", cov.part)),
                        as.formula(paste("pb_ln_c ~", cov.part)))

# Companion functions for use in the main functin ConstructStabilizedWeights().
TruncateWeights <- function(x, trunc){
  low <- quantile(x, 0 + trunc) # trunc
  upp <- quantile(x, 1 - trunc) # 1 - trunc
  message("Number of NaNs is ", sum(is.nan(x)))
  message("Number of x <= ", round(low, 2), " is ", sum(x[!is.nan(x)] <= low))
  message("Number of x > ", round(upp, 2), " is ", sum(x[!is.nan(x)] > upp))
	x[x <= low] <- low
	x[x > upp] <- upp
	return(x)
}

CalculatePropensityScores <- function(numerator, denominator) {
  scores <- vapply(1:nrow(bgd), function(i) {
    mixtools::dmvnorm(exposures[i, ], numerator[i, ], denominator)
  }, 0)
  return(scores)
}

ConstructStabilizedWeights <- function(method) {
  if (method == "lm") {
    num.fitted <- lm(num.lm.formula, data = bgd)
    den.fitted <- lm(den.lm.formula, data = bgd)
  } else if (method == "hlm") {
    num.fitted <- heavyLm(num.lm.formula, data = bgd,
                          family = Student(df = kDegreesOfFreedom))
    den.fitted <- heavyLm(den.lm.formula, data = bgd,
                          family = Student(df = kDegreesOfFreedom))
  } else if (method == "gam") {
    num.fitted <- gam(num.gam.formula, family = mvn(d = n.outcomes), data = bgd)
    den.fitted <- gam(den.gam.formula, family = mvn(d = n.outcomes), data = bgd)
  }

  num.mu <- num.fitted$fitted.values
  num.Sigma <- cov(residuals(num.fitted))
  num.ps <- CalculatePropensityScores(num.mu, num.Sigma)

  den.mu <- den.fitted$fitted.values
  den.Sigma <- cov(residuals(den.fitted))
  den.ps <- CalculatePropensityScores(den.mu, den.Sigma)

  weights <- num.ps / den.ps
  trunc.weights <- TruncateWeights(weights, kTruncateValue)

  return(list(weights = weights, trunc.weights = trunc.weights))
}

# Obtain stabilized weights from the three chosen multivariate models.
bgd$lm.weights <- ConstructStabilizedWeights("lm")$weights
bgd$lm.trunc.weights <- ConstructStabilizedWeights("lm")$trunc.weights
bgd$hlm.weights <- ConstructStabilizedWeights("hlm")$weights
bgd$hlm.trunc.weights <- ConstructStabilizedWeights("hlm")$trunc.weights
bgd$gam.weights <- ConstructStabilizedWeights("gam")$weights
bgd$gam.trunc.weights <- ConstructStabilizedWeights("gam")$trunc.weights

# Plot and compare the three stabailized weights.
par(mfrow = c(1, 3))
plot(density(bgd$lm.trunc.weights), lwd = 2)
plot(density(bgd$hlm.trunc.weights), lwd = 2)
plot(density(bgd$gam.trunc.weights), lwd = 2)

par(mfrow = c(1, 2))
plot(bgd$lm.trunc.weights, bgd$gam.trunc.weights)
plot(bgd$lm.trunc.weights, bgd$hlm.trunc.weights)

if (isTRUE(is.runGAM)) {
  # Begin analysis via gam().
  bgd$gps <- bgd$lm.trunc.weights # Use multivariate linear regression.

  outcome.part <- "ccs_z ~ s(as_ln_c, mn_ln_c, pb_ln_c) +"
  gam.standard <- gam(as.formula(paste(outcome.part, cov.part)), data = bgd)
  gam.gps <- gam(as.formula(paste(outcome.part, "s(gps, k = 4)")), data = bgd)

  plothFuncViaGam <- function(gam.fit, method, exposure.var) {
    is_Mn <- ifelse(exposure.var == "mn_ln_c", "y", "n")
    mn <- switch(is_Mn,
                 y = seq(min(bgd$mn_ln_c), max(bgd$mn_ln_c), length = 100),
                 n = mean(gam.fit$model$mn_ln_c))

    is_As <- ifelse(exposure.var == "as_ln_c", "y", "n")
    as <- switch(is_As,
                 y = seq(min(bgd$as_ln_c), max(bgd$as_ln_c), length = 100),
                 n = mean(gam.fit$model$as_ln_c))

    is_Pb <- ifelse(exposure.var == "pb_ln_c", "y", "n")
    pb <- switch(is_Pb,
                 y = seq(min(bgd$pb_ln_c), max(bgd$pb_ln_c), length = 100),
                 n = mean(gam.fit$model$pb_ln_c))

    if (method == "standard") {
      testdata <- data.frame(as_ln_c = as,
                             mn_ln_c = mn,
                             pb_ln_c = pb,
                             momage_c = mean(gam.fit$model$momage_c),
                             homescore_z = mean(gam.fit$model$homescore_z),
                             momIQ_z = mean(gam.fit$model$momIQ_z),
                             age_c = mean(gam.fit$model$age_c),
                             age_sq = mean(gam.fit$model$age_sq),
                             momeducd = 1,
                             smokenv = 1,
                             sex = 1,
                             protein = 1,
                             clinic = 1)
    } else if (method == "gps") {
      testdata <- data.frame(as_ln_c = as,
                             mn_ln_c = mn,
                             pb_ln_c = pb,
                             gps = mean(gam.fit$model$gps))
    }

    fit <- predict(gam.fit, newdata = testdata, type = "response", se = TRUE)
    predicts <- data.frame(testdata, fit)
    g <- ggplot(predicts, aes_string(x = exposure.var, y = "fit")) +
      geom_smooth(aes(ymin =  predicts$fit - 1.96 * (predicts$se.fit),
                      ymax = predicts$fit + 1.96 * predicts$se.fit),
                  fill = "gray80", size = 1, stat = "identity") +
      scale_x_continuous(limits = c(-2, 4)) +
      scale_y_continuous(limits = c(-2, 2)) +
      xlab(exposure.var) + ylab("estimate") +
      ggtheme.config("", 25, 0)
    return(g)
  }

  g1 <- plothFuncViaGam(gam.standard,"standard", "as_ln_c")
  g2 <- plothFuncViaGam(gps.gps, "gps", "as_ln_c")
  g3 <- plothFuncViaGam(gam.standard, "standard", "mn_ln_c")
  g4 <- plothFuncViaGam(gam.gps, "gps", "mn_ln_c")
  g5 <- plothFuncViaGam(gam.standard, "standard", "pb_ln_c")
  g6 <- plothFuncViaGam(gam.gps, "gps", "pb_ln_c")
  grid.arrange(g1, g3, g5, g2, g4, g6, ncol = 3)
}

if (isTRUE(is.runBKMR)) {
  # Begin analysis via bkrm().

  bgd$gps <- bgd$lm.trunc.weights # Use multivariate linear regression.

  # Set up response vector, and exposure and covariate matrices.
  response <- bgd$ccs_z
  covariates <- model.matrix(as.formula(paste("~ -1 +", cov.part)), data = bgd)
  gps.as.covariate <- bgd$gps

  # Set BKMR priors and control parameters.
  n.iterations <- n.bkmr.iterations
  control.params <- c(r.prior = "invunif", lambda.jump = 1,
                      list(r.params = c(list(r.a = 0, r.b = 1e2),
                                        list(r.jump = 0.07, r.jump1 = 2,
                                             r.jump2 = 0.1, r.muprop = 1))),
                      list(a.sigsq = 0.01, b.sigsq = 0.01))
  starting.value <- list(r = 0.02)

  # Fit different specifications of Bayesian kernel machine regression models.
  # Model 1: Adjusting confounders via the standard way.
  set.seed(500)
  y <- response
  Z <- exposures
  X <- covariates
  fitted.km.standard <- bkmr::kmbayes(y = y, Z = Z, X = X,
                                      iter = n.iterations,
                                      control.params = control.params,
                                      varsel = FALSE)
  save(fitted.km.standard, file = "fitted.km.standard.RData")

  # Model 2: Adjusting confounders via generalized propensity score regression.
  set.seed(500)
  y <- response
  Z <- exposures
  X <- gps.as.covariate
  fitted.km.gps <- bkmr::kmbayes(y = y, Z = Z, X = X,
                                 iter = n.iterations,
                                 control.params = control.params,
                                 varsel = FALSE)
  save(fitted.km.gps, file = "fitted.km.gps.RData")

  # Choose a model for result summary.
  bkmr.fit <- get(load("fitted.km.standard.RData"))

  filterData <- function(data.frame) {
    return(subset(data.frame, variable == "as_ln_c" | variable == "mn_ln_c" |
                    variable =="pb_ln_c"))
  }

  # Summarize MCMC results.
  samples <- bkmr::ExtractSamps(bkmr.fit, sel = NULL)
  par(mfrow = c(2, 2))
  bkmr::TracePlot(bkmr.fit, par = "r")
  bkmr::TracePlot(bkmr.fit, par = "lambda")
  bkmr::TracePlot(bkmr.fit, par = "beta")
  bkmr::TracePlot(bkmr.fit, par = "h")

  # Summarize posterior inclusion probabilities.
  pips <- bkmr::CalcPIPs(bkmr.fit)

  # Compare estimated h function when all predictors are at a particular quantile
  # to when all are at the 50th percentile.
  overall.risk <- bkmr::OverallRiskSummaries(fit = bkmr.fit, y = y, Z = Z, X = X)
  ggplot(overall.risk,
         aes(quantile, est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd)) +
    geom_pointrange() +
    ggtheme.config("", 25, 0)

  # Compute summaries of the risks associated with a change in a single variable
  # in Z from a 75th to a 25th percentile, for the other variables in Z fixed to a
  # specific level.
  single.risk <- bkmr::SingVarRiskSummaries(fit = bkmr.fit, y = y, Z = Z, X = X)
  ggplot(filterData(single.risk),
         aes(variable, est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd,
             col = q.fixed)) +
    geom_pointrange(position = position_dodge(width = 0.75)) +
    scale_y_continuous(limits = c(-0.6, 0.4)) +
    coord_flip() +
    ggtheme.config("", 25, 0)

  # Compare the single-predictor health risks when all of the other predictors in
  # Z are fixed to their 75th percentile to when all of the other predictors in Z
  # are fixed to their 25th percentile.
  interact.risk <- bkmr::SingVarIntSummaries(fit = bkmr.fit, y = y, Z = Z, X = X)
  ggplot(filterData(interact.risk),
         aes(variable, est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd)) +
    geom_hline(yintercept = 0, col = "brown", lty = 2) +
    geom_pointrange() +
    ggtheme.config("", 25, 0)

  # Plot univariate predictor-response function on a new grid of points.
  univariate.h <- PredictorResponseUnivar(fit = bkmr.fit, y = y, Z = Z, X = X)
  ggplot(filterData(univariate.h),
         aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
    geom_smooth(stat = "identity") +
    facet_wrap(~ variable) +
    ggtheme.config("", 25, 0)
}
