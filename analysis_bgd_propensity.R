# This script aims to re-analyse the Bangladesh neurodevelopment outcome via
# bkmr() and gam() by adjusting confounders using generalized propensity score
# regression. The generalized propensity scores are calculated using three
# methods including multivariate linear regression, lm() or heavyLm(), and
# multivariate normal additive model, mvn in gam().

# Load required libraries.
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(heavy))
suppressMessages(library(mgcv))
suppressMessages(library(mvtnorm))

n.outcomes <- 3 # Number of neurodevelopment outcomes.
kTruncateValue <- 0.05 # Value for truncating propensity scores.

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
bgd <- data %>%
  dplyr::select(ccs_z, as_ln_c, mn_ln_c, pb_ln_c, momage_c, homescore_z,
                momIQ_z, age_c, momeducd, smokenv, sex, protein, clinic) %>%
  na.omit

# Standardize exposure vectors.
bgd$as_ln_z <- with(bgd, (as_ln_c - mean(as_ln_c)) / sd(as_ln_c))
bgd$mn_ln_z <- with(bgd, (mn_ln_c - mean(mn_ln_c)) / sd(mn_ln_c))
bgd$pb_ln_z <- with(bgd, (pb_ln_c - mean(pb_ln_c)) / sd(pb_ln_c))

# Standardize continuous confounder vectors.
bgd$momage_z <- with(bgd, (momage_c - mean(momage_c)) / sd(momage_c))
bgd$age_z <- with(bgd, (age_c - mean(age_c)) / sd(age_c))
bgd$momeducd <- factor(bgd$momeducd)
bgd$smokenv <- factor(bgd$smokenv)
bgd$sex <- factor(bgd$sex)
bgd$protein <- factor(bgd$protein)
bgd$clinic <- factor(bgd$clinic)

# Finalize the analysis data.
bgd <- dplyr::select(bgd, ccs_z, as_ln_z, mn_ln_z, pb_ln_z, momage_z,
                     homescore_z, momIQ_z, age_z, momeducd, smokenv, sex,
                     protein, clinic)

# Set up exposure matrices and formulaes for models.
exposures <- cbind(bgd$as_ln_z, bgd$mn_ln_z, bgd$pb_ln_z)
num.lm.formula <- as.formula("cbind(as_ln_z, mn_ln_z, pb_ln_z) ~ 1")
cov.part <- paste(stringr::str_c("momage_z + homescore_z + momIQ_z + age_z +",
                                 "momeducd + smokenv + sex + protein + clinic"))
den.lm.formula <- as.formula(paste("cbind(as_ln_z, mn_ln_z, pb_ln_z) ~",
                                   cov.part))
num.gam.formula <- list(as_ln_z ~ 1, mn_ln_z ~ 1, pb_ln_z ~ 1)
den.gam.formula <- list(as.formula(paste("as_ln_z ~", cov.part)),
                        as.formula(paste("mn_ln_z ~", cov.part)),
                        as.formula(paste("pb_ln_z ~", cov.part)))

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
                          family = Student(df = 1e300))
    den.fitted <- heavyLm(den.lm.formula, data = bgd,
                          family = Student(df = 1e300))
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
par(mfrow = c(1, 2))
plot(bgd$lm.trunc.weights, bgd$gam.trunc.weights)
plot(bgd$lm.trunc.weights, bgd$hlm.trunc.weights)

# Begin analysis via gam().
Main <- function(weights) {
  outcome.part <- "ccs_z ~ s(as_ln_z, mn_ln_z, pb_ln_z) +"
  gam.standard <- gam(as.formula(paste(outcome.part, cov.part)), data = bgd)

  bgd$gps <- weights
  gam.regression <- gam(as.formula(paste(outcome.part, "s(gps, k = 3)")),
                        data = bgd)

  plothFuncViaGam <- function(gam.fit, method, exposure.var) {
    is_Mn <- ifelse(exposure.var == "mn_ln_z", "y", "n")
    mn <- switch(is_Mn,
                 y = seq(min(bgd$mn_ln_z), max(bgd$mn_ln_z), length = 100),
                 n = mean(gam.fit$model$mn_ln_z))

    is_As <- ifelse(exposure.var == "as_ln_z", "y", "n")
    as <- switch(is_As,
                 y = seq(min(bgd$as_ln_z), max(bgd$as_ln_z), length = 100),
                 n = mean(gam.fit$model$as_ln_z))

    is_Pb <- ifelse(exposure.var == "pb_ln_z", "y", "n")
    pb <- switch(is_Pb,
                 y = seq(min(bgd$pb_ln_z), max(bgd$pb_ln_z), length = 100),
                 n = mean(gam.fit$model$pb_ln_z))

    if (method == "standard") {
      testdata <- data.frame(as_ln_z = as,
                             mn_ln_z = mn,
                             pb_ln_z = pb,
                             momage_z = mean(gam.fit$model$momage_z),
                             homescore_z = mean(gam.fit$model$homescore_z),
                             momIQ_z = mean(gam.fit$model$momIQ_z),
                             age_z = mean(gam.fit$model$age_z),
                             momeducd = 1,
                             smokenv = 1,
                             sex = 1,
                             protein = 1,
                             clinic = 1)
    } else if (method == "regression") {
      testdata <- data.frame(as_ln_z = as,
                             mn_ln_z = mn,
                             pb_ln_z = pb,
                             gps = mean(gam.fit$model$gps))
    }

    fit <- predict(gam.fit, newdata = testdata, type = "response", se = TRUE)
    predicts <- data.frame(testdata, fit)
    g <- ggplot(predicts, aes_string(x = exposure.var, y = "fit")) +
      geom_smooth(aes(ymin =  predicts$fit - 1.96 * (predicts$se.fit),
                      ymax = predicts$fit + 1.96 * predicts$se.fit),
                  fill = "gray80", size = 1, stat = "identity") +
      scale_x_continuous(limits = c(-2, 4)) +
      scale_y_continuous(limits = c(-1, 1)) +
      xlab(exposure.var) + ylab("estimate") +
      ggtheme.config("", 25, 0)
    return(g)
  }

  g1 <- plothFuncViaGam(gam.standard,"standard", "as_ln_z")
  g2 <- plothFuncViaGam(gam.regression, "regression", "as_ln_z")
  g3 <- plothFuncViaGam(gam.standard, "standard", "mn_ln_z")
  g4 <- plothFuncViaGam(gam.regression, "regression", "mn_ln_z")
  g5 <- plothFuncViaGam(gam.standard, "standard", "pb_ln_z")
  g6 <- plothFuncViaGam(gam.regression, "regression", "pb_ln_z")
  grid.arrange(g1, g3, g5, g2, g4, g6, ncol = 3)
}
