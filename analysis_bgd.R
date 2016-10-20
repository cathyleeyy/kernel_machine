#!/usr/local/bin/Rscript

# Load required libraries.
suppressMessages(library(bkmr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(mgcv))

# Set FLAGS for analysis.
FLAGS <- data.frame(is_restrict = FALSE,
                    n_iterations = 100,
                    prior_option = 1)

# Set configuration for themes (non-data components of the plot).
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

# Set path and load data.
DATA.PATH <- stringr::str_c("/Users/cathylee/Documents/Harvard/Postdoc/",
                            "EarlySigns/earlyS.data/bgd_data.txt")
data <- read.table(DATA.PATH, header = TRUE)
str(data)

# Set higher-order terms for continuous variables.
data$age_sq <- data$age_c ^ 2
data$momIQ_sq <- data$momIQ_z ^ 2
data$approxage_sq < -data$approxage ^ 2
data$homescore_sq <- data$homescore_z ^ 2

# Restrict to the first clinic.
if (FLAGS$is_restrict) {
  data <- dplyr::filter(data, clinic == 0)
}
attach(data)

# Number of sample size: 825.
# Neurodevelopmental outcomes: ccs_z, lcs_z, mcs_z.
# Arsenic, Manganese and Lead exposure: as_ln_c, mn_ln_c, pb_ln_c.
# Covariates: momage_c, homescore_z, momIQ_z, age_c, approxage, momeducd,
#             smokenv, sex, protein, clinic.
df <- na.omit(data.frame(repro_id,
                         ccs_z, lcs_z, mcs_z,
                         as_ln_c, mn_ln_c, pb_ln_c,
                         momage_c, homescore_z, momIQ_z, age_c, approxage,
                         momeducd, smokenv, sex, protein, clinic))

# Declare factor variables.
df$momeducd <- as.factor(df$moneducd)
df$smokenv <- as.factor(df$smokenv)
df$sex <- as.factor(df$sex)
df$protein <- as.factor(df$protein)
df$clinic <- as.factor(df$clinic)

# Set up response vector, and exposure and covariate matrices.
response.names <- c("repro_id", "ccs_z")
exposure.names <- c("as_ln_c", "mn_ln_c", "pb_ln_c")
response <- df$ccs_z
exposures <- as.matrix(df[, exposure.names])

# Seperate continuous and categorical covariates.
continuous.names <- c("momage_c", "homescore_z", "momIQ_z", "age_c")
categorical.names <- c("momeducd", "smokenv", "sex", "protein", "clinic")
continuous <- as.matrix(df[, continuous.names])
categorical <- as.matrix(df[, categorical.names])
covariates <- cbind(continuous, categorical)
cate.matrix <- model.matrix(~ momeducd + smokenv + sex + protein + clinic,
                            data = data.frame(categorical))

# Set options for BKMR priors and control parameters.
n.iterations <- FLAGS$n_iterations # Default: 50000
r.prior.opts <- c(rep("invunif", 3), "gamma")
r.params.opts <- list(opt1 = list(r.a = 0, r.b = 1e2),
                      opt2 = list(r.a = 0, r.b = 1e3),
                      opt3 = list(r.a = 0, r.b = 1e4),
                      opt4 = list(mu.r = 0.01, sigma.r = 0.01))

# Set BKMR priors and control parameters. This is to check sensitivity of
# results according to various parameters.
i <- FLAGS$prior_option
control.params <- c(r.prior = r.prior.opts[i], lambda.jump = 1,
                    list(r.params = c(r.params.opts[[i]],
                                      list(r.jump = 0.07, r.jump1 = 2,
                                           r.jump2 = 0.1, r.muprop = 1))),
                    list(a.sigsq = 0.01, b.sigsq = 0.01))
starting.value <- list(r = 0.02)

# Fit different specifications of Bayesian kernel machine regression models.
# Model 1: Put all exposures and covariates into the h function, without any
# variable selection.
set.seed(500)
y <- response
Z <- cbind(exposures, continuous, cate.matrix)
km_hAll_no_varsel <- bkmr::kmbayes(y = y, Z = Z,
                                   iter = n.iterations,
                                   control.params = control.params,
                                   varsel = FALSE)
save(km_hAll_no_varsel, file = "km_hAll_no_varsel.RData")

# Model 2: Put all exposures and continuous covariates into the h function,
# without any variable selection.
set.seed(500)
y <- response
Z <- cbind(exposures, continuous)
X <- cate.matrix
km_hCont_no_varsel <- bkmr::kmbayes(y = y, Z = Z, X = X,
                                    iter = n.iterations,
                                    control.params = control.params,
                                    varsel = FALSE)
save(km_hCont_no_varsel, file = "km_hCont_no_varsel.RData")

################################################################################
# Choose a model for result summary.
X <- matrix(0, length(response), 1) # Model 1.
X <- cate.matrix                    # Model 2.
bkmr.fit <- get(load("km_hAll_no_varsel.RData"))

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

# Plot bivariate predictor-response function on a new grid of points.
bivariate.h <- PredictorResponseBivar(fit = bkmr.fit, y = y, Z = Z, X = X)
bivariate.h.levels <-
  PredictorResponseBivarLevels(pred.resp.df = bivariate.h, Z = Z)

bivariate.df <- subset(bivariate.h.levels,
                       (variable1 == "as_ln_c" & variable2 == "mn_ln_c") |
                       (variable1 == "as_ln_c" & variable2 == "pb_ln_c") |
                       (variable1 == "mn_ln_c" & variable2 == "as_ln_c") |
                       (variable1 == "mn_ln_c" & variable2 == "pb_ln_c") |
                       (variable1 == "as_ln_c" & variable2 == "age_c") |
                       (variable1 == "mn_ln_c" & variable2 == "age_c") |
                       (variable1 == "pb_ln_c" & variable2 == "age_c") |
                       (variable1 == "as_ln_c" & variable2 == "momage_c") |
                       (variable1 == "mn_ln_c" & variable2 == "momage_c") |
                       (variable1 == "pb_ln_c" & variable2 == "momage_c") |
                       (variable1 == "as_ln_c" & variable2 == "momIQ_z") |
                       (variable1 == "mn_ln_c" & variable2 == "momIQ_z") |
                       (variable1 == "pb_ln_c" & variable2 == "momIQ_z"))

bivariate.h.df %>% ggplot(aes(z1, est)) +
  geom_smooth(aes(col = quantile), stat = "identity") +
  facet_grid(variable1 ~ variable2) +
  ggtitle("f(z1 | quantiles of z2)") +
  ggtheme.config("", 25, 0)

################################################################################
# Set up data for GAM analysis.
# Model 1: Put only exposures into the s function.
gam_hZ_no_varsel <- gam(ccs_z ~ s(mn_ln_c, as_ln_c, pb_ln_c) + momage_c +
                          homescore_z + momIQ_z + age_c + momeducd + smokenv +
                          sex + protein + clinic, data = df)

# Model 2: Put all exposures and continuous covariates into the s function.
gam_hCont_no_varsel <- gam(ccs_z ~ s(mn_ln_c, as_ln_c, pb_ln_c, momage_c,
                                     homescore_z, momIQ_z, age_c) +
                             momeducd + smokenv + sex + protein + clinic,
                           data = df)

plothFuncViaGam <- function(gam.fit, exposure.var) {
  is_Mn <- ifelse(exposure.var == "mn_ln_c", "y", "n")
  mn <- switch(is_Mn,
               y = seq(min(df$mn_ln_c), max(df$mn_ln_c), length = 100),
               n = mean(gam.fit$model$mn_ln_c))

  is_As <- ifelse(exposure.var == "as_ln_c", "y", "n")
  as <- switch(is_As,
               y = seq(min(df$as_ln_c), max(df$as_ln_c), length = 100),
               n = mean(gam.fit$model$as_ln_c))

  is_Pb <- ifelse(exposure.var == "pb_ln_c", "y", "n")
  pb <- switch(is_Pb,
               y = seq(min(df$pb_ln_c), max(df$pb_ln_c), length = 100),
               n = mean(gam.fit$model$pb_ln_c))

  testdata <- data.frame(mn_ln_c = mn,
                         as_ln_c = as,
                         pb_ln_c = pb,
                         momage_c = mean(gam.fit$model$momage_c),
                         homescore_z = mean(gam.fit$model$homescore_z),
                         momIQ_z = mean(gam.fit$model$momIQ_z),
                         age_c = mean(gam.fit$model$age_c),
                         momeducd = 1,
                         smokenv = 1,
                         sex = 1,
                         protein = 1,
                         clinic = 1)

  fit <- predict(gam.fit, newdata = testdata, type = "response", se = TRUE)
  predicts <- data.frame(testdata, fit)
  g <- ggplot(predicts, aes_string(x = exposure.var, y = "fit")) +
    geom_smooth(aes(ymin =  predicts$fit - 1.96 * (predicts$se.fit),
                    ymax = predicts$fit + 1.96 * predicts$se.fit),
                fill = "gray80", size = 1, stat = "identity") +
    scale_x_continuous(limits = c(-2, 4)) +
    scale_y_continuous(limits = c(-1, 1)) +
    xlab(exposure.var) + ylab("h(x)") +
    ggtheme.config("", 25, 0)
  return(g)
}

g1 <- plothFuncViaGam(gam_hZ_no_varsel, "as_ln_c")
g2 <- plothFuncViaGam(gam_hCont_no_varsel, "as_ln_c")
g3 <- plothFuncViaGam(gam_hZ_no_varsel, "mn_ln_c")
g4 <- plothFuncViaGam(gam_hCont_no_varsel, "mn_ln_c")
g5 <- plothFuncViaGam(gam_hZ_no_varsel, "pb_ln_c")
g6 <- plothFuncViaGam(gam_hCont_no_varsel, "pb_ln_c")
grid.arrange(g1, g3, g5, g2, g4, g6, ncol = 3)
grid.arrange(g1, g3, g5, ncol = 3)
grid.arrange(g2, g4, g6, ncol = 3)


# par(mfrow = c(1, 3))
# vis.gam(gam_hAll_no_varsel, view = c("mn_ln_c", "as_ln_c"),
#         plot.type = "contour", xlab = "Mn", ylab = "As", main = "h(x)")
# vis.gam(gam_hAll_no_varsel, view = c("pb_ln_c", "mn_ln_c"),
#         plot.type = "contour", xlab = "Pb", ylab = "Mn", main = "h(x)")
# vis.gam(gam_hAll_no_varsel, view = c("as_ln_c", "pb_ln_c"),
#         plot.type = "contour", xlab = "As", ylab = "Pb", main = "h(x)")

# Test for arsenic and manganese interaction, we compare two specifications of
# GAM employing the tensor product smoother, ti(): one that includes the
# smoothed terms for the metals additively and the other that allows for a joint
# smoothed effect of arsenic and manganese.
# gam.nointeraction.fit <- gam(ccs_z ~ ti(mn_ln_c, k = 4) + ti(as_ln_c, k = 4) +
#                               ti(pb_ln_c, k = 4) + gender + age_c + age_sq +
#                               approxage + approxage_sq + momIQ_z + momIQ_sq +
#                               homescore_z + homescore_sq + momeducd + smokenv +
#                               protein, data = data)
#
# gam.interaction.fit <- gam(ccs_z ~ ti(mn_ln_c, k = 4) + ti(as_ln_c, k = 4) +
#                              ti(as_ln_c, mn_ln_c, k = 4) + ti(pb_ln_c, k = 4) +
#                              gender + age_c + age_sq + approxage +
#                              approxage_sq + momIQ_z + momIQ_sq + homescore_z +
#                              homescore_sq + momeducd + smokenv + protein,
#                            data = data)
#
# anova(gam.nointeraction.fit, gam.interaction.fit, test = "Chisq")
#
# par(mfrow = c(2, 2))
# vis.gam(gam.interaction.fit, view = c("mn_ln_c", "as_ln_c"),
#         plot.type = "contour", xlab = "Mn", ylab = "As", main = "h(Mn,As)")
# vis.gam(gam.nointeraction.fit, view = c("mn_ln_c", "as_ln_c"),
#         plot.type = "contour", xlab = "Mn", ylab = "As", main = "h(Mn) + h(As)")
# plot(gam.nointeraction.fit, ylim = c(-2, 1))






























