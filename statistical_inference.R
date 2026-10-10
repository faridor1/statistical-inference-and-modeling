# Statistical Inference in R
# Run from the repository root (the folder containing this file).
# Install dependencies once: install.packages(c("MASS", "cmprsk", "survival", "dplyr", "ggplot2", "emmeans"))

install.packages(c("MASS", "cmprsk", "survival", "dplyr", "ggplot2", "emmeans"))
library(MASS)
library(cmprsk)
library(survival)
library(dplyr)
library(ggplot2)
library(emmeans)


# ====================================================================
# 01 — Blood pressure: linear and logistic regression
# ====================================================================
setwd("...")
data1 <- read.table("data/data 1.txt", header = TRUE, sep = " ")

vars <- c("Blood_Pressure", "Age", "BMI", "Exercise_Frequency")

# Summarize each variable and construct a 95% t confidence interval for its
# population mean. The standard error is estimated using the sample variance.
summary_data1 <- do.call(rbind, lapply(vars, function(v) {
  x    <- data1[[v]]
  n    <- length(x)
  m    <- mean(x)
  vari <- var(x)
  ci   <- t.test(x, conf.level = 0.95)$conf.int
  data.frame(
    n = n,
    Variable = v,
    Mean     = round(m, 3),
    Variance = round(vari, 3),
    CI95_lower = round(ci[1], 3),
    CI95_upper = round(ci[2], 3)
  )
}))

print(summary_data1, row.names = FALSE)

# OLS estimates the conditional mean of blood pressure given the predictors.
# summary(model) reports coefficient standard errors, t-tests, and model fit.
model <- lm(Blood_Pressure ~ Age + BMI + Exercise_Frequency, data = data1)
summary(model)

# R-squared measures explained sample variation; adjusted R-squared
# additionally penalizes model complexity.
summary(model)$r.squared

summary(model)$adj.r.squared

# Define the binary outcome: blood pressure above 120 mmHg.
# Logistic regression models log-odds and estimates coefficients by MLE.
data1$BP_binary <- ifelse(data1$Blood_Pressure > 120, 1, 0)

# Each logistic coefficient is a conditional log-odds change; exp(beta)
# gives the corresponding odds ratio, holding other predictors fixed.
model_logit <- glm(BP_binary ~ Age + BMI + Exercise_Frequency, 
                   data = data1, family = binomial)
summary(model_logit)


# ====================================================================
# 02 — Beetle counts: Poisson and negative binomial regression
# ====================================================================

data2 <- read.table("data/data 2.txt", header = TRUE)

# Report observed means and variances (and invasive-species prevalence).
# These describe the sample; they are not conditional regression estimates.
summary_table <- data.frame(
  Beetle_Count = c(
    mean(data2$Beetle_Count, na.rm = TRUE),
    var(data2$Beetle_Count, na.rm = TRUE)
  ),
  Burn_Intensity = c(
    mean(data2$Burn_Intensity, na.rm = TRUE),
    var(data2$Burn_Intensity, na.rm = TRUE)
  ),
  Invasive_Species = c(
    mean(data2$Invasive_Species, na.rm = TRUE),  
    NA                                         
  )
)

rownames(summary_table) <- c("mean/prevalence", "variance")
summary_table

# Poisson GLM with log link: log(E[Y | X]) = X beta.
# Assumes conditional variance equals the conditional mean.
# Poisson regression
pois_model <- glm(Beetle_Count ~ Burn_Intensity + Invasive_Species,
                  family = poisson(link = "log"),
                  data = data2)

summary(pois_model)

# Pearson dispersion is approximately 1 under a well-fitting Poisson
# variance assumption; values substantially above 1 suggest overdispersion.
# Check overdispersion
dispersion <- sum(residuals(pois_model, type = "pearson")^2) /
              pois_model$df.residual
dispersion

# Negative binomial regression allows Var(Y | X) = mu + mu^2 / size,
# making it suitable for overdispersed counts.
# Negative Binomial regression
nb_model <- glm.nb(Beetle_Count ~ Burn_Intensity + Invasive_Species,
                   data = data2)

summary(nb_model)

# Custom MLE: estimate regression coefficients and the dispersion parameter
# by minimizing the negative log-likelihood with BFGS.
# Negative binomial MLE: log(size) parametrization enforces size > 0.
# Sum log probabilities rather than log(product(probabilities)) for stability.
log_likelihood <- function(par, Y, X1, X2) {
  mu <- exp(par[1] + par[2] * X1 + par[3] * X2)
  size <- exp(par[4])
  -sum(dnbinom(Y, mu = mu, size = size, log = TRUE))
}

initial_values <- c(0, 0, 0, log(0.5))
fit <- optim(
  par = initial_values, fn = log_likelihood, method = "BFGS",
  Y = data2$Beetle_Count, X1 = data2$Burn_Intensity,
  X2 = data2$Invasive_Species
)

# Coefficients and dispersion (size) on the original scale
beta_hat <- c(fit$par[1:3], size = exp(fit$par[4]))
beta_hat

# AIC = 2*k - 2*maximized_log_likelihood: lower values favor the
# better fit/complexity trade-off among these candidate models.
# Compare fit using AIC (for Poisson vs glm.nb)
AIC(pois_model, nb_model)

# Unconditional variance/mean ratio: descriptive only; unlike Pearson
# dispersion, this does not adjust for covariates.
# Check overdispersion in the outcome
var(data2$Beetle_Count) / mean(data2$Beetle_Count)


# ====================================================================
# 03 — Competing risks: cumulative incidence and hazard models
# ====================================================================

data3 <- read.table("data/data 3.txt", header = TRUE)
data3$Treatment <- factor(data3$Treatment)

# Status 0 = censored, 1 = recovery, 2 = death.
# The cumulative incidence function estimates P(T <= t, cause = j),
# accounting for the fact that competing events prevent the target event.
cif_fit <- cuminc(
  ftime   = data3$Time,
  fstatus = data3$Status,
  group   = data3$Treatment
)

# check order 
names(cif_fit)

# Rename CIF curve components
nm <- names(cif_fit)
nm[nm == "A 1"] <- "Treatment A - Recovery"
nm[nm == "B 1"] <- "Treatment B - Recovery"
nm[nm == "A 2"] <- "Treatment A - Death"
nm[nm == "B 2"] <- "Treatment B - Death"
names(cif_fit) <- nm

# Write the figure to a fixed-size image instead of the RStudio plot pane.
# This avoids "figure margins too large" when the pane is too small.
# The fitted CIFs and plotting arguments remain unchanged.
dir.create("figures", showWarnings = FALSE)
png("figures/competing_risks_cif.png", width = 1200, height = 800, res = 150)
plot(cif_fit,
     xlab = "Time (days)",
     ylab = "CIF",
     main = "Cumulative Incidence Functions for Competing Risks",
     lty  = c(1,2,3,4),
     lwd  = 2)
dev.off()

data3$Treatment <- relevel(data3$Treatment, ref = "A")

# Fine-Gray regression models the subdistribution hazard associated with
# the recovery cumulative incidence (event code 1).
# Fine-Gray model (subdistribution hazard)
fg_fit <- crr(
  ftime  = data3$Time,
  fstatus = data3$Status,
  cov1   = model.matrix(~ Treatment, data3)[, -1]
)

summary(fg_fit)

# Cause-specific Cox regression models the instantaneous recovery hazard
# among people still free of any event; deaths exit its risk set.
# Cause-specific Cox-model
cox_fit <- coxph(
  Surv(Time, Status == 1) ~ Treatment,
  data = data3
)

summary(cox_fit)


# ====================================================================
# 04 — Treatment comparisons: ANOVA and simultaneous contrasts
# ====================================================================

data4 <- read.table("data/data 4.txt", header = TRUE)

# Treatment as factor
data4$treatment_group <- factor(data4$treatment_group)

# Group means and individual 95% t confidence intervals; these are not
# adjusted for simultaneous coverage across multiple groups.
summary_stats <- data4 %>%
  group_by(treatment_group) %>%
  summarise(
    n = n(),
    mean = mean(blood_pressure_reduction),
    variance = var(blood_pressure_reduction),
    sd = sd(blood_pressure_reduction),
    se = sd / sqrt(n),
    lower_CI = mean + qt(0.025, df = n - 1) * se,
    upper_CI = mean + qt(0.975, df = n - 1) * se
  )

summary_stats

ggplot(summary_stats,
       aes(x = treatment_group, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
                width = 0.2) +
  labs(title = "Mean Blood Pressure reduction with 95% CI",
       x = "Treatment Group",
       y = "Mean reduction (mmHg)") +
  theme_minimal()

# One-way ANOVA tests H0: all five population treatment means are equal.
# A significant F-test alone does not identify which groups differ.
# Fit ANOVA
anova_model <- aov(blood_pressure_reduction ~ treatment_group, data = data4)
summary(anova_model)

# In this one-factor model, estimated marginal means equal group means.
# Compute estimated marginal means
emm <- emmeans(anova_model, "treatment_group")
summary(emm)

# Estimated means for each treatment group
emm <- emmeans(anova_model, ~ treatment_group)

# Each contrast estimates a treatment mean difference relative to A.
# Contrast coefficients sum to zero; e.g. (-1, 1, 0, 0, 0) = mu_B - mu_A.
# contrast matrix
custom_contrasts <- list(
  "B-A" = c(-1, 1, 0, 0, 0),
  "C-A" = c(-1, 0, 1, 0, 0),
  "D-A" = c(-1, 0, 0, 1, 0),
  "E-A" = c(-1, 0, 0, 0, 1)
)

# The mvt adjustment accounts for correlation among the four contrasts
# and targets simultaneous familywise inference rather than separate tests.
# Compute contrasts with multivariate-t simultaneous adjustment (max-t)
res <- contrast(emm, method = custom_contrasts, adjust = "mvt")
summary(res)

# Simultaneous confidence intervals
confint(res, adjust = "mvt")

