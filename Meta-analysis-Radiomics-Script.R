# Load the necessary packages
library("mada")
library("meta")
library(mvmeta)
library("readxl")
library(Matrix)
library(lme4)
library(lmtest)
library(usethis)
library(devtools)
library(MASS)
library(dplyr)
library(ggplot2)
library(metafor)

# Meta-analysis with subgroups
# Replace 'data$TP', 'data$TN', 'data$FP', 'data$FN' with standard variable names
# e.g., 'data$TP' -> 'TP', 'data$FN' -> 'FN'

sensitivity_model <- metaprop(TP, TP + FN,
                              random = TRUE,
                              fixed = FALSE, 
                              sm = "PLOGIT",
                              method.ci = "CP",
                              studlab = study_ID,
                              subgroup = study_group)

specificity_model <- metaprop(TN, TN + FP,
                              random = TRUE,
                              fixed = FALSE, 
                              sm = "PLOGIT",
                              method.ci = "CP",
                              studlab = study_ID,
                              subgroup = study_group)

# AUC Calculation using dmetatools package
# Make sure the standard variable names are used here as well
AUC_bootstrap <- dmetatools::AUC_boot(TP, FP, FN, TN)

# SROC Calculation
# SROC fitting using Reitsma model
fit_data <- reitsma(data)
summary(fit_data)

# SROC plot
plot(fit_data, xlim = c(0, 1), ylim = c(0.5, 1), main = "SROC Curve")
lines(sroc(fit_data), lty = 1)
points(fpr(data), sens(data), pch = 1, cex = .5)
legend("bottomright", c("Individual studies"), pch = 1)
legend("bottomleft", c("95% CI region"), lty = c(1))

# Deeks' DOR (Diagnostic Odds Ratio)
data$DOR <- (TP / FN) / (FP / TN)
data$logDOR <- log(data$DOR)
data$invSqrtVariance <- 1 / sqrt(1 / TP + 1 / FP + 1 / FN + 1 / TN)
data <- data[!is.na(data$logDOR) & !is.infinite(data$logDOR) & 
               !is.na(data$invSqrtVariance) & !is.infinite(data$invSqrtVariance), ]

# Leave-one-out Analysis
data_inf <- metainf(data, pooled = "random")
forest(data_inf,
       xlim = c(0.5, 1),
       leftcols = c("studlab"),
       rightcols = c("effect", "ci", "tau2", "tau", "I2"),
       rightlabs = c("Sensitivity", "95%CI"))

# Threshold effect calculation
data$threshold_sensitivity <- TP / (TP + FN)
data$one_minus_spec <- 1 - (TN / (TN + FP))

# Spearman correlation test for threshold effect
spearman_test <- cor.test(data$threshold_sensitivity, data$one_minus_spec, method = "spearman")

# Meta-regression with covariate
rma_sens <- rma.uni(yi = logit(Sensitivity), sei = SE_Sensitivity, mods = ~ Covariate, method = "REML", data = data)
rma_spec <- rma.uni(yi = logit(Specificity), sei = SE_Specificity, mods = ~ Covariate, method = "REML", data = data)

# Complementary Analysis using Hartung-Knapp-Sidik-Jonkman (HKSJ) method
# Logit transformation (PLO = proportion logit)
escalado_data <- escalc(measure = "PLO", xi = TP, ni = TP + FN, data = MRI)

# Random effects model using HKSJ
modelo_HKSJ <- rma(yi, vi, data = escalado_data, method = "DL", test = "knha")

# Model summary
summary(modelo_HKSJ)

# Backtransform to get the pooled sensitivity in original scale
predict(modelo_HKSJ, transf = transf.ilogit)
