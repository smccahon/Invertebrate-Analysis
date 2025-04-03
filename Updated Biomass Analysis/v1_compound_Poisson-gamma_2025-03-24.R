#------------------------------#
#    Invert Biomass Analysis   #
#     Created 2025-03-24       #
#    Modified 2025-03-27       #
#------------------------------#

# Load packages & data ----
library(tidyverse)
library(trtools)
library(MASS)
library(MuMIn)
library(cplm)
library(ggbreak) 
library(patchwork)
library(car)
library(data.table)
library(viridis)

invert <- read.csv("../data/Macroinverterbrate_Analysis_2025-03-06.csv")

# ...transform variables ----

invert$Season <- factor(invert$Season, levels = c("Spring", "Fall"))

invert$Permanence <- factor(invert$Permanence, 
                            levels = c("Temporary", "Seasonal",
                                       "Semipermanent", "Permanent"))

invert$AgCategory <- factor(invert$AgCategory,
                            levels = c("Low", "Moderate", "High"))

invert$DominantCrop <- factor(invert$DominantCrop,
                              levels = c("Grassland",
                                         "Canola",
                                         "Corn",
                                         "Soybean",
                                         "Wheat"))

invert <- invert %>% 
  mutate(PercentAg = PercentAg * 100)

# ...standardize continuous variables ----
invert.cs <- invert
invert.cs$PercentAg <- scale(invert.cs$PercentAg)
invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
invert.cs$pH <- scale(invert.cs$pH)
invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)

# ...correlations and covariate selection ----
# code from biomass_analysis_2025-03-06

# Season and permanence (-0.63)
# Distance to nearest crop and percent agriculture (-0.72)
# agricultural category and percent agriculture (0.89)
# Dominant crop and percent agriculture (0.63)
# agricultural category and nearest crop distance (-0.66)
# permanence and percent buffer (0.668)

### KATIE: If two variables are highly correlated, use AIC to tell you which one is better
### Can indicate all the variables you considered in methods, and say you dropped the one with less explanatory power
### Use season as informed null
### Test the effect of permanence on biomass AFTER separately (another analysis given the correlation with season)

# Agriculture
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ AgCategory + Season, link = "log", data = invert.cs)
m3 <- cpglm(Biomass ~ NearestCropDistance_m + Season, link = "log", 
            data = invert.cs)
m4 <- cpglm(Biomass ~ DominantCrop + Season, link = "log", 
            data = invert.cs)

### ...AIC 
models <- list(m1, m2, m3, m4)
model.sel(models)

# PercentAg is the best predictor with and without outliers

# Compound Poisson-gamma model ----
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", 
            data = invert.cs)
summary(m1)
plot(resid(m1))
plot(predict(m1),rstudent(m1)) # does not work

# ...R documentation code ----
# residual and qq plot
parold <- par(mfrow = c(2, 2), mar = c(5, 5, 2, 1))

# 1. regular plot
r1 <- resid(m1) / sqrt(m1$phi)
plot(r1 ~ fitted(m1), cex = 0.5)
qqnorm(r1, cex = 0.5)

# 2. quantile residual plot to avoid overlapping
u <- tweedie::ptweedie(m1$y, m1$p, fitted(m1), m1$phi)
u[m1$y == 0] <- runif(sum(m1$y == 0), 0, u[m1$y == 0])
r2 <- qnorm(u)
plot(r2 ~ fitted(m1), cex = 0.5)
qqnorm(r2, cex = 0.5)
par(parold)



# ...model comparison ----
invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m3 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m4 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season + 
              MaxBufferWidth_m, 
            link = "log", 
            data = invert.cs)
m5 <- cpglm(Biomass ~ pH + Season, link = "log", data = invert.cs)
m6 <- cpglm(Biomass ~ TDS_mg.L + Season, link = "log", data = invert.cs)
m7 <- cpglm(Biomass ~ pH + TDS_mg.L + Season, link = "log", data = invert.cs)
m8 <- cpglm(Biomass ~ Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m9 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)

# ...AIC ----
models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
model_comparison <- model.sel(models)
print(model_comparison)

# ...deviance ----
deviance(m1)
deviance(m2)
deviance(m3)
deviance(m4)
deviance(m5)
deviance(m6)
deviance(m7)
deviance(m8)
deviance(m9)

deviance_explained <- 1 - (deviance(m1) / deviance(m9))


# ...plot top models ----
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert)

# create data frame
d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                 max(invert$PercentAg), 
                                 length.out = 1000),
                 Season = unique(invert$Season))

ci <- trtools::lincon(m1, fcov=vcov)

### ...add confidence intervals ----

# Generate predictions on the log scale using type = "link"
d$yhat_log <- predict(m1, newdata = d, type = "link")

# Extract the covariance matrix from the model
vcov_matrix <- vcov(m1)

# Create the design matrix for the new data
X <- model.matrix(~ PercentAg + Season, data = d)

# Calculate the variance for each prediction (on the log scale)
pred_var <- diag(X %*% vcov_matrix %*% t(X))

# Calculate standard errors on the log scale
pred_se_log <- sqrt(pred_var)

# Calculate the 95% confidence intervals on the log scale
d$lower_CI_log <- d$yhat_log - 1.96 * pred_se_log
d$upper_CI_log <- d$yhat_log + 1.96 * pred_se_log

# Exponentiate the predictions and confidence intervals to get them 
# on the original scale
d$yhat <- exp(d$yhat_log)
d$lower_CI <- exp(d$lower_CI_log)
d$upper_CI <- exp(d$upper_CI_log)

# Plot predictions with confidence intervals
# Biomass ~ Ag + Season Plot ----
ggplot(d, aes(x = PercentAg, y = yhat, col = Season)) +  
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = Season), 
              alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) +
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = PercentAg, y = Biomass),
             size = 1.5) +
  scale_y_break(c(1.00,1.40), scales =0.25, 
                ticklabels = c(1.50, 3.00),
                space = 0.5)

# ...double check confidence intervals ----
# Check the manually exponentiated predictions and compare them to the 
# ones in the data
expected_yhat <- exp(d$yhat_log)
all.equal(d$yhat, expected_yhat)  # Should return TRUE if everything is correct

# Check if the manually calculated confidence intervals match
expected_lower_CI <- exp(d$lower_CI_log)
expected_upper_CI <- exp(d$upper_CI_log)

all.equal(d$lower_CI, expected_lower_CI)  # Should return TRUE
all.equal(d$upper_CI, expected_upper_CI)  # Should return TRUE

# Manually calculate the lower and upper confidence intervals on the log scale

# Formula: CI_lower = yhat_log - 1.96 * SE_log
#          CI_upper = yhat_log + 1.96 * SE_log

expected_lower_CI_log <- d$yhat_log - 1.96 * pred_se_log
expected_upper_CI_log <- d$yhat_log + 1.96 * pred_se_log

# Check if these manually calculated values match the ones in d
all.equal(d$lower_CI_log, expected_lower_CI_log)  # Should return TRUE
all.equal(d$upper_CI_log, expected_upper_CI_log)  # Should return TRUE

# ...assumptions ----
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)


# Plot the deviance residuals to check for any patterns
plot(residuals(m1), main = "Deviance Residuals", ylab = "Residuals", 
     xlab = "Index") # no clear pattern


# Chi-square test for goodness of fit (m1) ----
# Calculate the expected values for the biomass response
expected_values <- predict(m1, type = "response")

# Calculate the Pearson residuals
pearson_residuals <- (invert.cs$Biomass - expected_values) / sqrt(expected_values)

# Calculate the Pearson chi-square statistic
pearson_chi_square <- sum(pearson_residuals^2)

df <- 77-3
p_value <- 1 - pchisq(27.29, df)
p_value # 0.99

# Plot residuals vs fitted values ----
# generally uniform around 0
plot(fitted(m1), residuals(m1), main = "Residuals vs Fitted Values", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Q-Q plot of residuals
# mostly normally distributed
qqnorm(residuals(m1))
qqline(residuals(m1), col = "red")

# Overdispersion check: residual deviance / degrees of freedom
overdispersion_stat <- deviance(m1) / df.residual(m1) # 1.29

# Variance prop. to the square of the expected response ----

# Get the expected values (means) from the model
expected_values <- predict(m1, type = "response")

# Calculate the residuals (observed - predicted)
residuals <- invert.cs$Biomass - expected_values

# Calculate the observed variance (squared residuals)
observed_variance <- residuals^2

# Now compare the observed variance to the square of the expected values
comparison <- data.frame(expected = expected_values, observed_variance = observed_variance)

# Plot the expected vs observed variance
plot(comparison$expected^2, comparison$observed_variance, 
     main = "Observed Variance vs Square of Expected Value", 
     xlab = "Square of Expected Value", ylab = "Observed Variance")
abline(0, 1, col = "red")  # Add a 45-degree line for reference

# Calculate the variance-to-mean ratio
variance_to_mean_ratio <- observed_variance / expected_values

# Plot the variance-to-mean ratio vs the expected values
plot(expected_values, variance_to_mean_ratio, 
     main = "Variance-to-Mean Ratio vs Expected Value", 
     xlab = "Expected Value", ylab = "Variance-to-Mean Ratio")

# Plot parameter estimates ----
# Plot parameter estimates of top model
ci <- trtools::lincon(m1, fcov=vcov)
ci <- as.data.frame(ci)
ci$coefficient <- rownames(ci)

ggplot(data = ci, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = NULL) +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.9)) +
  scale_x_discrete(limits = c("PercentAg","SeasonFall","(Intercept)"),
                   labels = c("PercentAg" = "% Surrounding Ag", 
                              "SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

# ...influence of outliers ----
m <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert)



# Wetland permanence ----
# Compare to other models ####
m1 <- cpglm(Biomass ~ Permanence, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)
m3 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data= invert.cs)
m4 <- cpglm(Biomass ~ PercentAg + Permanence, link = "log", data = invert.cs)
m5 <- cpglm(Biomass ~ 1, link = "log", data = invert.cs)
  
# AIC MODEL SELECTION
models <- list(m1, m2, m3, m4)
model.sel(models)

deviance(m1)
deviance(m2)
deviance(m3)
deviance(m4)

summary(m1)
ci <- trtools::lincon(m1, fcov=vcov)
# Wetland permanence performs significantly better than informed null
# Significant difference in biomass between temporary and permanent wetlands
# Surrounding Ag + Season performs significantly better than wetland permanence

# plot parameter estimates and graph
ci <- trtools::lincon(m1, fcov=vcov)
ci <- as.data.frame(ci)
ci$coefficient <- rownames(ci)
ggplot(data = ci, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = "Wetland Permanence") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.9)) +
  scale_x_discrete(limits = c("PermanencePermanent","PermanenceSemipermanent",
                              "PermanenceSeasonal", "(Intercept)"),
                   labels = c("PermanencePermanent" = "Permanent",
                              "PermanenceSemipermanent" = "Semipermanent", 
                              "PermanenceSeasonal" = "Seasonal",
                              "(Intercept)" = "Temporary")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

# Plot permanence univariate model ----
m1 <- cpglm(Biomass ~ Permanence, link = "log", data = invert)

# create data frame
d <- expand.grid(Permanence = unique(invert$Permanence))

ci <- trtools::lincon(m1, fcov=vcov)

### ...add confidence intervals ----

# Generate predictions on the log scale using type = "link"
d$yhat_log <- predict(m1, newdata = d, type = "link")

# Extract the covariance matrix from the model
vcov_matrix <- vcov(m1)

# Create the design matrix for the new data
X <- model.matrix(~ Permanence, data = d)

# Calculate the variance for each prediction (on the log scale)
pred_var <- diag(X %*% vcov_matrix %*% t(X))

# Calculate standard errors on the log scale
pred_se_log <- sqrt(pred_var)

# Calculate the 95% confidence intervals on the log scale
d$lower_CI_log <- d$yhat_log - 1.96 * pred_se_log
d$upper_CI_log <- d$yhat_log + 1.96 * pred_se_log

# Exponentiate the predictions and confidence intervals to get them 
# on the original scale
d$yhat <- exp(d$yhat_log)
d$lower_CI <- exp(d$lower_CI_log)
d$upper_CI <- exp(d$upper_CI_log)

# Plot predictions with confidence intervals
# Biomass ~ Ag + Season Plot ----
ggplot(d, aes(x = Permanence, y = yhat)) +  
  geom_point(data = invert, aes(x = Permanence, y = Biomass, col = Permanence),
             size = 3)+
  geom_point(size = 4, col = "black") +  
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.1,
                col = "black",
                size = 1) +
  theme_classic() +
  labs(x = "Wetland Permanence", 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 12)),
        axis.title.y = element_text(size = 20, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_y_break(c(0.65,1.40), scales =0.25, 
                ticklabels = c(1.50, 3.00),
                space = 0.5) +
  scale_color_viridis_d(begin = 0.2, end = 0.95)
  
# remove outliers and rerun analysis ----
invert <- subset(invert, Biomass < 1.4)
invert.cs <- invert
invert.cs$PercentAg <- scale(invert.cs$PercentAg)
invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
invert.cs$pH <- scale(invert.cs$pH)
invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)

invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m3 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m4 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season + 
              MaxBufferWidth_m, 
            link = "log", 
            data = invert.cs)
m5 <- cpglm(Biomass ~ pH + Season, link = "log", data = invert.cs)
m6 <- cpglm(Biomass ~ TDS_mg.L + Season, link = "log", data = invert.cs)
m7 <- cpglm(Biomass ~ pH + TDS_mg.L + Season, link = "log", data = invert.cs)
m8 <- cpglm(Biomass ~ Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m9 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)

m10 <- cpglm(Biomass ~ PercentAg + MaxBufferWidth_m + PercentLocalVeg_50m +
              pH + TDS_mg.L + Dist_Closest_Wetland_m + Season, data = invert.cs,
             link = "log")

# ...AIC ----
models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
model_comparison <- model.sel(models)
print(model_comparison)

trtools::lincon(m7, fcov=vcov)



# ...deviance ----
deviance(m1)
deviance(m2)
deviance(m3)
deviance(m4)
deviance(m5)
deviance(m6)
deviance(m7)
deviance(m8)
deviance(m9)
deviance(m10)

# ...plot parameter estimates from 4 top models ----
ci1 <- trtools::lincon(m1, fcov=vcov)
ci1 <- as.data.frame(ci1)
ci1$coefficient <- rownames(ci1)

ci3 <- trtools::lincon(m3, fcov=vcov)
ci3 <- as.data.frame(ci3)
ci3$coefficient <- rownames(ci3)

ci4 <- trtools::lincon(m4, fcov=vcov)
ci4 <- as.data.frame(ci4)
ci4$coefficient <- rownames(ci4)

ci8 <- trtools::lincon(m8, fcov=vcov)
ci8 <- as.data.frame(ci8)
ci8$coefficient <- rownames(ci8)

p1 <- ggplot(data = ci1, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = NULL) +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 12)),
        axis.title.y = element_text(size = 18, margin = margin(r = 12)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(limits = c("PercentAg","SeasonFall","(Intercept)"),
                   labels = c("PercentAg" = "% Surrounding Ag", 
                              "SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))


p3 <- ggplot(data = ci3, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = NULL) +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 12)),
        axis.title.y = element_text(size = 18, margin = margin(r = 12)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(limits = c("PercentLocalVeg_50m","SeasonFall","(Intercept)"),
                   labels = c("PercentLocalVeg_50m" = "Local Vegetation Cover", 
                              "SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

p4 <- ggplot(data = ci4, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = NULL) +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 12)),
        axis.title.y = element_text(size = 18, margin = margin(r = 12)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(limits = c("PercentLocalVeg_50m","MaxBufferWidth_m","SeasonFall","(Intercept)"),
                   labels = c("PercentLocalVeg_50m" = "Local Vegetation Cover",
                              "MaxBufferWidth_m" = "Maximum Buffer Cover",
                              "SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

p8 <- ggplot(data = ci8, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = NULL) +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 12)),
        axis.title.y = element_text(size = 18, margin = margin(r = 12)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_x_discrete(limits = c("Dist_Closest_Wetland_m","SeasonFall","(Intercept)"),
                   labels = c("Dist_Closest_Wetland_m" = "Dist. to Closest Wetland",
                              "MaxBufferWidth_m" = "Maximum Buffer Cover",
                              "SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

p3 + p4 + p1 + p8

# ...is permanence still important without outliers? ----
m1 <- cpglm(Biomass ~ Permanence, link = "log", data = invert)
summary(m1)

# Compare to other models ####
m1 <- cpglm(Biomass ~ Permanence, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)
m3 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data= invert.cs)
m4 <- cpglm(Biomass ~ PercentAg + Permanence, link = "log", data = invert.cs)
m5 <- cpglm(Biomass ~ 1, link = "log", data = invert.cs)

# AIC MODEL SELECTION
models <- list(m1, m2)
model.sel(models)

deviance(m1)
deviance(m2)
deviance(m3)
deviance(m4)

summary(m1)
ci <- trtools::lincon(m1, fcov=vcov)
# Wetland permanence performs significantly better than informed null
# Significant difference in biomass between temporary and permanent wetlands
# Surrounding Ag + Season performs significantly better than wetland permanence

# plot parameter estimates and graph
ci <- trtools::lincon(m1, fcov=vcov)
ci <- as.data.frame(ci)
ci$coefficient <- rownames(ci)
ggplot(data = ci, aes(x = coefficient, y = estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", 
             lwd = 1.5) +  # Horizontal line at y=0
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "black", width = 0, 
                lwd = 1.5) +  # Confidence intervals as error bars
  coord_flip() +  # Flip axes to make the plot horizontal
  geom_point(size = 3) +  # Plot the point estimates
  theme_classic() +  # Clean theme
  labs(y = "Parameter Estimate", x = "Wetland Permanence") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.9)) +
  scale_x_discrete(limits = c("PermanencePermanent","PermanenceSemipermanent",
                              "PermanenceSeasonal", "(Intercept)"),
                   labels = c("PermanencePermanent" = "Permanent",
                              "PermanenceSemipermanent" = "Semipermanent", 
                              "PermanenceSeasonal" = "Seasonal",
                              "(Intercept)" = "Temporary")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

# ...outliers change inference...what should I do? ----


# models with informative parameters: stage 2 ----
invert <- subset(invert, Biomass < 1.4)
invert.cs <- invert
invert.cs$PercentAg <- scale(invert.cs$PercentAg)
invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
invert.cs$pH <- scale(invert.cs$pH)
invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)

invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")

# model selection
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m3 <- cpglm(Biomass ~ Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m4 <- cpglm(Biomass ~ MaxBufferWidth_m + Season, link = "log", 
                  data = invert.cs)
m5 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)

# two-way additive combinations
m6 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m7 <- cpglm(Biomass ~ PercentAg + Dist_Closest_Wetland_m + Season, link = "log", 
          data = invert.cs)
m8 <- cpglm(Biomass ~ PercentAg + MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m9 <- cpglm(Biomass ~ PercentLocalVeg_50m + Dist_Closest_Wetland_m + Season, 
            link = "log", data = invert.cs)
m10 <- cpglm(Biomass ~ PercentLocalVeg_50m + MaxBufferWidth_m + Season, 
            link = "log", data = invert.cs)
m11 <- cpglm(Biomass ~ Dist_Closest_Wetland_m + MaxBufferWidth_m + Season, 
            link = "log", data = invert.cs)

# three-way additive combinations
m12 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m + 
               Dist_Closest_Wetland_m + Season, link = "log", data = invert.cs)
m13 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m + 
               MaxBufferWidth_m + Season, link = "log", data = invert.cs)
m14 <- cpglm(Biomass ~ PercentAg + Dist_Closest_Wetland_m + 
               MaxBufferWidth_m + Season, link = "log", data = invert.cs)
m15 <- cpglm(Biomass ~ PercentLocalVeg_50m + Dist_Closest_Wetland_m + 
               MaxBufferWidth_m + Season, link = "log", data = invert.cs)

# four-way additive combinations
m16 <- cpglm(Biomass ~ PercentLocalVeg_50m + Dist_Closest_Wetland_m + 
               MaxBufferWidth_m + PercentAg + Season, link = "log", 
             data = invert.cs)


# does model performance improve with interactions? 
# ag + local veg
m1 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m2 <- cpglm(Biomass ~ PercentAg * PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)

models <- list(m1, m2)
model.sel(models) # model without interaction is better with and without outliers

# ag + buffer
m1 <- cpglm(Biomass ~ PercentAg + MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m2 <- cpglm(Biomass ~ PercentAg * MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)

models <- list(m1, m2)
model.sel(models) # model without interaction is better with and without outliers

# local veg + buffer
m1 <- cpglm(Biomass ~ PercentLocalVeg_50m + MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m2 <- cpglm(Biomass ~ PercentLocalVeg_50m * MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)

models <- list(m1, m2)
model.sel(models) # model with interaction is better with outliers
                  # model without interaction is better without outliers

# dist to wetland + ag
m1 <- cpglm(Biomass ~ PercentAg + Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m2 <- cpglm(Biomass ~ PercentAg * Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)

models <- list(m1, m2)
model.sel(models) # model without interaction is better with and without outliers

# dist to wetland + local veg
m1 <- cpglm(Biomass ~ PercentLocalVeg_50m + Dist_Closest_Wetland_m + Season, 
            link = "log", 
            data = invert.cs)
m2 <- cpglm(Biomass ~ PercentLocalVeg_50m * Dist_Closest_Wetland_m + Season, 
            link = "log", 
            data = invert.cs)

models <- list(m1, m2)
model.sel(models)  # model without interaction is better with and without outliers

# dist to wetland + buffer
m1 <- cpglm(Biomass ~ MaxBufferWidth_m + Dist_Closest_Wetland_m + Season, 
            link = "log", 
            data = invert.cs)
m2 <- cpglm(Biomass ~ MaxBufferWidth_m * Dist_Closest_Wetland_m + Season, 
            link = "log", 
            data = invert.cs)

models <- list(m1, m2)
model.sel(models) # model with interaction is better with outliers
                  # model without interaction is better without outliers



# compare informative parameters with relevant interactions

m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m3 <- cpglm(Biomass ~ Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m4 <- cpglm(Biomass ~ MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m5 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)

# two-way additive combinations
m6 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m7 <- cpglm(Biomass ~ PercentAg + Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m8 <- cpglm(Biomass ~ PercentAg + MaxBufferWidth_m + Season, link = "log", 
            data = invert.cs)
m9 <- cpglm(Biomass ~ PercentLocalVeg_50m + Dist_Closest_Wetland_m + Season, 
            link = "log", data = invert.cs)
m10 <- cpglm(Biomass ~ PercentLocalVeg_50m * MaxBufferWidth_m + Season, 
             link = "log", data = invert.cs)
m11 <- cpglm(Biomass ~ Dist_Closest_Wetland_m * MaxBufferWidth_m + Season, 
             link = "log", data = invert.cs)

# three-way additive combinations
m12 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m + 
               Dist_Closest_Wetland_m + Season, link = "log", data = invert.cs)
m13 <- cpglm(Biomass ~ PercentAg + PercentLocalVeg_50m * MaxBufferWidth_m +
               Season, link = "log", data = invert.cs)
m14 <- cpglm(Biomass ~ PercentAg + Dist_Closest_Wetland_m * MaxBufferWidth_m + 
               Season, link = "log", data = invert.cs)

# four-way additive combinations
m15 <- cpglm(Biomass ~ PercentLocalVeg_50m * MaxBufferWidth_m + 
               Dist_Closest_Wetland_m * MaxBufferWidth_m + PercentAg + 
               Season, link = "log", 
             data = invert.cs)

# ...AIC ----
models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, 
               m14, m15)

model_comparison <- model.sel(models)
print(model_comparison)

trtools::lincon(m8, fcov=vcov)
trtools::lincon(m14, fcov=vcov)
trtools::lincon(m13, fcov=vcov)
trtools::lincon(m1, fcov=vcov)
trtools::lincon(m15, fcov=vcov)
trtools::lincon(m7, fcov=vcov)
trtools::lincon(m6, fcov=vcov)
trtools::lincon(m12, fcov=vcov)
trtools::lincon(m11, fcov=vcov)
trtools::lincon(m10, fcov=vcov)
trtools::lincon(m5, fcov=vcov)
trtools::lincon(m4, fcov=vcov)
trtools::lincon(m2, fcov=vcov)
trtools::lincon(m3, fcov=vcov)
trtools::lincon(m9, fcov=vcov)

# Biomass ~ Local Veg + Season -------------------------------------------------
# outliers removed
 invert <- subset(invert, Biomass < 1.4)
# invert.cs <- invert
# invert.cs$PercentAg <- scale(invert.cs$PercentAg)
# invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
# invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
# invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
# invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
# invert.cs$pH <- scale(invert.cs$pH)
# invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
# invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)
# 
# invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")
m1 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season, link = "log", data = invert)

# create data frame
d <- expand.grid(PercentLocalVeg_50m = seq(min(invert$PercentLocalVeg_50m), 
                                 max(invert$PercentLocalVeg_50m), 
                                 length.out = 1000),
                 Season = unique(invert$Season))

ci <- trtools::lincon(m1, fcov=vcov)

### ...add confidence intervals ----

# Generate predictions on the log scale using type = "link"
d$yhat_log <- predict(m1, newdata = d, type = "link")

# Extract the covariance matrix from the model
vcov_matrix <- vcov(m1)

# Create the design matrix for the new data
X <- model.matrix(~ PercentLocalVeg_50m + Season, data = d)

# Calculate the variance for each prediction (on the log scale)
pred_var <- diag(X %*% vcov_matrix %*% t(X))

# Calculate standard errors on the log scale
pred_se_log <- sqrt(pred_var)

# Calculate the 95% confidence intervals on the log scale
d$lower_CI_log <- d$yhat_log - 1.96 * pred_se_log
d$upper_CI_log <- d$yhat_log + 1.96 * pred_se_log

# Exponentiate the predictions and confidence intervals to get them 
# on the original scale
d$yhat <- exp(d$yhat_log)
d$lower_CI <- exp(d$lower_CI_log)
d$upper_CI <- exp(d$upper_CI_log)

# Plot predictions with confidence intervals
# Biomass ~ Ag + Season Plot ----
ggplot(d, aes(x = PercentLocalVeg_50m, y = yhat, col = Season)) +  
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = Season), 
              alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Local Vegetation Cover", 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) +
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = PercentLocalVeg_50m, y = Biomass),
             size = 1.5) +
  scale_y_break(c(1.00,1.40), scales =0.25, 
                ticklabels = c(1.50, 3.00),
                space = 0.5)

# Biomass ~ PercentAg + Max Buffer Width + Season ------------------------------
# outliers removed

m1 <- cpglm(Biomass ~ PercentAg + MaxBufferWidth_m + Season, 
            link = "log", data = invert)

# create data frame
d <- expand.grid(MaxBufferWidth_m = seq(min(invert$MaxBufferWidth_m), 
                                           max(invert$MaxBufferWidth_m), 
                                           length.out = 1000),
                 PercentAg = mean(invert$PercentAg),
                 Season = unique(invert$Season))

ci <- trtools::lincon(m1, fcov=vcov)

### ...add confidence intervals ----

# Generate predictions on the log scale using type = "link"
d$yhat_log <- predict(m1, newdata = d, type = "link")

# Extract the covariance matrix from the model
vcov_matrix <- vcov(m1)

# Create the design matrix for the new data
X <- model.matrix(~ PercentAg + MaxBufferWidth_m + Season, data = d)

# Calculate the variance for each prediction (on the log scale)
pred_var <- diag(X %*% vcov_matrix %*% t(X))

# Calculate standard errors on the log scale
pred_se_log <- sqrt(pred_var)

# Calculate the 95% confidence intervals on the log scale
d$lower_CI_log <- d$yhat_log - 1.96 * pred_se_log
d$upper_CI_log <- d$yhat_log + 1.96 * pred_se_log

# Exponentiate the predictions and confidence intervals to get them 
# on the original scale
d$yhat <- exp(d$yhat_log)
d$lower_CI <- exp(d$lower_CI_log)
d$upper_CI <- exp(d$upper_CI_log)

# Plot predictions with confidence intervals
# Biomass ~ Ag + Season Plot ----
ggplot(d, aes(x = MaxBufferWidth_m, y = yhat, col = Season)) +  
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = Season), 
              alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "Maximum Buffer Width (m)", 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) +
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = MaxBufferWidth_m, y = Biomass),
             size = 1.5) +
  scale_y_break(c(1.60,2.00), scales =0.10, 
                ticklabels = c(1.75, 3.00),
                space = 0.5)

summary(m1)

# Biomass ~ PercentAg + LocalVeg + Dist to Wetland + MaxBufferWidth + Season----
# outliers removed
invert <- subset(invert, Biomass < 1.4)

# invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")
m1 <- cpglm(Biomass ~ PercentLocalVeg_50m + PercentAg +
            Dist_Closest_Wetland_m + MaxBufferWidth_m + Season, 
            link = "log", data = invert)

# create data frame
d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                           max(invert$PercentAg), 
                                           length.out = 1000),
                 MaxBufferWidth_m = mean(invert$MaxBufferWidth_m),
                 PercentLocalVeg_50m = mean(invert$PercentLocalVeg_50m),
                 Dist_Closest_Wetland_m = mean(invert$Dist_Closest_Wetland_m),
                 Season = unique(invert$Season))

ci <- trtools::lincon(m1, fcov=vcov)

### ...add confidence intervals ----

# Generate predictions on the log scale using type = "link"
d$yhat_log <- predict(m1, newdata = d, type = "link")

# Extract the covariance matrix from the model
vcov_matrix <- vcov(m1)

# Create the design matrix for the new data
X <- model.matrix(~ PercentLocalVeg_50m + PercentAg +
                    Dist_Closest_Wetland_m + MaxBufferWidth_m + Season, 
                  data = d)

# Calculate the variance for each prediction (on the log scale)
pred_var <- diag(X %*% vcov_matrix %*% t(X))

# Calculate standard errors on the log scale
pred_se_log <- sqrt(pred_var)

# Calculate the 95% confidence intervals on the log scale
d$lower_CI_log <- d$yhat_log - 1.96 * pred_se_log
d$upper_CI_log <- d$yhat_log + 1.96 * pred_se_log

# Exponentiate the predictions and confidence intervals to get them 
# on the original scale
d$yhat <- exp(d$yhat_log)
d$lower_CI <- exp(d$lower_CI_log)
d$upper_CI <- exp(d$upper_CI_log)

# Plot predictions with confidence intervals
# Biomass ~ Ag + Season Plot ----
ggplot(d, aes(x = PercentAg, y = yhat, col = Season)) +  
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = Season), 
              alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) +
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = PercentAg, y = Biomass),
             size = 1.5)

