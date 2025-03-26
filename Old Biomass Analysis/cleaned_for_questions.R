#------------------------------#
#   Order Diversity Analysis   #
#     Created 2025-03-19       #
#    Modified 2025-03-24       #
#------------------------------#

# Load packages & data ----
library(tidyverse)
library(AICcmodavg)
library(glmmTMB) # my data set has 14% zeros
library(trtools)
library(MASS)
library(pscl)

invert <- read.csv("processed/Macroinverterbrate_Analysis_2025-03-06.csv")
invert <- read.csv("Macroinverterbrate_Analysis_2025-03-06.csv")
 
# Data Manipulation & Visualization ----

# Data Distribution
hist(invert$Diversity, breaks = 20)


ggplot(invert, aes(x = Diversity, fill = Season)) + 
  geom_histogram(binwidth = 0.5, position = "stack", alpha = 0.5,
                 col = "black") +
  theme_classic() +
  labs(x = "Macroinvertebrate Diversity (# Orders)", 
       y = "Count") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen")) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(breaks = seq(0, 13, by = 2))
  


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
invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
invert.cs$pH <- scale(invert.cs$pH)
invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)

# Subset spring and fall
spring <- subset(invert.cs, Season == "Spring")
fall <- subset(invert.cs, Season == "Fall")

# Modeling ----

# ...Poisson ----
m <- glm(Diversity ~ Season + PercentAg, data = invert, family = 'poisson')


m <- lm(Biomass ~ PercentAg, data = invert)
# identity doesn't work, but log might?
m <- glm(Biomass ~ PercentAg, 
         family = quasi(link = "identity", variance = "mu^2"), 
         data = invert, start = c(0.27,0))
plot(predict(m),rstudent(m))

m <- glm(Biomass ~ PercentAg, 
         family = quasi(link = "log", variance = "mu^2"), 
         data = invert)
summary(m)
# summary
summary(m)
plot(predict(m),rstudent(m))

# dispersion parameter
deviance(m) / df.residual(m) # 1.38

# plotting
d <- expand.grid(Season = c("Spring",
                            "Fall"),
                 PercentAg = seq(min(invert$PercentAg),
                                 max(invert$PercentAg),
                                 length = 1000))

# type = "response" is required
d$yhat <- predict(m, newdata = d, type = "response")

# add confidence intervals
d <- cbind(d, glmint(m, newdata = d))

ggplot(invert, aes(x = PercentAg, y = Diversity, col = Season)) +
  geom_point(size = 3) +  
  geom_line(aes(y = yhat), data = d, size = 1) +
  geom_ribbon(aes(ymin = low, ymax = upp, y = NULL, fill = Season),
              color = NA, alpha = 0.25, data = d) + 
  theme_classic() +
  labs(x = "% Surrounding Agriculture", 
       y = "Macroinverterbate Diversity (# Orders)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) +
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))
### ...residuals ----

# Poisson
m.res <- resid(m)
plot(fitted(m), m.res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')
abline(0,0)

# quasi-Poisson ----
m1 <- glm(Diversity ~ Season + PercentAg, data = invert, family = 'quasipoisson')

# summary
summary(m1)

### ...residuals ----
m1.res <- resid(m1)
plot(fitted(m1), m1.res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Quasi-Poisson')
abline(0,0)

### ...QAICc ----

# Extract log likelihood
m1 <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = poisson)
m2 <- glm(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs, 
          family = poisson)
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs, 
          family = poisson)
m4 <- glm(Diversity ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
          data = invert.cs, family = poisson)
m5 <- glm(Diversity ~ pH + Season, data = invert.cs, family = poisson)
m6 <- glm(Diversity ~ TDS_mg.L + Season, data = invert.cs, family = poisson)
m7 <- glm(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs, family = poisson)
m8 <- glm(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs, 
          family = poisson)
m9 <- glm(Diversity ~ Season, data = invert.cs, family = poisson)

model_names <- paste0("m", 1:9)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m3)

# Manual AICc calculation
# Extract log-likelihood, number of parameters, and sample size
logLik_m1 <- logLik(m1)
k <- length(coef(m1))  # Number of parameters
n <- nrow(invert.cs)   # Number of observations
chat <- deviance(m1) / df.residual(m1)  # Dispersion parameter

# Calculate QAIC and QAICc
QAIC_m1 <- -2 * logLik_m1 / chat + 2 * k
QAICc_m1 <- QAIC_m1 + (2 * k * (k + 1)) / (n - k - 1)

# For loop 
# Initialize a data frame to store results
results <- data.frame(Model = paste0("m", 1:9), 
                      QAIC = NA, 
                      QAICc = NA, 
                      Delta_QAICc = NA, 
                      Weight = NA, 
                      LogLik = NA, 
                      k = NA)

# Loop through models m1 to m9
for (i in 1:9) {
  # Extract the model object
  model <- get(paste0("m", i))  # Retrieves the model object by its name (m1 to m9)
  
  # Extract log-likelihood, number of parameters, and sample size
  logLik_m <- logLik(model)
  k <- length(coef(model))  # Number of parameters
  n <- nrow(invert.cs)  # Number of observations
  chat <- deviance(model) / df.residual(model)  # Dispersion parameter
  
  # Calculate QAIC and QAICc
  QAIC_m <- -2 * logLik_m / chat + 2 * k
  QAICc_m <- QAIC_m + (2 * k * (k + 1)) / (n - k - 1)
  
  # Store the results in the data frame
  results$QAIC[i] <- QAIC_m
  results$QAICc[i] <- QAICc_m
  results$LogLik[i] <- logLik_m
  results$k[i] <- k
}

# Calculate delta QAICc and weights
min_QAICc <- min(results$QAICc, na.rm = TRUE) 
results$Delta_QAICc <- results$QAICc - min_QAICc 
results$Weight <- exp(-0.5 * results$Delta_QAICc) / 
  sum(exp(-0.5 * results$Delta_QAICc))

# Print
results <- results %>% 
  arrange(QAICc) %>% 
  print

# Manual QAICc calculation to double check (m1 as an example)
# Extract log-likelihood, number of parameters, and sample size
logLik_m1 <- logLik(m1)
k <- length(coef(m1))  # Number of parameters
n <- nrow(invert.cs)   # Number of observations
chat <- deviance(m1) / df.residual(m1)  # Dispersion parameter

# Calculate QAIC and QAICc
QAIC_m1 <- -2 * logLik_m1 / chat + 2 * k
QAICc_m1 <- QAIC_m1 + (2 * k * (k + 1)) / (n - k - 1)

# manual and computer calculation are equivalent

# quasi-poisson
m1 <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = quasipoisson)
m2 <- glm(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs, 
          family = quasipoisson)
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs, 
          family = quasipoisson)
m4 <- glm(Diversity ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
          data = invert.cs, family = quasipoisson)
m5 <- glm(Diversity ~ pH + Season, data = invert.cs, family = quasipoisson)
m6 <- glm(Diversity ~ TDS_mg.L + Season, data = invert.cs, family = quasipoisson)
m7 <- glm(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs, family = quasipoisson)
m8 <- glm(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs, 
          family = quasipoisson)
m9 <- glm(Diversity ~ Season, data = invert.cs, family = quasipoisson)

# ...model summaries ----
options(digits = 3)
cbind(summary(m3)$coefficients, confint(m3))



# Biomass Analysis ----
invert <- invert %>% 
  mutate(LogBiomass = log(invert$Biomass + 0.001))

hist(invert$LogBiomass, breaks = 20,
     xlab = "ln(Biomass + 0.001)",
     ylab = "Frequency",
     main = NULL)

# Negative Binomial ----


### Agricultural Models ####
m1 <- glm.nb(Diversity ~ PercentAg + Season, data = invert.cs)

### Vegetation Models ####
m2 <- glm.nb(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs)
m3 <- glm.nb(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs)
m4 <- glm.nb(Diversity ~ PercentBufferAroundWetland + PercentLocalVeg_50m + Season, 
         data = invert.cs)

### Water Quality Models #### 
m5 <- glm.nb(Diversity ~ pH + Season, data = invert.cs)
m6 <- glm.nb(Diversity ~ TDS_mg.L + Season, data = invert.cs)
m7 <- glm.nb(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Model ####
m8 <- glm.nb(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs)

### Informed Null
m9 <- glm.nb(Diversity ~ Season + PercentAg, data = invert.cs)

## ...by season ----

### Agricultural Models ####
m1 <- glm.nb(Diversity ~ PercentAg, data = fall)

### Vegetation Models ####
m2 <- glm.nb(Diversity ~ PercentBufferAroundWetland, data = fall)
m3 <- glm.nb(Diversity ~ PercentLocalVeg_50m, data = fall)
m4 <- glm.nb(Diversity ~ PercentBufferAroundWetland + PercentLocalVeg_50m, 
             data = fall)

### Water Quality Models #### 
m5 <- glm.nb(Diversity ~ pH, data = fall)
m6 <- glm.nb(Diversity ~ TDS_mg.L, data = fall)
m7 <- glm.nb(Diversity ~ pH + TDS_mg.L, data = fall)

### Hydroperiod Model ####
m8 <- glm.nb(Diversity ~ Dist_Closest_Wetland_m, data = fall)

### Null
m9 <- glm.nb(Diversity ~ 1, data = fall)


model_names <- paste0("m", 1:9)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)
summary(m1)

plot(resid(m1))

boxplot(fall$Diversity)
boxplot(fall$PercentAg)

# Zero-inflated models ----

zinb_model <- zeroinfl(Diversity ~ PercentAg + Season | PercentAg + Season, 
                       data = invert.cs, 
                       dist = "negbin")

summary(zinb_model)

AIC(zinb_model, m1, mp)
BIC(zinb_model, m1)


# ...versus Poisson
mp <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = "poisson")

# Compound Poisson-gamma model ----
library(cplm)
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", 
           data = invert.cs)
summary(m1)
plot(resid(m1))

# residual and qq plot
parold <- par(mfrow = c(2, 2), mar = c(5, 5, 2, 1))
# 1. regular plot
r1 <- resid(m1) / sqrt(m1$phi)
plot(r1 ~ fitted(m1), cex = 0.5)
qqnorm(r1, cex = 0.5)



m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", data = invert.cs)
m2 <- cpglm(Biomass ~ PercentBufferAroundWetland + Season, link = "log", 
            data = invert.cs)
m3 <- cpglm(Biomass ~ PercentLocalVeg_50m + Season, link = "log", 
            data = invert.cs)
m4 <- cpglm(Biomass ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
            link = "log", 
            data = invert.cs)
m5 <- cpglm(Biomass ~ pH + Season, link = "log", data = invert.cs)
m6 <- cpglm(Biomass ~ TDS_mg.L + Season, link = "log", data = invert.cs)
m7 <- cpglm(Biomass ~ pH + TDS_mg.L + Season, link = "log", data = invert.cs)
m8 <- cpglm(Biomass ~ Dist_Closest_Wetland_m + Season, link = "log", 
            data = invert.cs)
m9 <- cpglm(Biomass ~ Season, link = "log", data = invert.cs)

# ...AIC ----
library(MuMIn)
models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
model_comparison <- model.sel(models)
print(model_comparison)

# Percent Ag and Season have an extremely strong effect

summary(m6)


predictions <- predict(m1, newdata = invert.cs, type = "response", 
                       se.fit = TRUE)

# Calculate the confidence intervals for predictions
lower_bound <- predictions$fit - 1.96 * predictions$se.fit
upper_bound <- predictions$fit + 1.96 * predictions$se.fit

# Combine the predictions with the confidence intervals
predictions_with_CI <- data.frame(
  predicted = predictions$fit,
  lower = lower_bound,
  upper = upper_bound
)


d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                 max(invert$PercentAg), 
                                 length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m1, newdata = d, se.fit = TRUE, type = "response")

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentAg, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Agriculture Surrounding Wetland", 
       y = "Macroinvertebrate Diversity (# Orders)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen")) +
  geom_point(data = invert, aes(x = PercentAg, y = Diversity,
                                col = Season), size = 2)