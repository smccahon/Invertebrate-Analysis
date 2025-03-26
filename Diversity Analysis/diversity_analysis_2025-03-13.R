#------------------------------#
#   Order Diversity Analysis   #
#     Created 2025-03-06       #
#    Modified 2025-03-24       #
#------------------------------#

# Load packages
library(ggplot2)
library(dplyr)
library(car)
library(AICcmodavg)
library(glmmTMB) 
library(MASS)
library(MuMIn)

# Read in data
setwd("processed")
invert <- read.csv("Macroinverterbrate_Analysis_2025-03-06.csv")

# Transform data
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


# Standardize continuous variables
invert.cs <- invert
invert.cs$PercentAg <- scale(invert.cs$PercentAg)
invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
invert.cs$pH <- scale(invert.cs$pH)
invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)

# Models
invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")

### Agricultural Models ####
summary(zinb_model)
m1 <- glm.nb(Diversity ~ PercentAg + Season, data = invert.cs, link = "identity")

### Vegetation Models ####
m2 <- glm.nb(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs)
m3 <- glm.nb(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs)
m4 <- glm.nb(Diversity ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
         data = invert.cs)

### Water Quality Models #### 
m5 <- glm.nb(Diversity ~ pH + Season, data = invert.cs)
m6 <- glm.nb(Diversity ~ TDS_mg.L + Season, data = invert.cs)
m7 <- glm.nb(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Model ####
m8 <- glm.nb(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs)

### Informed Null Model ####
m9 <- glm.nb(Diversity ~ Season, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:9)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Season is an issue due to skewed data...use quasi-Poisson distribution

# Extract log likelihood
m1 <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = poisson)
m2 <- glm(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs, family = poisson)
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs, family = poisson)
m4 <- glm(Diversity ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
             data = invert.cs, family = poisson)
m5 <- glm(Diversity ~ pH + Season, data = invert.cs, family = poisson)
m6 <- glm(Diversity ~ TDS_mg.L + Season, data = invert.cs, family = poisson)
m7 <- glm(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs, family = poisson)
m8 <- glm(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs, family = poisson)
m9 <- glm(Diversity ~ Season, data = invert.cs, family = poisson)

model_names <- paste0("m", 1:9)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Manual
chat <- deviance(m1)/df.residual(m1)
logLik <- logLik(m1)
-2*logLik/chat + 2*3*(77/(77-3-1))

# Extract log-likelihood, number of parameters, and sample size
logLik_m1 <- logLik(m1)
k <- length(coef(m1))  # Number of parameters
n <- nrow(invert.cs)   # Number of observations
chat <- deviance(m1) / df.residual(m1)  # Dispersion parameter

# Calculate QAIC and QAICc
QAIC_m1 <- -2 * logLik_m1 / chat + 2 * k
QAICc_m1 <- QAIC_m1 + (2 * k * (k + 1)) / (n - k - 1)

# For loop ####
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
min_QAICc <- min(results$QAICc, na.rm = TRUE)  # Get the best (lowest) QAICc value
results$Delta_QAICc <- results$QAICc - min_QAICc  # Delta QAICc for each model

# Calculate weights: exp(-0.5 * delta_QAICc) / sum(exp(-0.5 * delta_QAICc))
results$Weight <- exp(-0.5 * results$Delta_QAICc) / sum(exp(-0.5 * results$Delta_QAICc))

# Print the updated results
results <- results %>% 
  arrange(QAICc) %>% 
  print

confint(m3)

# Modeling with quasipoisson ####
m1 <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = quasipoisson)
m2 <- glm(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs, family = quasipoisson)
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs, family = quasipoisson)
m4 <- glm(Diversity ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
          data = invert.cs, family = quasipoisson)
m5 <- glm(Diversity ~ pH + Season, data = invert.cs, family = quasipoisson)
m6 <- glm(Diversity ~ TDS_mg.L + Season, data = invert.cs, family = quasipoisson)
m7 <- glm(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs, family = quasipoisson)
m8 <- glm(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs, family = quasipoisson)
m9 <- glm(Diversity ~ Season, data = invert.cs, family = quasipoisson)


confint(m2)
confint(m7)
confint(m6)
confint(m5)
confint(m9)
confint(m8)
confint(m1)
confint(m4)
confint(m3)

model_names <- paste0("m", 1:9)
models <- mget(model_names)
model_avg <- model.avg(models)
summary(model_avg)


# Modeling with Poisson ####
m1 <- glm(Diversity ~ PercentAg + Season + -1, data = invert.cs, family = poisson)
m2 <- glm(Diversity ~ PercentBufferAroundWetland + Season, data = invert.cs, family = poisson)
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs, family = poisson)
m4 <- glm(Diversity ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
          data = invert.cs, family = poisson)
m5 <- glm(Diversity ~ pH + Season, data = invert.cs, family = poisson)
m6 <- glm(Diversity ~ TDS_mg.L + Season, data = invert.cs, family = poisson)
m7 <- glm(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs, family = poisson)
m8 <- glm(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs, family = poisson)
m9 <- glm(Diversity ~ Season, data = invert.cs, family = poisson)

model_names <- paste0("m", 1:9)
aictab(models, modnames = model_names)

summary(m1)

# ...model average ----
models <- mget(model_names)
model_avg <- model.avg(models)


summary(model_avg)

# 