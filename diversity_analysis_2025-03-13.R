#------------------------------#
#   Order Diversity Analysis   #
#     Created 2025-03-06       #
#    Modified 2025-03-13       #
#------------------------------#

# Load packages
library(ggplot2)
library(dplyr)
library(car)
library(AICcmodavg)
library(glmmTMB) # my data set has 14% zeros
library(MASS)

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
m1 <- glm.nb(Diversity ~ PercentAg + Season, data = invert.cs)

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

# If season is still an issue, maybe do poisson regression with a robust standard error (quasipoisson)

model_quasi <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = quasipoisson)
summary(model_quasi)
confint(model_quasi)
