#------------------------------#
#   Biomass & Neonic Analysis  #
#     Created 2025-03-17       #
#    Modified 2025-03-17       #
#------------------------------#

# Load packages
library(ggplot2)
library(dplyr)
library(car)
library(AICcmodavg)
library(glmmTMB) # my data set has 14% zeros
library(MuMIn)
library(data.table)

# Read in data
setwd("processed")
invert <- read.csv("Macroinverterbrate_Analysis_2025-03-06.csv")


invert <- invert %>% 
  mutate(LogBiomass = log(invert$Biomass + 0.001))

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
invert.cs$NeonicInvert_ng.g <- scale(invert.cs$NeonicInvert_ng.g)
invert.cs$NeonicWater_ng.L <- scale(invert.cs$NeonicWater_ng.L)
invert.cs$PesticideInvert_ng.g <- scale(invert.cs$PesticideInvert_ng.g)

# Correlations ####
# subset data
sample <- invert[, c("Season",
                     "PercentAg",
                     "TDS_mg.L",
                     "Permanence",
                     "NearestCropDistance_m",
                     "MaxBufferWidth_m",
                     "PercentLocalVeg_50m",
                     "AgCategory",
                     "pH_probe",
                     "DominantCrop",
                     "PercentBufferAroundWetland",
                     "Dist_Closest_Wetland_m",
                     "NeonicWater_ng.L",
                     "WaterDetection")]

sample.reduced <- invert[, c("Season",
                             "PercentAg",
                             "TDS_mg.L",
                             "Permanence",
                             "NearestCropDistance_m",
                             "MaxBufferWidth_m",
                             "PercentLocalVeg_50m",
                             "AgCategory",
                             "pH_probe",
                             "DominantCrop",
                             "PercentBufferAroundWetland",
                             "Dist_Closest_Wetland_m",
                             "NeonicInvert_ng.g",
                             "NeonicWater_ng.L",
                             "PesticideInvert_ng.g",
                             "WaterDetection",
                             "InvertPesticideDetection")]

# Convert categorical to numeric ####
sample$Season <- as.numeric(sample$Season)
sample$Permanence <- as.numeric(sample$Permanence)
sample$AgCategory <- as.numeric(sample$AgCategory)
sample$DominantCrop <- as.numeric(sample$DominantCrop)
sample$WaterDetection <- ifelse(sample$WaterDetection == "Y", 1, 
                                ifelse(sample$WaterDetection == "N", 0, NA))

sample.reduced$Season <- as.numeric(sample.reduced$Season)
sample.reduced$Permanence <- as.numeric(sample.reduced$Permanence)
sample.reduced$AgCategory <- as.numeric(sample.reduced$AgCategory)
sample.reduced$DominantCrop <- as.numeric(sample.reduced$DominantCrop)
sample.reduced$WaterDetection <- ifelse(sample.reduced$WaterDetection == "Y", 1, 
                                ifelse(sample.reduced$WaterDetection == "N", 0, NA))
sample.reduced$InvertPesticideDetection <- ifelse(sample.reduced$InvertPesticideDetection == 
                                                    "Y", 1, 
                                ifelse(sample.reduced$InvertPesticideDetection == "N", 0, NA))

# Correlations with full dataset (n = 77)
cor(sample)

# Correlations with reduced dataset with invert detections (n = 51)
reduced <- sample.reduced[!is.na(sample.reduced$InvertPesticideDetection),]

cor(reduced)

### High Correlations > 0.6 ####
# Full dataset (water neonic)

### Permanence & Season
# Percent Ag & Nearest Crop Distance
# Percent Ag & Ag Category
# Percent Ag & Dominant Crop
# Percent Buffer & Permanence
# Ag Category & Nearest Crop Distance
# Percent Buffer & Buffer Width

### Reduced dataset (invert neonic)
# Water detection & water neonic conc
# Percent Ag & Nearest Crop Distance
# Percent Ag & Ag Category
# Percent Ag & Dominant Crop
# Percent Buffer & Buffer Width
# Ag Category + Distance to nearest crop


# Which neonic variables should I use? ####

### Reduced dataset ####
# *if you want to use full dataset you must use water detections


I LEFT OFF HERE, ADD LOGBIOMASS TO REDUCED DATASET




# Top Models from Previous Analysis ####
### Informed Null ####
invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")
m1 <- lm(LogBiomass ~ Season, data = invert.cs)

### Single Models ####
m2 <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)
m3 <- lm(LogBiomass ~ PercentLocalVeg_50m + Season, data = invert.cs)
m4 <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert.cs)
m5 <- lm(LogBiomass ~ NeonicWater_ng.L)

### Additive Models ####
m5 <- lm(LogBiomass ~ Season + PercentAg + PercentLocalVeg_50m, data = invert.cs) # ag and veg
m6 <- lm(LogBiomass ~ Season + PercentAg + Dist_Closest_Wetland_m, data  = invert.cs) # ag and hydro
m7 <- lm(LogBiomass ~ Season + PercentAg + PercentLocalVeg_50m + Dist_Closest_Wetland_m, 
         data = invert.cs) # ag, veg, and hydro
m8 <- lm(LogBiomass ~ Season + PercentLocalVeg_50m + Dist_Closest_Wetland_m, 
         data = invert.cs) # hydro + veg

# AIC MODEL SELECTION
model_names <- paste0("m", 1:8)

models <- mget(model_names)

aictab(models, modnames = model_names)

cor(invert$NeonicWater_ng.L, )













