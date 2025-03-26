#------------------------------#
#    Invert Biomass Analysis   #
#     Created 2025-03-24       #
#    Modified 2025-03-24       #
#------------------------------#

# Load packages & data ----
library(tidyverse)
library(AICcmodavg)
library(trtools)
library(MASS)
library(pscl)
library(MuMIn)
library(cplm)

invert <- read.csv("processed/Macroinverterbrate_Analysis_2025-03-06.csv")
invert <- read.csv("Macroinverterbrate_Analysis_2025-03-06.csv")

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


# Compound Poisson-gamma model ----
m1 <- cpglm(Biomass ~ PercentAg + Season, link = "log", 
            data = invert.cs)
summary(m1)
plot(resid(m1))

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
