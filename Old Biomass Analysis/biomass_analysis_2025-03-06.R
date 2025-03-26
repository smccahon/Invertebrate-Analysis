#------------------------------#
# Biomass & Diversity Analysis #
#     Created 2025-03-06       #
#    Modified 2025-03-14       #
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

# Data Transformation ####
# hist(invert$Biomass, breaks = 20) # left  skewed distribution
# hist(log(invert$Biomass + 0.001), breaks = 20) # normally distributed

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

# Subset spring and fall
spring <- subset(invert.cs, Season == "Spring")
fall <- subset(invert.cs, Season == "Fall")

# Data Exploration ####

## Date & Season ####
ggplot(invert, aes(x = Season, y = Biomass)) + geom_boxplot() +
  theme_light() + 
  labs(x = "Season",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  ylim(0, 0.7)


ggplot(invert, aes(x = Julian, y = Biomass)) + geom_point() +
  theme_light() + 
  labs(x = "Julian Date",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  ylim(0, 0.7)

## Ag & Biomass ####
ggplot(invert, aes(x = SPEI, y = Biomass, col = Season)) + 
  geom_point(size = 2.5) +
  theme_classic() + 
  labs(x = "Drought Condition Index (SPEI)",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))

ggplot(invert, aes(x = NearestCropDistance_m, y = Biomass, col = Season)) + 
  geom_point(size = 2.5) +
  theme_classic() + 
  labs(x = "Nearest Crop Distance (m)",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))

ggplot(invert, aes(x = DominantCrop, y = Biomass, fill = Season)) + 
  geom_boxplot(alpha = 0.4) +
  theme_classic() + 
  labs(x = "Dominant Land Use Type",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
   ylim(0, 0.7) +
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))

ggplot(invert, aes(x = DominantCrop, y = Biomass)) + 
  geom_boxplot(alpha = 0.4) +
  theme_classic() + 
  labs(x = "Dominant Land Use Type",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))

ggplot(invert, aes(x = AgCategory, y = Biomass, fill = Season)) + 
  geom_boxplot(alpha = 0.4) +
  theme_classic() + 
  labs(x = "Agricultural Category",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen")) +
  scale_x_discrete(labels = c("Low" = "Low (< 25%)",
                              "Moderate" = "Moderate (25 - 75%)",
                              "High" = "High (> 75%)"))

ggplot(invert, aes(x = AgCategory, y = Biomass)) + 
  geom_boxplot(alpha = 0.4) +
  theme_classic() + 
  labs(x = "Agricultural Category",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_x_discrete(labels = c("Low" = "Low (< 25%)",
                              "Moderate" = "Moderate (25 - 75%)",
                              "High" = "High (> 75%)"))

## Buffer & Biomass ####
ggplot(invert, aes(x = PercentBufferAroundWetland, y = Biomass, col = Season)) + 
  geom_point(size = 2.5) +
  theme_classic() + 
  labs(x = "Vegetation Buffer Cover (%)",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))

ggplot(invert, aes(x = MaxBufferWidth_m, y = Biomass, col = Season)) + 
  geom_point(size = 2.5) +
  theme_classic() + 
  labs(x = "Maximum Buffer Width at Wetland (m)",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.7) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))

## Permanence & Biomass ####
ggplot(invert, aes(x = Permanence, y = Biomass, fill = Season)) + 
  geom_boxplot(alpha=0.4) +
  theme_classic() + 
  labs(x = "Wetland Permanence",
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  ylim(0, 0.6) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen"))


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
                     "Dist_Closest_Wetland_m")]

# Convert categorical to numeric
sample$Season <- as.numeric(sample$Season)
sample$Permanence <- as.numeric(sample$Permanence)
sample$AgCategory <- as.numeric(sample$AgCategory)
sample$DominantCrop <- as.numeric(sample$DominantCrop)

# Season and permanence can't be in the same model
m <- glm(Season ~ Permanence, data = sample, family = "binomial")
summary(m)
confint(m)
anova(m)
kruskal.test(Season~Permanence, data = invert)

ggplot(sample, aes(x = factor(Season), y = Permanence)) + geom_boxplot()

# Dominant crop and percent agriculture can't be in the same model
m <- lm(PercentAg ~ DominantCrop, data = invert)
summary(m)
confint(m)
anova(m)
ggplot(invert, aes(x = DominantCrop, y = PercentAg)) + geom_boxplot()

# Agricultural category and nearest crop distance can't be in the same model
m <- lm(NearestCropDistance_m ~ AgCategory, data = invert)
summary(m)
confint(m)

# Percent Buffer & Permanence can't be in the same model
m <- lm(PercentBufferAroundWetland ~ Permanence, data = invert)
confint(m)
ggplot(invert, aes(x = Permanence, y = PercentBufferAroundWetland)) +
  geom_boxplot()

# Percent Buffer & Max Buffer Width can't be in the same model

cor(sample)

### High Correlations ####

# Season and permanence (-0.63)
# Distance to nearest crop and percent agriculture (-0.72)
# agricultural category and percent agriculture (0.89)
# Dominant crop and percent agriculture (0.63)
# agricultural category and nearest crop distance (-0.66)
# max buffer width and percent buffer (0.77)
# permanence and percent buffer (0.668)

### KATIE: If two variables are highly correlated, use AIC to tell you which one is better
### Can indicate all the variables you considered in methods, and say you dropped the one with less explanatory power
### Also use season as null
### Test the effect of permanence on biomass AFTER separately (another analysis given the correlation with season)

# Model selection (a priori) BIOMASS ####

### INFORMED NULL ####
m17 <- lm(LogBiomass ~ Season, data = invert.cs)

### Agricultural Models ####

m1 <- lm(LogBiomass ~ PercentAg, data = invert.cs)
m2 <- lm(LogBiomass ~ NearestCropDistance_m, data = invert.cs)
m3 <- lm(LogBiomass ~ DominantCrop, data = invert.cs)
m4 <- lm(LogBiomass ~ AgCategory, data = invert.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Best agricultural model is percent agriculture

### Vegetation Models ####
m1 <- lm(LogBiomass ~ PercentBufferAroundWetland, data = invert.cs)
m2 <- lm(LogBiomass ~ MaxBufferWidth_m, data = invert.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Best is percent buffer around wetland

### Water Quality Models #### 
m12 <- lm(LogBiomass ~ pH + Season, data = invert.cs)
m13 <- lm(LogBiomass ~ TDS_mg.L + Season, data = invert.cs)
m14 <- lm(LogBiomass ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Models ####
m15 <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert.cs)

# Null Model
m16 <- lm(LogBiomass ~ 1, data = invert.cs)

#
m1 <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)
m2 <- lm(LogBiomass ~ Season, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

## Top Model Summaries ####
par(mfrow = c(2,2))

summary(m15)
confint(m4) # percent ag has a significant negative effect
vif(m4)
plot(m4)

summary(m13) # season has a significant effect
confint(m13)
vif(m13)
plot(m13)

summary(m1) # season has a significant effect
confint(m1)
vif(m1)
plot(m1)

summary(m7)
confint(m7) # season has a significant effect
vif(m7)
plot(m7)

summary(m12)
confint(m12) # season has a significant effect
vif(m12)
plot(m12)

summary(m15)
confint(m15) # season has a significant effect
vif(m15)

summary(m14)
confint(m14) # season has a significant effect
vif(m14)

summary(m6)
confint(m6) # maximum buffer width does not have a significant effect, season does
vif(m6)

summary(m5)
confint(m5) # season has a significant effect
vif(m5)

summary(m8)
confint(m8) #season has a significant effect
vif(m8)

summary(m10)
confint(m10) # season has a significant effect
vif(m10)

summary(m2)
confint(m2) # season has a significant effect
vif(m2)

summary(m9)
confint(m9) # season has a significant effect
vif(m9)

summary(m3)
confint(m3) # season has a significant effect

summary(m11)
confint(m11) # season has a significant effect

# Graph Results ####
## LOCAL VEG MODEL ####
m <- lm(LogBiomass ~ PercentLocalVeg_50m + Season, data = invert)

d <- expand.grid(PercentLocalVeg_50m = seq(min(invert$PercentLocalVeg_50m), 
                                           max(invert$PercentLocalVeg_50m), 
                                           length.out = 1000),  # Correct sequence for PercentLocalVeg_50m
                 Season = unique(invert$Season))

# d$PercentLocalVeg_50m <- scale(d$PercentLocalVeg_50m, 
#                                center = mean(invert.cs$PercentLocalVeg_50m), 
#                                scale = sd(invert.cs$PercentLocalVeg_50m))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentLocalVeg_50m, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "Local Vegetation Cover (%)", 
       y = "Log(Macroinvertebrate Biomass)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.1,0.9)) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen")) +
  geom_point(data = invert, aes(x = PercentLocalVeg_50m, y = LogBiomass,
                         col = Season), size = 2)


## DISTANCE TO WETLAND MODEL ####
m <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert)

d <- expand.grid(Dist_Closest_Wetland_m = seq(min(invert$Dist_Closest_Wetland_m), 
                                           max(invert$Dist_Closest_Wetland_m), 
                                           length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = Dist_Closest_Wetland_m, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "Distance to Nearest Wetland (m)", 
       y = "Log(Macroinvertebrate Biomass)") +
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
  geom_point(data = invert, aes(x = Dist_Closest_Wetland_m, y = LogBiomass,
                                col = Season), size = 2)

## PERCENT AGRICULTURE MODEL ####
m <- lm(LogBiomass ~ PercentAg + Season + NearestCropDistance_m, data = invert.cs)

d <- expand.grid(NearestCropDistance_m = seq(min(invert$NearestCropDistance_m), 
                                              max(invert$NearestCropDistance_m), 
                                              length.out = 1000),  
                 Season = unique(invert$Season),
                 PercentAg = mean(invert$PercentAg))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = NearestCropDistance_m, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "Distance to Nearest Crop (m)", 
       y = "Log(Macroinvertebrate Biomass)") +
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
  geom_point(data = invert, aes(x = NearestCropDistance_m, y = LogBiomass,
                                col = Season), size = 2)

## AG CATEGORY MODELS ####
m <- lm(LogBiomass ~ AgCategory + Season, data = invert.cs)

d <- expand.grid(AgCategory = unique(invert$AgCategory), 
                 Season = c("Spring"))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lwr <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upr <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = AgCategory, y = predicted_Mass)) +
  geom_point(size = 3) +  
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 1) +
  theme_classic() +
  labs(x = "Agricultural Category", 
       y = "Model Estimated Log(Biomass)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_x_discrete(labels = c("Low" = "Low (< 25%)",
                              "Moderate" = "Moderate (25 - 75%)",
                              "High" = "High (> 75%)"))


# Model selection with interactions ####

m1 <- lm(LogBiomass ~ Season * Permanence, data = invert.cs)
m2 <- lm(LogBiomass ~ Season + Permanence, data = invert.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)


### Agricultural Models ####
m1 <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)
m2 <- lm(LogBiomass ~ NearestCropDistance_m + Season, data = invert.cs)
m3 <- lm(LogBiomass ~ DominantCrop + Season, data = invert.cs)
m4 <- lm(LogBiomass ~ PercentAg + NearestCropDistance_m + Season, data = invert.cs)

### Vegetation Models ####
m5 <- lm(LogBiomass ~ PercentBufferAroundWetland + Season, data = invert.cs)
m6 <- lm(LogBiomass ~ MaxBufferWidth_m + Season, data = invert.cs)
m7 <- lm(LogBiomass ~ PercentLocalVeg_50m + Season, data = invert.cs)
m8 <- lm(LogBiomass ~ PercentBufferAroundWetland + Season +
           MaxBufferWidth_m + PercentLocalVeg_50m, data = invert.cs)

### Water Quality Models #### 
m9 <- lm(LogBiomass ~ pH + Season, data = invert.cs)
m10 <- lm(LogBiomass ~ TDS_mg.L + Season, data = invert.cs)
m11 <- lm(LogBiomass ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Models ####
m12 <- lm(LogBiomass ~ Permanence * Season, data = invert.cs)
m13 <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert.cs)
m14 <- lm(LogBiomass ~ SPEI + Season, data = invert.cs)
m15 <- lm(LogBiomass ~ Permanence * Season + Dist_Closest_Wetland_m + SPEI,
          data = invert.cs)

# extras
m16 <- lm(LogBiomass ~ Permanence + Season + Dist_Closest_Wetland_m + SPEI,
          data = invert.cs)
m17 <- lm(LogBiomass ~ Permanence * Season + Dist_Closest_Wetland_m + SPEI * Season,
          data = invert.cs)
m18 <- lm(LogBiomass ~ SPEI * Season, data = invert.cs)
m19 <- lm(LogBiomass ~ pH * Season, data = invert.cs)
m20 <- lm(LogBiomass ~ TDS_mg.L * Season, data = invert.cs)
m21 <- lm(LogBiomass ~ pH * Season + TDS_mg.L * Season, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:21)

models <- mget(model_names)

aictab(models, modnames = model_names)

### INTERACTIONS DID NOT IMPROVE INFERENCE ####

# Stratify analysis by season ####

### Agricultural Models ####
m1 <- lm(LogBiomass ~ PercentAg, data = spring)
m2 <- lm(LogBiomass ~ NearestCropDistance_m, data = spring)
m3 <- lm(LogBiomass ~ DominantCrop, data = spring)
m4 <- lm(LogBiomass ~ PercentAg + NearestCropDistance_m, data = spring)

### Vegetation Models ####
m5 <- lm(LogBiomass ~ PercentBufferAroundWetland, data = spring)
m6 <- lm(LogBiomass ~ MaxBufferWidth_m, data = spring)
m7 <- lm(LogBiomass ~ PercentLocalVeg_50m, data = spring)
m8 <- lm(LogBiomass ~ PercentBufferAroundWetland +
           MaxBufferWidth_m + PercentLocalVeg_50m, data = spring)

### Water Quality Models #### 
m9 <- lm(LogBiomass ~ pH, data = spring)
m10 <- lm(LogBiomass ~ TDS_mg.L, data = spring)
m11 <- lm(LogBiomass ~ pH + TDS_mg.L, data = spring)

### Hydroperiod Models ####
m12 <- lm(LogBiomass ~ Permanence, data = spring)
m13 <- lm(LogBiomass ~ Dist_Closest_Wetland_m, data = spring)
m14 <- lm(LogBiomass ~ SPEI, data = spring)
m15 <- lm(LogBiomass ~ Permanence + Dist_Closest_Wetland_m + SPEI,
          data = spring)
m16 <- lm(LogBiomass ~ 1, data = spring)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:16)

models <- mget(model_names)

aictab(models, modnames = model_names)

# For fall, biomass higher in permanent wetlands compared to temporary wetlands and lower in high ag wetlands
# For spring, no covariates are significant


# Model selection (a priori) DIVERSITY ####
## CORRELATIONS ####


### Agricultural Models ####
m1 <- lm(LogDiversity ~ PercentAg + Season, data = invert.cs)
m2 <- lm(LogDiversity ~ NearestCropDistance_m + Season, data = invert.cs)
m3 <- lm(LogDiversity ~ DominantCrop + Season, data = invert.cs)
m4 <- lm(LogDiversity ~ PercentAg + NearestCropDistance_m + Season, data = invert.cs)

### Vegetation Models ####
m5 <- lm(LogDiversity ~ PercentBufferAroundWetland + Season, data = invert.cs)
m6 <- lm(LogDiversity ~ MaxBufferWidth_m + Season, data = invert.cs)
m7 <- lm(LogDiversity ~ PercentLocalVeg_50m + Season, data = invert.cs)
m8 <- lm(LogDiversity ~ PercentBufferAroundWetland + Season +
           MaxBufferWidth_m + PercentLocalVeg_50m, data = invert.cs)

### Water Quality Models #### 
m9 <- lm(LogDiversity ~ pH + Season, data = invert.cs)
m10 <- lm(LogDiversity ~ TDS_mg.L + Season, data = invert.cs)
m11 <- lm(LogDiversity ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Models ####
m12 <- lm(LogDiversity ~ Permanence + Season, data = invert.cs)
m13 <- lm(LogDiversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs)
m14 <- lm(LogDiversity ~ SPEI + Season, data = invert.cs)
m15 <- lm(LogDiversity ~ Permanence + Dist_Closest_Wetland_m + SPEI +
            Season, data = invert.cs)

# Null Model
m16 <- lm(LogDiversity ~ 1, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:16)

models <- mget(model_names)

aictab(models, modnames = model_names)

hist(invert$Diversity)

# Top Model Summaries ####
par(mfrow = c(2,2))
summary(m7)
confint(m7) # local veg cover does not have a significant effect (season does)
plot(m7)

summary(m12)
confint(m12) # permanence does not have a significant effect (season does)

summary(m1)
confint(m1) # percent ag does not have a significant effect (season does)

summary(m13)
confint(m13) # distance to closest wetland does not have a significant effect (season does)

summary(m8)
confint(m8) # buffer width, local veg, and buffer cover do not have a significant effect (season does)

summary(m6)
confint(m6) # buffer width does not have a significant effect (season does)

summary(m4)
confint(m4) # percent ag and nearest crop distance do not have a significant effect (season does)

summary(m9)
confint(m9) # pH does not have a significant effect (season does)

summary(m14)
confint(m14) # SPEI does not have a significant effect (season does)

summary(m10) 
confint(m10) # TDS does not have a significant effect (season does)

summary(m2)
confint(m2) # nearest crop distance does not have a significant effect (season does)

summary(m5)
confint(m5) # buffer cover does not have a significant effect (season does)

summary(m11)
confint(m11) # pH and TDS do not have a significant effect (season does)

summary(m15)
confint(m15) # permanence, dist to closest wetland, and SPEI are not significant (Season is)

summary(m3)
confint(m3) # dominant crop is not significant (season is)

# Graph Results ####

### LOCAL VEG MODEL (BEST) DIVERSITY ####
m <- lm(Diversity ~ PercentLocalVeg_50m + Season, data = invert)

d <- expand.grid(PercentLocalVeg_50m = seq(min(invert$PercentLocalVeg_50m), 
                                           max(invert$PercentLocalVeg_50m), 
                                           length.out = 1000),  # Correct sequence for PercentLocalVeg_50m
                 Season = unique(invert$Season))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentLocalVeg_50m, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "Local Vegetation Cover (%)", 
       y = "Macroinvertebrate Diversity (# Orders)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.1,0.9)) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen")) +
  geom_point(data = invert, aes(x = PercentLocalVeg_50m, y = Diversity,
                                col = Season), size = 2)

### PERCENT AG MODEL DIVERSITY ####
m <- lm(Diversity ~ PercentAg + Season, data = invert)

d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                 max(invert$PercentAg), 
                                 length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m, newdata = d, se.fit = TRUE)

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


## AG CATEGORY MODELS ####
m <- lm(LogBiomass ~ AgCategory + Season, data = invert)

d <- expand.grid(AgCategory = unique(invert$AgCategory), 
                 Season = c("Spring"))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lwr <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upr <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = AgCategory, y = predicted_Mass)) +
  geom_point(size = 3) +  
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 1) +
  theme_classic() +
  labs(x = "Agricultural Category", 
       y = "Model Estimated Biomass") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_x_discrete(labels = c("Low" = "Low (< 25%)",
                              "Moderate" = "Moderate (25 - 75%)",
                              "High" = "High (> 75%)"))



## MODELS WITHOUT CORRELATIONS ####
cor(invert$NearestCropDistance_m, invert$PercentAg) # -0.724
cor(invert$PercentBufferAroundWetland, invert$MaxBufferWidth_m) # 0.777
cor(invert$PercentBufferAroundWetland, invert$PercentLocalVeg_50m) # -0.01
cor(invert$MaxBufferWidth_m, invert$PercentLocalVeg_50m) # -0.0271
cor(invert$pH, invert$TDS_mg.L) # 0.269
cor(invert$SPEI, invert$Dist_Closest_Wetland_m) # -0.0485

### Agricultural Models ####
m1 <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)
m2 <- lm(LogBiomass ~ NearestCropDistance_m + Season, data = invert.cs)
m3 <- lm(LogBiomass ~ DominantCrop + Season, data = invert.cs)

### Vegetation Models ####
m4 <- lm(LogBiomass ~ PercentBufferAroundWetland + Season, data = invert.cs)
m5 <- lm(LogBiomass ~ MaxBufferWidth_m + Season, data = invert.cs)
m6 <- lm(LogBiomass ~ PercentLocalVeg_50m + Season, data = invert.cs)
m7 <- lm(LogBiomass ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
         data = invert.cs)
m8 <- lm(LogBiomass ~ MaxBufferWidth_m + Season + PercentLocalVeg_50m, 
         data = invert.cs)

### Water Quality Models #### 
m9 <- lm(LogBiomass ~ pH + Season, data = invert.cs)
m10 <- lm(LogBiomass ~ TDS_mg.L + Season, data = invert.cs)
m11 <- lm(LogBiomass ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Models ####
m12 <- lm(LogBiomass ~ Permanence + Season, data = invert.cs)
m13 <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert.cs)
m14 <- lm(LogBiomass ~ SPEI + Season, data = invert.cs)
m15 <- lm(LogBiomass ~ Permanence + Dist_Closest_Wetland_m + SPEI +
            Season, data = invert.cs)

# Null Model
m16 <- lm(LogBiomass ~ 1, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:16)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m11)


m <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)

d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                 max(invert$PercentAg), 
                                 length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentAg, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Agriculture Surrounding Wetland", 
       y = "Log(Macroinvertebrate Biomass)") +
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
  geom_point(data = invert, aes(x = PercentAg, y = LogBiomass,
                                col = Season), size = 2)
confint(m)

summary(m)

## PCA
# Check correlation
cor(invert$PercentAg, invert$NearestCropDistance_m)

# Scatter plot to visualize the relationship
ggplot(invert, aes(x = PercentAg, y = NearestCropDistance_m, col= Season)) + 
  geom_point(size = 2) +
  theme_classic() + 
  labs(x = "% Agriculture Surrounding Wetland",
       y = "Distance to Nearest Crop (m)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 13)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 13)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen"))

m <- lm(LogBiomass ~ AgCategory + Season, data = invert.cs)
confint(m)


# RERUN ANALYSIS WITH SEASON AS INFORMED NULL ####
invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")

### Agricultural Models ####
m1 <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)

### Vegetation Models ####
m2 <- lm(LogBiomass ~ PercentBufferAroundWetland + Season, data = invert.cs)
m3 <- lm(LogBiomass ~ PercentLocalVeg_50m + Season, data = invert.cs)
m4 <- lm(LogBiomass ~ PercentBufferAroundWetland + Season + PercentLocalVeg_50m, 
         data = invert.cs)

### Water Quality Models #### 
m5 <- lm(LogBiomass ~ pH + Season, data = invert.cs)
m6 <- lm(LogBiomass ~ TDS_mg.L + Season, data = invert.cs)
m7 <- lm(LogBiomass ~ pH + TDS_mg.L + Season, data = invert.cs)

### Hydroperiod Model ####
m8 <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert.cs)

### Informed Null Model ####
m9 <- lm(LogBiomass ~ Season, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:9)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m8)
summary(m8)
confint(m1)
summary(m1)
confint(m3)
summary(m3)
summary(m4)
confint(m4)
summary(m9)
confint(m9)
summary(m2)
confint(m2)
summary(m6)
confint(m6)
summary(m5)
confint(m5)
summary(m7)
confint(m7)

par(mfrow = c(2,2))
plot(m2)

# MODEL AVERAGE ####
model_avg <- model.avg(models)
summary(model_avg)
confint(model_avg)

# Plot model averaged results ####
ma <- summary(model_avg)
df <- as.data.frame(ma$coefmat.full)
ci <- as.data.frame(confint(model_avg, full = T))
df$low <-ci$`2.5 %`
df$upp <-ci$`97.5 %`
setDT(df, keep.rownames = "coefficient")
names(df) <- gsub(" ", "", names(df))
df

# add spring data
spring <- data.table(
  "coefficient" = "SeasonSpring",  
  "Estimate" = -2.117356046 ,         
  "Std.Error" = 0.42804650 ,       
  "AdjustedSE" = 0.43458554 ,        
  "zvalue" = 4.87212719 ,         
  "Pr(>|z|)" = 0.0000011 ,             
  "low" = -2.96912806 ,        
  "upp" = -1.26558403               
)

# Add the "SeasonSpring" row to the original df
df <- rbind(df, spring)

ggplot(data = df[2:9,], aes(x = coefficient, y = Estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", lwd = 1.5) + 
  geom_errorbar(aes(ymin = Estimate - AdjustedSE, ymax = Estimate + AdjustedSE), 
                color="black", width = 0, lwd = 1.5) +
  coord_flip() +  # flipping x and y axes
  geom_point(size=4) +
  theme_classic() +
  labs(y = "Model-Averaged Parameter Estimate", 
       x = NULL) +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.9)) +
  scale_x_discrete(limits = c("pH", "TDS_mg.L", "PercentBufferAroundWetland", 
                              "PercentLocalVeg_50m", "PercentAg", "Dist_Closest_Wetland_m", 
                              "SeasonFall", "SeasonSpring"),
                   labels = c("pH" = "pH",
                              "TDS_mg.L" = "Total Dissolved Solids",
                              "PercentBufferAroundWetland" = "Buffer Cover",
                              "PercentLocalVeg_50m" = "Local Vegetation Cover", 
                              "PercentAg" = "% Surrounding Ag", 
                              "Dist_Closest_Wetland_m" = "Dist. to Closest Wetland",
                              "SeasonFall" = "Season (Fall)", 
                              "SeasonSpring" = "Season (Spring)"))


## SECOND STAGE: MULTIPLE CAUSATION ####
# Covariates: season + distance to wetland + surrounding agriculture + local veg cover

### Informed Null ####
invert.cs$Season <- relevel(invert.cs$Season, ref = "Spring")
m1 <- lm(LogBiomass ~ Season, data = invert.cs)

### Single Models ####
m2 <- lm(LogBiomass ~ PercentAg + Season, data = invert.cs)
m3 <- lm(LogBiomass ~ PercentLocalVeg_50m + Season, data = invert.cs)
m4 <- lm(LogBiomass ~ Dist_Closest_Wetland_m + Season, data = invert.cs)

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

summary(m1)
confint(m1)
plot(m1)

# MODEL AVERAGE ####
model_avg <- model.avg(models)
summary(model_avg)
confint(model_avg)

# Plot model averaged results ####
ma <- summary(model_avg)
df <- as.data.frame(ma$coefmat.full)
ci <- as.data.frame(confint(model_avg, full = T))
df$low <-ci$`2.5 %`
df$upp <-ci$`97.5 %`
setDT(df, keep.rownames = "coefficient")
names(df) <- gsub(" ", "", names(df))
df

# add spring data
spring <- data.table(
  "coefficient" = "SeasonSpring",  
  "Estimate" = -2.1597627,         
  "Std.Error" = 0.3918080,       
  "AdjustedSE" = 0.3983169,        
  "zvalue" = 5.4222218,         
  "Pr(>|z|)" = 0.0000001,             
  "low" = -2.9404495,        
  "upp" = -1.3790759               
)

# Add the "SeasonSpring" row to the original df
df <- rbind(df, spring)

df <- df[c(1, 2, 6, 3, 4, 5), ]

ggplot(data = df[2:6,], aes(x = coefficient, y = Estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", lwd = 1.5) + 
  geom_errorbar(aes(ymin = Estimate - AdjustedSE, ymax = Estimate + AdjustedSE), 
                color="black",
                width = 0, lwd = 1.5) +
  coord_flip()+ # flipping x and y axes
  geom_point(size=5) +
  theme_classic() +
  labs(y = "Model-Averaged Parameter Estimate", 
       x = NULL) +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = c(0.9,0.9)) +
  scale_x_discrete(labels = c("SeasonFall" = "Season (Fall)",
                              "SeasonSpring" = "Season (Spring)",
                              "PercentLocalVeg_50m" = "Local Vegetation Cover",
                              "PercentAg" = "% Surrounding Agriculture",
                              "Dist_Closest_Wetland_m" = "Dist. to Closest Wetland"))


## DISTANCE TO WETLAND & AGRICULTURE MODEL ####
m <- lm(LogBiomass ~ Dist_Closest_Wetland_m + PercentAg + Season, data = invert)

d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                              max(invert$PercentAg), 
                                              length.out = 1000),  
                 Dist_Closest_Wetland_m = mean(invert$Dist_Closest_Wetland_m),
                 Season = unique(invert$Season))

# d$Dist_Closest_Wetland_m <- scale(d$Dist_Closest_Wetland_m, 
#                                 center = mean(invert$Dist_Closest_Wetland_m), 
#                                 scale = sd(invert$Dist_Closest_Wetland_m))
# 
# d$PercentAg <- scale(d$PercentAg, 
#                                   center = mean(invert.cs$PercentAg), 
#                                   scale = sd(invert.cs$PercentAg))
 
predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentAg, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "ln(Macroinvertebrate Biomass)") +
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
   geom_point(data = invert, aes(x = PercentAg, y = LogBiomass,
                                 col = Season), size = 2)

## DISTANCE TO WETLAND & VEG MODEL ####
m <- lm(Biomass ~ PercentLocalVeg_50m + Dist_Closest_Wetland_m + Season, data = invert)

d <- expand.grid(PercentLocalVeg_50m = seq(min(invert$PercentLocalVeg_50m), 
                                 max(invert$PercentLocalVeg_50m), 
                                 length.out = 1000),  
                 Dist_Closest_Wetland_m = mean(invert$Dist_Closest_Wetland_m),
                 Season = unique(invert$Season))

# d$Dist_Closest_Wetland_m <- scale(d$Dist_Closest_Wetland_m, 
#                                 center = mean(invert$Dist_Closest_Wetland_m), 
#                                 scale = sd(invert$Dist_Closest_Wetland_m))
# 
# d$PercentAg <- scale(d$PercentAg, 
#                                   center = mean(invert.cs$PercentAg), 
#                                   scale = sd(invert.cs$PercentAg))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentLocalVeg_50m, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Local Vegetation Cover within 50 m", 
       y = "ln(Macroinvertebrate Biomass)") +
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
  geom_point(data = invert, aes(x = PercentLocalVeg_50m, y = LogBiomass,
                                col = Season), size = 2)

### EFFECT OF WETLAND PERMANENCE? ####
# Compare to the informed null ####
m1 <- lm(LogBiomass ~ Permanence, data = invert.cs)
m2 <- lm(LogBiomass ~ Season, data = invert.cs)

# AIC MODEL SELECTION
model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m1)
confint(m1)
par(mfrow = c(2,2))
