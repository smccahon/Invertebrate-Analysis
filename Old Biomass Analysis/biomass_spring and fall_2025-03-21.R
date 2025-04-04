#------------------------------#
#   Biomass (S and F) Analysis #
#     Created 2025-03-21       #
#    Modified 2025-03-27       #
#------------------------------#

# Load packages
library(ggplot2)
library(dplyr)
library(car)
library(AICcmodavg)
library(glmmTMB) # my data set has 14% zeros
library(MuMIn)
library(ggbreak) 
library(patchwork)

# Read in data
invert <- read.csv("Macroinverterbrate_Analysis_2025-03-06.csv")

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

invert$BareAg <- factor(invert$BareAg, levels = c("Y", "N"),
                        labels = c("Bare Agriculture",
                                   "Vegetated"))

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



### Raw Data for Biomass Coded by Season ----

ggplot(invert, aes(x = AgCategory, y = Biomass, fill = Season)) + 
  geom_boxplot(alpha = 0.5) +
  theme_classic() +
  labs(x = "% Agriculture Surrounding Wetland", 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "right") +
  scale_color_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen")) + 
  scale_fill_manual(values = c("Spring" = "skyblue", 
                               "Fall" = "darkolivegreen")) +
  scale_y_break(c(0.5,1.40), scales = 0.25,
                c(1.5, 3))

## bare ag ----
bare <- subset(invert, BareAg == "Y")
ggplot(bare, aes(x = Biomass)) + geom_boxplot()

### Raw Data for Biomass Coded by BareAg ----

# all
ggplot(invert, aes(x = BareAg, y = Biomass, fill = BareAg, col = BareAg)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 3) +
  theme_classic() +
  labs(x = NULL, 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = c("Bare Agriculture" = "skyblue", 
                               "Vegetated" = "darkolivegreen")) +
  scale_color_manual(values = c("Bare Agriculture" = "skyblue", 
                                "Vegetated" = "darkolivegreen")) +
  scale_y_break(c(0.5,1.40), scales = 0.25,
                c(1.5, 3))



ggplot(invert, aes(x = BareAg, y = Biomass, fill = BareAg, col = BareAg)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 3) +
  theme_classic() +
  labs(x = NULL, 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = c("Bare Agriculture" = "skyblue", 
                               "Vegetated" = "darkolivegreen")) +
  scale_color_manual(values = c("Bare Agriculture" = "skyblue", 
                               "Vegetated" = "darkolivegreen")) +
  ylim(0, 0.04)


# spring
ggplot(spring, aes(x = BareAg, y = Biomass, fill = BareAg, col = BareAg)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 3) +
  theme_classic() +
  labs(x = NULL, 
       y = "(Spring) Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = c("Bare Agriculture" = "skyblue", 
                               "Vegetated" = "darkolivegreen")) +
  scale_color_manual(values = c("Bare Agriculture" = "skyblue", 
                                "Vegetated" = "darkolivegreen"))



# low veg ag compared to high veg ag
lowveg <- subset(invert, PercentLocalVeg_50m < 15)

invert <- invert %>% 
  mutate(lowveg = case_when(
    PercentLocalVeg_50m < 25 & PercentAg >= 75 ~ "Low Vegetation Cover (< 25%)",
    PercentLocalVeg_50m > 25 & PercentAg >= 75 ~ "High Vegetation Cover (> 25%)",
    TRUE ~ NA_character_
  ))


ggplot(invert %>% filter(!is.na(lowveg)), aes(x = lowveg, y = Biomass, 
                                              fill = lowveg, col = lowveg)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 3) +
  theme_classic() +
  labs(x = NULL, 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = c("Low Vegetation Cover (< 25%)" = "skyblue", 
                               "High Vegetation Cover (> 25%)" = "darkolivegreen")) +
  scale_color_manual(values = c("Low Vegetation Cover (< 25%)" = "skyblue", 
                                "High Vegetation Cover (> 25%)" = "darkolivegreen"))

invert <- invert %>% 
  mutate(lowveg = case_when(
    PercentLocalVeg_50m < 25 & PercentAg >= 25 ~ "Low Vegetation Cover (< 25%)",
    PercentLocalVeg_50m > 25 & PercentAg >= 25 ~ "High Vegetation Cover (> 25%)",
    TRUE ~ NA_character_
  ))


ggplot(invert %>% filter(!is.na(lowveg)), aes(x = lowveg, y = Biomass, 
                                              fill = lowveg, col = lowveg)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 3) +
  theme_classic() +
  labs(x = NULL, 
       y = "Macroinvertebrate Biomass (g)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = c("Low Vegetation Cover (< 25%)" = "skyblue", 
                               "High Vegetation Cover (> 25%)" = "darkolivegreen")) +
  scale_color_manual(values = c("Low Vegetation Cover (< 25%)" = "skyblue", 
                                "High Vegetation Cover (> 25%)" = "darkolivegreen"))

wilcox.test(invert$Biomass ~ invert$lowveg, data = na.omit(invert))
