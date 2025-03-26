#------------------------------#
#   Order Diversity Analysis   #
#     Created 2025-03-26       #
#    Modified 2025-03-26       #
#------------------------------#

# load and transform data ----
# load packages
library(ggplot2)
library(dplyr)
library(car)
library(AICcmodavg)
library(MASS)
library(MuMIn)
library(pscl)

# read in data
invert <- read.csv("../data/Macroinverterbrate_Analysis_2025-03-06.csv")

# transform data
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


# standardize continuous variables
invert.cs <- invert
invert.cs$PercentAg <- scale(invert.cs$PercentAg)
invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
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
# max buffer width and percent buffer (0.77)
# permanence and percent buffer (0.668)

# Agriculture
m1 <- glm(Diversity ~ PercentAg + Season, family = "poisson", data = invert.cs)
m2 <- glm(Diversity ~ AgCategory + Season, family = "poisson", data = invert.cs)
m3 <- glm(Diversity ~ NearestCropDistance_m + Season, family = "poisson", 
            data = invert.cs)
m4 <- glm(Diversity ~ DominantCrop + Season, family = "poisson", 
            data = invert.cs)

### ...AIC 
models <- list(m1, m2, m3, m4)
model.sel(models)

# PercentAg is the best predictor

# Vegetation
m1 <- glm(Diversity ~ MaxBufferWidth_m + Season, family = "poisson", 
            data = invert.cs)
m2 <- glm(Diversity ~ PercentBufferAroundWetland + Season, family = "poisson",
            data = invert.cs)

### ...AIC
models <- list(m1, m2)
model.sel(models)

# Maximum buffer width is the best predictor

# ...model selection ----
m1 <- glm(Diversity ~ PercentAg + Season, data = invert.cs, family = "poisson")
m2 <- glm(Diversity ~ MaxBufferWidth_m + Season, data = invert.cs, family = "poisson")
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert.cs, family = "poisson")
m4 <- glm(Diversity ~ MaxBufferWidth_m + Season + PercentLocalVeg_50m, 
          data = invert.cs, family = "poisson")
m5 <- glm(Diversity ~ pH + Season, data = invert.cs, family = "poisson")
m6 <- glm(Diversity ~ TDS_mg.L + Season, data = invert.cs, family = "poisson")
m7 <- glm(Diversity ~ pH + TDS_mg.L + Season, data = invert.cs, family = "poisson")
m8 <- glm(Diversity ~ Dist_Closest_Wetland_m + Season, data = invert.cs, family = "poisson")
m9 <- glm(Diversity ~ Season, data = invert.cs, family = "poisson")

models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
model.sel(models)

deviance(m1)
deviance(m2)
deviance(m3)
deviance(m4)
deviance(m5)
deviance(m6)
deviance(m7)
deviance(m8)
deviance(m9)

# plot parameter estimates ----
# model averaged results (95% confidence set)
models.95 <- list(m1, m2, m3, m4, m5, m6, m8, m9)
model_avg <- model.avg(models.95)
summary(model_avg)
confint(model_avg, method = "conditional")

ma <- summary(model_avg)
df <- as.data.frame(ma$coefmat.full)
ci <- as.data.frame(confint(model_avg, full = F))
df$low <-ci$`2.5 %`
df$upp <-ci$`97.5 %`
setDT(df, keep.rownames = "coefficient")
names(df) <- gsub(" ", "", names(df))
df

ggplot(data = df, aes(x = coefficient, y = Estimate)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", lwd = 1.5) + 
  geom_errorbar(aes(ymin = low, ymax = upp), 
                color="black", width = 0, lwd = 1) +
  coord_flip() +  # flipping x and y axes
  geom_point(size=4) +
  theme_classic() +
  labs(y = "Model-Averaged Parameter Estimate (conditional)", 
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
  scale_x_discrete(limits = c("pH", "TDS_mg.L", "MaxBufferWidth_m", 
                              "PercentLocalVeg_50m", "PercentAg", "Dist_Closest_Wetland_m", 
                              "SeasonFall", "(Intercept)"),
                   labels = c("pH" = "pH",
                              "TDS_mg.L" = "Total Dissolved Solids",
                              "PercentBufferAroundWetland" = "Buffer Cover",
                              "MaxBufferWidth_m" = "Max Buffer Width",
                              "PercentLocalVeg_50m" = "Local Vegetation Cover", 
                              "PercentAg" = "% Surrounding Ag", 
                              "Dist_Closest_Wetland_m" = "Dist. to Closest Wetland",
                              "SeasonFall" = "Season (Fall)", 
                              "(Intercept)" = "Season (Spring)"))

# ...top four models ----
ci3 <- trtools::lincon(m3, fcov=vcov)
ci3 <- as.data.frame(ci3)
ci3$coefficient <- rownames(ci3)

ci4 <- trtools::lincon(m4, fcov=vcov)
ci4 <- as.data.frame(ci4)
ci4$coefficient <- rownames(ci4)

ci9 <- trtools::lincon(m9, fcov=vcov)
ci9 <- as.data.frame(ci9)
ci9$coefficient <- rownames(ci9)

ci1 <- trtools::lincon(m1, fcov=vcov)
ci1 <- as.data.frame(ci1)
ci1$coefficient <- rownames(ci1)

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
  scale_x_discrete(limits = c("PercentLocalVeg_50m","MaxBufferWidth_m", "SeasonFall","(Intercept)"),
                   labels = c("PercentLocalVeg_50m" = "Local Vegetation Cover", 
                              "MaxBufferWidth_m" = "Max Buffer Width",
                              "SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

p9 <- ggplot(data = ci9, aes(x = coefficient, y = estimate)) +
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
  scale_x_discrete(limits = c("SeasonFall","(Intercept)"),
                   labels = c("SeasonFall" = "Season (Fall)",
                              "(Intercept)" = "Season (Spring)")) +
  scale_y_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3))

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

p3 + p4 + p9 + p1

# plot top model ----
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert, 
          family = "poisson")

d <- expand.grid(PercentLocalVeg_50m = seq(min(invert$PercentLocalVeg_50m), 
                                 max(invert$PercentLocalVeg_50m), 
                                 length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m3, newdata = d, se.fit = TRUE, type = "response")

d$predicted_Mass <- predictions$fit

d$lower_CI <- d$predicted_Mass - 1.96 * predictions$se.fit
d$upper_CI <- d$predicted_Mass + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentLocalVeg_50m, y = predicted_Mass, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Local Vegetation Cover", 
       y = "Macroinvertebrate Diversity (# Orders/Sample)") +
  theme(axis.title.x = element_text(size = 19,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 19,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) + 
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = PercentLocalVeg_50m, y = Diversity,
                                col = Season), size = 2)

# plot ag model ----
m3 <- glm(Diversity ~ PercentAg + Season, data = invert, 
          family = "poisson")

d <- expand.grid(PercentAg = seq(min(invert$PercentAg), 
                                           max(invert$PercentAg), 
                                           length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m3, newdata = d, se.fit = TRUE, type = "link")

d$yhat <- predictions$fit

d$lower_CI <- d$yhat - 1.96 * predictions$se.fit
d$upper_CI <- d$yhat + 1.96 * predictions$se.fit

d$yhat <- exp(d$yhat)
d$lower_CI <- exp(d$lower_CI)
d$upper_CI <- exp(d$upper_CI)

ggplot(d, aes(x = PercentAg, y = yhat, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Macroinvertebrate Diversity (# Orders/Sample)") +
  theme(axis.title.x = element_text(size = 19,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 19,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) + 
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = PercentAg, y = Diversity,
                                col = Season), size = 2)

# assumptions ----
m <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert, 
          family = "poisson")
par(mfrow = c(1,1))

# check residuals
plot(predict(m),rstudent(m)) # residuals do not show clear patterns

# Pearson residuals
# slight funnel shape but mostly uniform
pearson_residuals <- residuals(m, type = "pearson")
plot(pearson_residuals, main = "Pearson Residuals", ylab = "Residuals", xlab = "Index")

# Deviance residuals
# slight funnel shape but mostly uniform
deviance_residuals <- residuals(m, type = "deviance")
plot(deviance_residuals, main = "Deviance Residuals", ylab = "Residuals", xlab = "Index")

# Chi-square test for goodness of fit ----
m <- glm(Diversity ~ PercentLocalVeg_50m + Season, data = invert, 
          family = "quasipoisson")
# Calculate the expected values for the biomass response
expected_values <- predict(m, type = "response")

# Calculate the Pearson residuals
pearson_residuals <- (invert$Diversity - expected_values) / sqrt(expected_values)

# Calculate the Pearson chi-square statistic
pearson_chi_square <- sum(pearson_residuals^2)

df <- 77-3
p_value <- 1 - pchisq(92.9, df)
p_value # 0.0679
# Model not a great fit but also not horrible

# negative binomial? ----
m <- glm.nb(Diversity ~ PercentLocalVeg_50m + Season, data = invert)
summary(m) # theta extremely large
d <- expand.grid(PercentLocalVeg_50m = seq(min(invert$PercentLocalVeg_50m), 
                                 max(invert$PercentLocalVeg_50m), 
                                 length.out = 1000),  
                 Season = unique(invert$Season))

predictions <- predict(m, newdata = d, se.fit = TRUE, type = "link")

d$yhat <- predictions$fit

d$lower_CI <- d$yhat - 1.96 * predictions$se.fit
d$upper_CI <- d$yhat + 1.96 * predictions$se.fit

d$yhat <- exp(d$yhat)
d$lower_CI <- exp(d$lower_CI)
d$upper_CI <- exp(d$upper_CI)

ggplot(d, aes(x = PercentLocalVeg_50m, y = yhat, col = Season)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Season), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Local Vegetation Cover", 
       y = "Macroinvertebrate Diversity (# Orders/Sample)") +
  theme(axis.title.x = element_text(size = 19,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 19,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  scale_color_manual(values = c("Spring" = "skyblue3", 
                                "Fall" = "darkorange3")) + 
  scale_fill_manual(values = c("Spring" = "skyblue3", 
                               "Fall" = "darkorange3")) +
  geom_point(data = invert, aes(x = PercentLocalVeg_50m, y = Diversity,
                                col = Season), size = 2)

m <- zeroinfl(Diversity ~ PercentLocalVeg_50m + Season | PercentLocalVeg_50m +
                Season, data = invert, dist = "negbin")
summary(m)
# Season:Fall has a very large SE...model not fitting appropriately
# Resort to Poisson

# compare models to wetland permanence
# Wetland permanence ----
# Compare to other models ####
m1 <- glm(Diversity ~ Permanence, family = "poisson", data = invert.cs)
m2 <- glm(Diversity ~ Season, family = "poisson", data = invert.cs)
m3 <- glm(Diversity ~ PercentLocalVeg_50m + Season, family = "poisson", 
          data= invert.cs)
m4 <- glm(Diversity ~ PercentLocalVeg_50m + Permanence, family = "poisson", 
          data = invert.cs)
m5 <- glm(Diversity ~ 1, family = "poisson", data = invert.cs)

# AIC MODEL SELECTION
models <- list(m1, m2, m3, m4, m5)
model.sel(models)

deviance(m1)
deviance(m2)
deviance(m3)
deviance(m4)






