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
invert <- read.csv("processed/Macroinverterbrate_Analysis_2025-03-06.csv")

# Data Distribution
hist(invert$Diversity, breaks = 20) # left  skewed distribution
 
 # Data Exploration ####
 ## Ag & Diversity ####
 ggplot(invert, aes(x = PercentAg, y = Diversity, col = Season)) + 
   geom_point(size = 2.5) +
   theme_classic() + 
   labs(x = "% Surrounding Agriculture",
        y = "Macroinvertebrate Diversity (# Orders)") +
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
                                 "Fall" = "darkolivegreen")) + 
   scale_fill_manual(values = c("Spring" = "skyblue", 
                                "Fall" = "darkolivegreen"))

# Overdispersion
m <- mean(invert$Diversity)
var <- var(invert$Diversity)

ratio <- var/m
print(ratio) # mild overdispersion: 1.59

# Poisson & Negative Binomial Model
m1 <- glm(Diversity ~ PercentAg, data = invert, family = 'poisson')
m2 <- glm.nb(Diversity ~ PercentAg, data = invert)

poisson_dispersion <- sum(residuals(m1, type = "pearson")^2) / m1$df.residual

# Standardized residuals for overdispersion
par(mfcol = c(1,3))
plot(predict(m1), rstandard(m1, type = "pearson"), ylim = c(-5,5), main = "Pearson")
abline(h = c(-2,2), lty = 2)
plot(predict(m1), rstandard(m1, type = "deviance"), ylim = c(-5,5), main = "Deviance")
abline(h = c(-2,2), lty = 2)
plot(predict(m1), rstudent(m1), ylim = c(-5,5), main = "Studentized")
abline(h = c(-2,2), lty = 2)

summary(m)

# overdispersion parameter (residual deviance / df)
op <- 101.59/74
op

1 - pchisq(101.59, df = 74) # significant mild overdispersion

# Residual plot for Poisson regression
m1.res <- resid(m1)
plot(fitted(m1), m1.res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')
abline(0,0)

# Residual plot for negative binomial regression
m2.res <- resid(m2)
plot(fitted(m2), m2.res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Negative Binomial')
abline(0,0)

# Likelihood Ratio Test
pchisq(2 * (logLik(m2) - logLik(m1)), df = 1, lower.tail = FALSE)

# Other residual plots
par(mfrow = c(2,2))
plot(m1)
plot(m2)

# Get the AIC values for both models
aic_m1 <- AIC(m1)
aic_m2 <- AIC(m2)

# Create a data frame to compare AIC values
aic_comparison <- data.frame(
  Model = c("Poisson", "Negative Binomial"),
  AIC = c(aic_m1, aic_m2)
)

# Print the AIC comparison
print(aic_comparison)


# Transform variables

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
invert.cs$SPEI <- scale(invert.cs$SPEI)
invert.cs$PercentAg <- scale(invert.cs$PercentAg)
invert.cs$NearestCropDistance_m <- scale(invert.cs$NearestCropDistance_m)
invert.cs$PercentBufferAroundWetland <- scale(invert.cs$PercentBufferAroundWetland)
invert.cs$MaxBufferWidth_m <- scale(invert.cs$MaxBufferWidth_m)
invert.cs$PercentLocalVeg_50m <- scale(invert.cs$PercentLocalVeg_50m)
invert.cs$pH <- scale(invert.cs$pH)
invert.cs$TDS_mg.L <- scale(invert.cs$TDS_mg.L)
invert.cs$Dist_Closest_Wetland_m <- scale(invert.cs$Dist_Closest_Wetland_m)

# Model Visualization
# Tim's example modified
m <- glm(Diversity ~ Season + PercentAg, data = invert, family = 'poisson')

d <- expand.grid(Season = c("Spring",
                            "Fall"),
                 PercentAg = seq(min(invert$PercentAg),
                                 max(invert$PercentAg),
                                 length = 1000))

d$yhat <- predict(m, newdata = d, type = "response") # type = "response" is required

# Add confidence intervals
library(trtools)
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
  
















