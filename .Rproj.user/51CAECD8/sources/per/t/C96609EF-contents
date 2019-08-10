# Developer: Brady Lange
# Date: 03/30/2019
# Description: MUDAC 2019 - MN water quality

# Set-up workspace
graphics.off()
rm(list = ls())
setwd("C:/Users/brady/Downloads/MUDAC 2019/Problem 1")

# Load libraries
library(tidyverse)
library(readxl)
library(lmtest)
library(DAAG)
library(ISLR)
library(MASS)
library(fmsb)
library(leaps)

# Problem 1
# ----------------------------------------------------------------------------
# Load data
wtr_qlty <- read_excel(path = "MUDAC_data_Problem1.xlsx", sheet = "data", 
                       skip = 2, n_max = 50)
d <- wtr_qlty %>%
    dplyr::select(-c(`Watershed/Monitoring Site Location`, Lat., Long.))

# Explore data
head(d)
tail(d)
names(d)

# Correlation Heatmap
plot(d)
cor(d)
cor_mat <- round(cor(d), 2)
head(cor_mat)
library(reshape2)
melted_cor_mat <- melt(cor_mat)
head(melted_cor_mat)
# Sets lower triangle to NA's
get.upper.tri <- function(cor_mat)
{
    cor_mat[lower.tri(cor_mat)] <- NA
    return(cor_mat)
}
upper_tri <- get.upper.tri(cor_mat)
melted_cor_mat <- melt(upper_tri, na.rm = T)
ggplot(data = melted_cor_mat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle("Correlation Matrix") +
    coord_fixed()

# Leaps
df <- as.data.frame(d2)
leaps(x = df[ , 1:15], y = df[ , 16])

# ---------------------------------------------------------------------------
# Nitrate
# ---------------------------------------------------------------------------
# Nitrate Model
nit_mod <- lm(Nitrate ~ Forest + `Grass and Hay` + Wetlands + Cropland 
              + `Developed/Urban` + Sand + Silt + Clay + `Organic matter` 
              + `Land slope` + Lakes + `Lake interception`, data = d)
summary(nit_mod)
# VIF
VIF(nit_mod)
# AIC
AIC(nit_mod)
# BIC
BIC(nit_mod)
# BP
bptest(nit_mod)
# PRESS
press(nit_mod)

# Shapiro tests
shapiro.test(d2$Forest)
shapiro.test(d2$`Grass and Hay`)
shapiro.test(d2$Wetlands)
shapiro.test(d2$Cropland)
shapiro.test(d2$`Developed/Urban`)
shapiro.test(d2$Sand)  # Good
shapiro.test(d2$`Organic matter`)
shapiro.test(d2$Lakes)

best_d <- d %>%
    dplyr::select(Nitrate, `Grass and Hay`, Wetlands, 
                  `Developed/Urban`, Sand, Lakes)
plot(best_d)
title("Not Transformed", line = 3)

best_d <- d %>%
    dplyr::select(Nitrate, `Grass and Hay`, Wetlands,
                  `Developed/Urban`, Sand, Lakes) %>%
    mutate(Nitrate = I(log(Nitrate)))
plot(best_d)
title("Transformed Data", line = 3)

# Step AIC - find the best model
stepAIC(nit_mod, direction = "both")

# Best model
mod <- lm(Nitrate ~ Forest + `Grass and Hay` + Wetlands + 
                    Cropland + `Developed/Urban` + Sand + `Organic matter` + 
                    Lakes, data = d)
summary(mod)

# Transform model
trans_mod <- lm(I(log(Nitrate)) ~ `Grass and Hay` + Wetlands 
                + `Developed/Urban` + Sand
                + Lakes, data = d)
summary(trans_mod)
bptest(trans_mod) 
VIF(trans_mod)
# AIC
AIC(trans_mod)
# BIC
BIC(trans_mod)
# PRESS
press(trans_mod)

# Correlation Heatmap
cor_mat <- round(cor(best_d), 2)
melted_cor_mat <- melt(cor_mat)
upper_tri <- get.upper.tri(cor_mat)
melted_cor_mat <- melt(upper_tri, na.rm = T)
ggplot(data = melted_cor_mat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle("Nitrate - Correlation Matrix") +
    coord_fixed()

melt_d <- melt(best_d, id.vars = "Nitrate")
ggplot(melt_d) + 
    geom_point(aes(value, Nitrate, color = variable)) + 
    geom_smooth(aes(value, Nitrate, color = variable), method = "lm") +
    facet_wrap(~variable, scales = "free_x") 

# ---------------------------------------------------------------------------
# TSS
# ---------------------------------------------------------------------------
# TSS Model
tss_mod <- lm(`Total Suspended Solid (TSS)` ~ Forest + `Grass and Hay` 
              + Wetlands + Cropland + `Developed/Urban` + Sand + Silt 
              + Clay + `Organic matter` + `Land slope` + Lakes 
              + `Lake interception`, data = d)
summary(tss_mod)
# VIF
VIF(tss_mod)
# AIC
AIC(tss_mod)
# BIC
BIC(tss_mod)
# BP test
bptest(tss_mod)
# PRESS
press(tss_mod)

# Shapiro tests
shapiro.test(d2$Forest)
shapiro.test(d2$`Grass and Hay`)
shapiro.test(d2$Wetlands)
shapiro.test(d2$Cropland)
shapiro.test(d2$`Developed/Urban`)
shapiro.test(d2$Sand) 
shapiro.test(d2$Silt) 
shapiro.test(d2$Clay)
shapiro.test(d2$Lakes)
shapiro.test(d2$`Lake interception`)

best_d <- d %>%
    dplyr::select(`Total Suspended Solid (TSS)`, Forest, `Grass and Hay`, 
                  Wetlands, Cropland, `Developed/Urban`, Sand, Silt, Clay,
                  Lakes, `Lake interception`)
plot(best_d)

best_d <- d %>%
    dplyr::select(`Total Suspended Solid (TSS)`, Forest, `Grass and Hay`, 
                  Wetlands, Cropland, `Developed/Urban`, Sand, Silt, Clay, 
                  Lakes, `Lake interception`) %>%
    mutate(`Total Suspended Solid (TSS)` = I(log(`Total Suspended Solid (TSS)`)))
plot(best_d)

# Step AIC - find the best model
stepAIC(tss_mod, direction = "both")

# Model selection
mod <- lm(`Total Suspended Solid (TSS)` ~ Forest + `Grass and Hay` 
          + Wetlands + Cropland + `Developed/Urban` + Sand + Silt 
          + Clay + Lakes + `Lake interception`, data = d)
summary(mod)

# Transform model
trans_mod <- lm(I(log(`Total Suspended Solid (TSS)`)) ~ Forest + `Grass and Hay` 
                + Wetlands + Cropland + `Developed/Urban` + Sand + Silt 
                + Clay + Lakes + `Lake interception`, data = d)
summary(trans_mod)
bptest(trans_mod) 
VIF(trans_mod)
# AIC
AIC(trans_mod)
# BIC
BIC(trans_mod)
# PRESS
press(trans_mod)

# Correlation Heatmap
cor_mat <- round(cor(best_d), 2)
melted_cor_mat <- melt(cor_mat)
upper_tri <- get.upper.tri(cor_mat)
melted_cor_mat <- melt(upper_tri, na.rm = T)
ggplot(data = melted_cor_mat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle("Total Suspended Solid - Correlation Matrix") +
    coord_fixed()

melt_d <- melt(best_d, id.vars = "Total Suspended Solid (TSS)")
ggplot(melt_d) + 
    geom_point(aes(value, `Total Suspended Solid (TSS)`, color = variable)) + 
    geom_smooth(aes(value, `Total Suspended Solid (TSS)`, color = variable), method = "lm") +
    facet_wrap(~variable, scales = "free_x") 

# ---------------------------------------------------------------------------
# Combined - Nitrate & TSS (combined row averages)
# ---------------------------------------------------------------------------
# Combined Model
# Scale dependent variables
d2 <- d
d2$Nitrate <- scale(d$Nitrate)
d2$`Total Suspended Solid (TSS)` <- scale(d$`Total Suspended Solid (TSS)`)
# Average of dependent variables
d2$Nit_TSS <- rowMeans(d2[ , c(16, 17)])
nit_tss_mod <- lm(Nit_TSS ~ Forest + `Grass and Hay` 
                  + Wetlands + Cropland + `Developed/Urban` + Sand + Silt 
                  + Clay + `Organic matter` + `Land slope` + Lakes 
                  + `Lake interception`, data = d2)
summary(nit_tss_mod)
# VIF 
VIF(nit_tss_mod)
# AIC
AIC(nit_tss_mod)
# BIC
BIC(nit_tss_mod)
# BP test
bptest(nit_tss_mod)
# PRESS
press(nit_tss_mod)

# Shapiro tests
shapiro.test(d2$Forest)
shapiro.test(d2$`Grass and Hay`)
shapiro.test(d2$Wetlands)
shapiro.test(d2$Cropland)
shapiro.test(d2$`Developed/Urban`)
shapiro.test(d2$Sand) 
shapiro.test(d2$Silt)
shapiro.test(d2$Clay)
shapiro.test(d2$`Organic matter`)
shapiro.test(d2$`Land slope`)
shapiro.test(d2$Lakes)
shapiro.test(d2$`Lake interception`)

best_d <- d2 %>%
    dplyr::select(Nit_TSS, Forest, `Grass and Hay`, Wetlands, Cropland, 
                  `Developed/Urban`, Sand, Silt, Clay, `Organic matter`, 
                  `Land slope`, Lakes, `Lake interception`)
plot(best_d)

best_d <- d2 %>%
    dplyr::select(Nit_TSS, Forest, `Grass and Hay`, Wetlands, Cropland, 
                  `Developed/Urban`, Sand, Silt, Clay, `Organic matter`, 
                  `Land slope`, Lakes, `Lake interception`) %>%
    mutate(Cropland = I(1/log(Cropland)), Clay = I(1/log(Clay)), 
           `Developed/Urban` = I(exp(`Developed/Urban`)),
           `Organic matter` = I(1/(`Organic matter`)), Lakes = I(log(Lakes)), 
           Forest = I(log(Forest)), Wetlands = I(log(Wetlands)))
plot(best_d)

# Step AIC - find the best model
stepAIC(nit_tss_mod, direction = "both")

# Model selection
mod <- lm(Nit_TSS ~ Forest + `Grass and Hay` + Wetlands 
          + Cropland + `Developed/Urban` + Sand + Silt + Clay 
          + `Organic matter` + `Land slope` + Lakes 
          + `Lake interception`, data = d2)
summary(mod)

# Transform model
trans_mod <- lm(Nit_TSS ~ I(log(Forest)) + `Grass and Hay` + I(log(Wetlands))
                + I(1/log(Cropland)) + `Developed/Urban` + Sand + Silt 
                + I(1/log(Clay))
                + I(1/(`Organic matter`)) + `Land slope` + I(log(Lakes))
                + `Lake interception`, data = d2)
summary(trans_mod)
VIF(trans_mod)
bptest(trans_mod)
# AIC
AIC(trans_mod)
# BIC
BIC(trans_mod)
# PRESS
press(trans_mod)

# Correlation Heatmap
cor_mat <- round(cor(best_d), 2)
melted_cor_mat <- melt(cor_mat)
upper_tri <- get.upper.tri(cor_mat)
melted_cor_mat <- melt(upper_tri, na.rm = T)
ggplot(data = melted_cor_mat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle("Nitrate/TSS - Correlation Matrix") +
    coord_fixed()