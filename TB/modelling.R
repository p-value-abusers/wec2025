# File Description -------------------------------------------------------------
# Modelling file for WEC2025
# team: p-value abusers
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
options(scipen=0)
Sys.setenv(LANG = "en")

# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(corrplot)

# Data -------------------------------------------------------------------------
df <- read.table("../merged_dfs/df_all.csv", sep=',', head=TRUE)
dim(df)

df$working_hours <- as.integer(df$hour %in% 8:17)

# Mod0 -------------------------------------------------------------------------
mod0 = lm(sqrt(traffic) ~
            factor(station_id) +
            working_hours +
            is_holiday +
            weekedn
            , data = df)
summary(mod0)
par(mfrow = c(2, 2))
plot(mod0)

# Mod1 -------------------------------------------------------------------------
df['hour_weekend_0'] = df['hour'] * (df['weekedn'] == 0)
df['hour_weekend_1'] = df['hour'] * (df['weekedn'] == 1)
df['hour_weekend'] = df['hour'] * df['weekedn']
df['weekend'] = df['weekedn']

mod1 = gam(traffic ~
             # s(signal_ctd_m_100_0) + 
             signal_ctd_m_100_0 +
             pmin(signal_ctd_m_100_0 - 2, 0) +
             s(hour) +
             # automotive_5 + 
             # roads_intensity_1_0_2 + 
             # factor(station_id) +
             # working_hours + 
             # s(automotive_5) +
             # s(roads_intensity_1_0_2) +
             weekedn,
           data = df)
summary(mod1)
par(mfrow = c(1, 2))
plot(mod1)
gam.check(mod1)
sqrt(mean((mod1$residuals)**2))

df['resid'] = mod1$residuals

df_station <- df %>%
  select(date_station, hour, resid) %>%
  pivot_wider(names_from = hour, values_from = resid) %>% 
  select(-date_station)

df_hour <- df %>%
  select(date_hour, station_id, resid) %>%
  pivot_wider(names_from = station_id, values_from = resid) %>% 
  select(-date_hour)

par(mfrow = c(1, 1))
corrplot(cor(df_station))

par(mfrow = c(1, 1))
corrplot(cor(df_hour))

par(mfrow = c(1, 1))
plot(mod1$fitted.values, mod1$y)

gam.check(mod1)

par(mfrow = c(1, 1))
plot(df$signal_km_1, df$signal_rho_95)

# Mod2 -------------------------------------------------------------------------
# looking at specific stations
df2 = df %>% filter(station_id == 2301)
df2 = df %>% filter(station_id == 3332)
df2 = df %>% filter(station_id == 3324)
df2 = df %>% filter(station_id == 2309)
df2 = df %>% filter(station_id == 2307)
dim(df2)
# check out 2301 !!!

mod1 = gam(traffic ~
             s(hour, by = signal_ctd_m_100_0) +
             signal_ctd_m_100_0 +
             s(hour) +
             s(hour_weekend) +
             s(hour_weekend_0) +
             s(hour_weekend_1) +
             automotive_5 +
             roads_intensity_1_0_2 +
             factor(station_id) +
             working_hours +
             s(automotive_5) +
             s(roads_intensity_1_0_2)
           ,data = df2)
summary(mod1)
par(mfrow = c(1, 2))
plot(mod1)
gam.check(mod1)
sqrt(mean((mod1$residuals)**2))

df2['resid'] = mod1$residuals
df2 %>% select(resid) %>% pull() %>% order()

# Mod3 - STATIC -------------------------------------------------------------------
# Modelling of the base level
df_static <- df %>%
  group_by(station_id, date) %>%  
  summarise(
    traffic_95 = quantile(traffic, 0.95, na.rm = TRUE), 
    .groups = "drop"
  )
df_filtered <- df_static %>%
  filter(!station_id %in% c('3310', '3379', '5342', '5368'))
df2 <- df %>%
  left_join(df_filtered, by = c("station_id", "date"))
print(dim(df2))
df2 <- df2 %>%
  filter(!station_id %in% c('3310', '3379', '5342', '5368'))
print(dim(df2))
df2_full <- df2 %>% mutate(y_scaled = traffic / traffic_95)
dim(df2_full)
df2_full = df2_full %>% filter(traffic > 5)

# considering specific stations
df2 = df2_full %>% filter(station_id == 5347)
# df2 = df2_full %>% filter(station_id == 2301)
# df2 = df2_full %>% filter(station_id == 9398)
# df2 = df2_full %>% filter(station_id == 3332)
# df2 = df2_full %>% filter(station_id == 3324)
# df2 = df2_full %>% filter(station_id == 2309)
# df2 = df2_full %>% filter(station_id == 2307)

mod3 = gam(y_scaled ~
           s(hour_weekend_0) +
           s(hour_weekend_1) +
           weekend
           ,data = df2)
summary(mod3)
par(mfrow = c(1, 2))
plot(mod3, residuals = T, rug = T)
par(mfrow = c(2, 2))
gam.check(mod3)
sqrt(mean((mod3$residuals)**2))

df2['resid'] = mod3$residuals

numeric_columns <- df2 %>% select(where(is.numeric))
correlations <- cor(numeric_columns, use = "complete.obs")
cor_with_resid <- correlations[, "resid"]
round(sort(abs(cor_with_resid)), 3)

df2_weekday = df2 %>% filter(weekedn == 0)
mod4 = gam(y_scaled ~
             s(hour)
           ,data = df2_weekday)
summary(mod4)
par(mfrow = c(1, 1))
plot(mod4, residuals = T, rug = T)

x_1 = 5
x_2 = 8
x_3 = 15

# Mod5 - Parametrization of the hourly pattern ---------------------------------
mod5 = lm(y_scaled ~
             hour + I(hour^2) +
            pmax(hour - x_1, 0) +
            pmax(hour - x_2, 0) +
            pmax(hour - x_3, 0) +
            I(pmax((hour - x_1), 0)^2) +
            I(pmax((hour - x_2), 0)^2) +
            I(pmax((hour - x_3), 0)^2)
            # I(pmax((hour - x_2)^2, 0) + pmax((hour - x_3)^2, 0))
           ,data = df2_weekday)
summary(mod5)
par(mfrow = c(2, 2))
plot(mod5)

df2_weekday['y_hat'] = mod5$fitted.values
par(mfrow = c(1, 1))
plot(df2_weekday$hour, df2_weekday$y_scaled, main = "Station_id: 5347")
points(df2_weekday$hour, df2_weekday$y_hat, col = "red", pch = 16)
# lines(df2_weekday$hour, df2_weekday$y_hat, col = "red", lwd = 2)

## Mod7 - Extending the hourly pattern model -----------------------------------
x_1 = 5
x_2 = 8
x_3 = 15

df2_full['health_and_medical_wd1_SQRT'] = sqrt(df2_full['health_and_medical_wd1'])
df2_full_weekday = df2_full %>% filter(weekedn == 0) %>% filter(is_holiday == 0)

mod7 = lm(y_scaled ~
           hour + I(hour^2) +
           pmax(hour - x_1, 0) +
           pmax(hour - x_2, 0) +
           pmax(hour - x_3, 0) +
           I(pmax((hour - x_1), 0)^2) +
           I(pmax((hour - x_2), 0)^2) +
           I(pmax((hour - x_3), 0)^2) +
           
           I((Residential_2 / 10^6) * pmax(hour - x_1, 0)) +
           I((Residential_2 / 10^6) * pmax(hour - x_2, 0)) +
           I((Residential_2 / 10^6) * pmax(hour - x_3, 0)) +
           I((Residential_2 / 10^6) * pmax((hour - x_1), 0)^2) +
           I((Residential_2 / 10^6) * pmax((hour - x_2), 0)^2) +
           I((Residential_2 / 10^6) * pmax((hour - x_3), 0)^2) +
           
           I((health_and_medical_wd1_SQRT) * pmax(hour - x_1, 0)) +
           I((health_and_medical_wd1_SQRT) * pmax(hour - x_2, 0)) +
           I((health_and_medical_wd1_SQRT) * pmax(hour - x_3, 0)) +
           I((health_and_medical_wd1_SQRT) * pmax((hour - x_1), 0)^2) +
           I((health_and_medical_wd1_SQRT) * pmax((hour - x_2), 0)^2) +
           I((health_and_medical_wd1_SQRT) * pmax((hour - x_3), 0)^2)

         ,data = df2_full_weekday)
summary(mod7)
par(mfrow = c(2, 2))
plot(mod7)

df2_full_weekday['y_hat'] = mod7$fitted.values

df2_full_plot_1 = df2_full_weekday %>% filter(station_id == 9398)
df2_full_plot_2 = df2_full_weekday %>% filter(station_id == 5347)
dim(df2_full_plot_1)
dim(df2_full_plot_2)

par(mfrow = c(1, 2))
plot(df2_full_plot_1$hour, df2_full_plot_1$y_scaled, main = "Station_id: 1365")
points(df2_full_plot_1$hour, df2_full_plot_1$y_hat, col = "red", pch = 16)

plot(df2_full_plot_2$hour, df2_full_plot_2$y_scaled, main = "Station_id: 2307")
points(df2_full_plot_2$hour, df2_full_plot_2$y_hat, col = "red", pch = 16)


par(mfrow = c(1, 1))
plot(df2_full$hour, df2_full$y_scaled)
# lines(df2_full$hour, df2_full$y_hat, col = "red")
points(df2_full$hour, df2_full$y_hat, col = "red", pch = 16)