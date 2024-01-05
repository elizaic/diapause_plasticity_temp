# Packages
library(tidyverse)
library(readxl)
library(glmmTMB)
library(emmeans)
library(car)
library(DHARMa)
library(performance)
library(ggpubr)
library(survival)
library(survminer)
library(MASS)
library(viridis)
library(lmtest)

# CDL data ----
oviposition_data <- read_excel("./plasticity_experiment_data/diapause_datasheet_final_data.xlsx", col_types = c("guess","guess","guess","guess","guess", "date","date","date","guess","date", "guess", "guess"))

oviposition_data$population <- factor(oviposition_data$population, levels = c("De", "Sg", "Ci"))
oviposition_data$daylength <- factor(oviposition_data$daylength)
oviposition_data$temperature <- factor(oviposition_data$temperature, levels = c("38/23", "28/13"))
oviposition_data <- oviposition_data %>% mutate(time = end_date - eclosion_date_F) %>%
  mutate(daylength_hr = if_else(daylength == 1030, 10.5,
                 if_else(daylength==1125, 11.42,
                         if_else(daylength==1220, 12.33,
                                 if_else(daylength==1315, 13.25,
                                         if_else(daylength==1410, 14.17,
                                                 if_else(daylength==1505, 15.08,111111111)))))))
str(oviposition_data)

# Data exploration ----

summary_ovi <- oviposition_data %>%
  group_by (daylength,temperature,population) %>%
  summarise(
    n = n(),
    reproductive = sum(eggs_present, na.rm = T),
    not_reproductive = n - reproductive,
    reproductive_prop = reproductive/n,
    daylength_hr = mean(daylength_hr))
summary_ovi$temperature <- factor(summary_ovi$temperature, levels = c("38/23", "28/13"), labels = c("38", "28"))


## By population OLD VERSION
ggplot(summary_ovi, aes(x = daylength, y = reproductive_prop, color = temperature, group = temperature)) +
  geom_point(size = 3) +
  geom_line(lwd = 1.2) +
  facet_wrap( ~ population, nrow = 3) +
  theme_bw()

pop_labels_2020 <- c("De" = "Delta\n(39°N)", "Sg" = "St. George\n(37°N)", "Ci" = "Cibola\n(33°N)")
CDL_est_daylength_hr <- CDL_est %>% rename(daylength_hr = diap_50_est) %>%
  mutate(reproductive_prop = rep(-Inf,6))
CDL_est_daylength_hr$temperature <- factor(CDL_est_daylength_hr$temperature, levels = c("38", "28"))

## By population NEW VERSION
ggplot(summary_ovi, aes(x = daylength_hr, y = reproductive_prop, color = temperature, group = temperature)) +
  geom_vline(data = CDL_est_daylength_hr, aes(xintercept = daylength_hr, color = temperature), size = .8) +
  geom_rect(data = CDL_est_daylength_hr,
            aes(xmin=daylength_hr-SE, xmax=daylength_hr+SE, ymin=reproductive_prop, ymax=Inf, fill = temperature),
            alpha = 0.2, color = NA) +
  geom_point(size = 3, position = position_dodge(width = 0.05)) +
  geom_line(lwd = 1.2,  position = position_dodge(width = 0.05), aes(group = temperature)) +
  facet_wrap( ~ population, nrow = 4, strip.position = "right", labeller = labeller(population = pop_labels_2020)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = 'none') +
  labs(y = "Proportion Reproductive", x = "Daylength (hr)", color = "Temperature\nRegime") +
  theme_bw(base_size = 15) +
  theme(strip.background=element_rect(color = 'black', fill="white", size = 1))



## By photoperiod
ggplot(data = summary_ovi) +
  geom_point(aes(x = temperature, y = reproductive_prop,
                 shape = population, color = population), alpha = .8, size = 4) +
  geom_path(aes(x = temperature, y = reproductive_prop, group = population, color = population), size = 1) +
  theme_bw() +
  facet_grid(rows = vars(daylength))

## Plasticity
cool <- summary_ovi %>%
  filter(temperature == "28/13") %>%
  ungroup() %>%
  select(population, daylength, reproductive_prop) %>%
  rename(cool_reproductive_prop=reproductive_prop)

warm <- summary_ovi %>%
  filter(temperature == "38/23") %>%
  ungroup() %>%
  select(population, daylength, reproductive_prop) %>%
  rename(warm_reproductive_prop=reproductive_prop)

plasticity <- full_join(cool, warm) %>%
  mutate(temp_diff = warm_reproductive_prop - cool_reproductive_prop)

ggplot(plasticity) +
  geom_col(aes(x = daylength, y = temp_diff,
               fill = population), position = "dodge", alpha = .8, size = 4) +
  ylab("proportion reproproductive in warm - cool treamtements") +
  theme_bw()

ggplot(plasticity) +
  geom_hline(yintercept =0) +
  geom_point(aes(x = daylength, y = temp_diff,
                 color = population, shape = population), alpha = .8, size = 4) +
  geom_path(aes(x = daylength, y = temp_diff, group = population, color = population), size = 1) +
  ylab("proportion reproproductive in warm - cool treamtements") +
  theme_bw()

# Survival Analysis ----

SurvObj <- Surv(oviposition_data$time, oviposition_data$eggs_present)

model_surv1 <- survfit(SurvObj ~ population + temperature + daylength, data = oviposition_data)
model_surv1
summary(model_surv1)
print(survdiff(SurvObj ~ population, data = oviposition_data), digits = 6)

ggsurvplot(model_surv1, data = oviposition_data,
           surv.median.line = "v",
           conf.int = T)


# Formal Modeling ----
str(oviposition_data)
model1 <- glmmTMB(eggs_present ~ daylength_hr * temperature * population, data = oviposition_data,
                  family = binomial)

check_collinearity(model1)                         # check for collinearity (with no interactions in model)
simulateResiduals(model1) %>% plot()               # QQ plot & resid vs predicted

summary(model1)                                    # model results
Anova(model1, type = 3)
lrtest(model1, model1a_2020)
emmeans(model1, pairwise ~ temperature|daylength_hr|population, type = "response")

#likelihood ratio tests



model1a_2020 <- glm(eggs_present ~ daylength_hr * temperature * population, data = oviposition_data,
                    family = binomial)
simulateResiduals(model1a_2020) %>% plot()               # QQ plot & resid vs predicted
summary(model1a_2020)
Anova(model1a_2020, type = 3)

# CDL ----
#DE cool
DeCool_glm <- glm(eggs_present ~ daylength_hr,
                 data = oviposition_data %>% filter(population == "De" & temperature == "28/13"),
                 family = binomial)

#DE warm
DeWarm_glm <- glm(eggs_present ~ daylength_hr,
                  data = oviposition_data %>% filter(population == "De" & temperature == "38/23"),
                  family = binomial)
#SG cool
SgCool_glm <- glm(eggs_present ~ daylength_hr,
                  data = oviposition_data %>% filter(population == "Sg" & temperature == "28/13"),
                  family = binomial)
#SG warm
SgWarm_glm <- glm(eggs_present ~ daylength_hr,
                  data = oviposition_data %>% filter(population == "Sg" & temperature == "38/23"),
                  family = binomial)
#CI cool
CiCool_glm <- glm(eggs_present ~ daylength_hr,
                  data = oviposition_data %>% filter(population == "Ci" & temperature == "28/13"),
                  family = binomial)
#CI warm
CiWarm_glm <- glm(eggs_present ~ daylength_hr,
                  data = oviposition_data %>% filter(population == "Ci" & temperature == "38/23"),
                  family = binomial)
coef(DeCool_glm)

# Calculate CDLs, SE, and CIs
DeCool_cdl <- dose.p(DeCool_glm, cf = 1:2, p= 0.5)
DeCool_cdl_CI <- DeCool_cdl + attr(DeCool_cdl, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

DeWarm_cdl <- dose.p(DeWarm_glm, cf = c(1,2), p= 0.5)
DeWarm_cdl_CI <- DeWarm_cdl + attr(DeWarm_cdl, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

SgCool_cdl <- dose.p(SgCool_glm, cf = c(1,2), p= 0.5)
SgCool_cdl_CI <- SgCool_cdl + attr(SgCool_cdl, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

SgWarm_cdl <- dose.p(SgWarm_glm, cf = c(1,2), p= 0.5)
SgWarm_cdl_CI <- SgWarm_cdl + attr(SgWarm_cdl, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

CiCool_cdl <- dose.p(CiCool_glm, cf = c(1,2), p= 0.5)
CiCool_cdl_CI <- CiCool_cdl + attr(CiCool_cdl, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

CiWarm_cdl <- dose.p(CiWarm_glm, cf = c(1,2), p= 0.5)
CiWarm_cdl_CI <- CiWarm_cdl + attr(CiWarm_cdl, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

#Other way to calculate CIs
#Wald (same caclulation as above)
DeCool_V <- vcov(DeCool_glm)
DeCool_cdl2 <- -DeCool_glm$coefficients[1]/DeCool_glm$coefficients[2]

DeCool_MoE <- qnorm(1-0.05/2)*sqrt( (DeCool_cdl2^2) *
                        ((DeCool_V[1,1]/DeCool_glm$coefficients[1]^2) +
                        (DeCool_V[2,2]/DeCool_glm$coefficients[2]^2) -
                        (2*DeCool_V[1,2]/(DeCool_glm$coefficients[1]*DeCool_glm$coefficients[2]))))
c(DeCool_cdl2 - DeCool_MoE, DeCool_cdl2 + DeCool_MoE)

#Fieller's Method (definitely worse and not recommended by two papers)
(-((DeCool_V[1,2] * (qnorm(1-0.05/2))^2) + (DeCool_glm$coefficients[1]*DeCool_glm$coefficients[2])) +
  sqrt( (DeCool_V[1,2] * (qnorm(1-0.05/2))^2 + (DeCool_glm$coefficients[1]*DeCool_glm$coefficients[2]))^2 -
        (DeCool_glm$coefficients[2]^2 - (DeCool_V[2,2] * (qnorm(1-0.05/2))^2)) *
        (DeCool_glm$coefficients[1]^2 - (DeCool_V[1,1] * (qnorm(1-0.05/2))^2)) ) ) /
  DeCool_glm$coefficients[2]^2 - (DeCool_V[2,2] * (qnorm(1-0.05/2))^2)

(-((DeCool_V[1,2] * (qnorm(1-0.05/2))^2) + (DeCool_glm$coefficients[1]*DeCool_glm$coefficients[2])) -
    sqrt( (DeCool_V[1,2] * (qnorm(1-0.05/2))^2 + (DeCool_glm$coefficients[1]*DeCool_glm$coefficients[2]))^2 -
            (DeCool_glm$coefficients[2]^2 - (DeCool_V[2,2] * (qnorm(1-0.05/2))^2)) *
            (DeCool_glm$coefficients[1]^2 - (DeCool_V[1,1] * (qnorm(1-0.05/2))^2)) ) ) /
  DeCool_glm$coefficients[2]^2 - (DeCool_V[2,2] * (qnorm(1-0.05/2))^2)



# CDL summary
CDL_est <- rbind(
  cbind(DeCool_cdl, SE = attributes(DeCool_cdl)$SE[,1], lower.CI = DeCool_cdl_CI[,1], upper.CI = DeCool_cdl_CI[,2], population = "De", temperature = "28"),
  cbind(DeWarm_cdl, SE = attributes(DeWarm_cdl)$SE[,1], lower.CI = DeWarm_cdl_CI[,1], upper.CI = DeWarm_cdl_CI[,2], population = "De", temperature = "38"),
  cbind(CiCool_cdl, SE = attributes(CiCool_cdl)$SE[,1], lower.CI = CiCool_cdl_CI[,1], upper.CI = CiCool_cdl_CI[,2], population = "Ci", temperature = "28"),
  cbind(CiWarm_cdl, SE = attributes(CiWarm_cdl)$SE[,1], lower.CI = CiWarm_cdl_CI[,1], upper.CI = CiWarm_cdl_CI[,2], population = "Ci", temperature = "38"),
  cbind(SgWarm_cdl, SE = attributes(SgWarm_cdl)$SE[,1], lower.CI = SgWarm_cdl_CI[,1], upper.CI = SgWarm_cdl_CI[,2], population = "Sg", temperature = "38"),
  cbind(SgCool_cdl, SE = attributes(SgCool_cdl)$SE[,1], lower.CI = SgCool_cdl_CI[,1], upper.CI = SgCool_cdl_CI[,2], population = "Sg", temperature = "28"))
CDL_est <- as.data.frame(CDL_est)
rownames(CDL_est) <- NULL
colnames(CDL_est)[1] <- "diap_50_est"
CDL_est$diap_50_est <- as.numeric(as.character(CDL_est$diap_50_est))
CDL_est$SE <- as.numeric(as.character(CDL_est$SE))
CDL_est$lower.CI <- as.numeric(as.character(CDL_est$lower.CI))
CDL_est$upper.CI <- as.numeric(as.character(CDL_est$upper.CI))
CDL_est$population = factor(CDL_est$population, levels=c("De", "Sg", "Ci"))
CDL_est

# CDL plasticity plot
ggplot(CDL_est, aes(x = temperature, y = diap_50_est, color = population)) +
  geom_point(aes(shape = population), size = 3) +
  geom_line(aes(group = population)) +
  geom_errorbar(aes(x = temperature, ymin = lower.CI, ymax = upper.CI), width = .2) +
  # scale_color_manual(values = c("#9e9ac8", "#756bb1", "#54278f")) +
  ylab("Estimated CDL") +
  theme_bw()

# Formal modeling of CDL ----
lm(diap_50_est ~ temperature, data = CDL_est %>% filter(population == 'De'))
lm(diap_50_est ~ temperature, data = CDL_est %>% filter(population == 'Sg'))
lm(diap_50_est ~ temperature, data = CDL_est %>% filter(population == 'Ci'))

cdl2020_mod <- lm(diap_50_est ~ population * temperature, data = CDL_est)
summary(cdl2020_mod)
Anova(cdl2020_mod, type = 3)

# 2022 Experimental Results ----
# data

oviposition_data2022 <- read_excel("./plasticity_experiment_data/2022_temp_photoperiod_data.xlsx")

oviposition_data2022$population <- factor(oviposition_data2022$population, levels = c("De", "Sg", "Bi", "Ci"))
oviposition_data2022$daylength <- factor(oviposition_data2022$daylength)
oviposition_data2022$daylength_block <- factor(oviposition_data2022$daylength_block)
oviposition_data2022$temperature <- factor(oviposition_data2022$temperature, levels = c("38", "28"))
str(oviposition_data2022)

pop_labels <- c("De" = "Delta\n(39°N)", "Sg" = "St. George\n(37°N)", "Bi" = "Big Bend\n(35°N)", "Ci" = "Cibola\n(33°N)")
temp_labels = c("38" = "38°", "28" = "28°")

#Data exploration/visualization
## By population OLD VERSION
ggplot(oviposition_data2022, aes(x = daylength, y = prop_repro, color = temperature, group = temperature)) +
  geom_point(size = 3) +
  geom_line(lwd = 1.2) +
  facet_wrap( ~ population, nrow = 4) +
  theme_bw()

#wrangle estimated CDL data so it fits into the plot below
CDL_est_22_daylength_hr <- CDL_est_22 %>% rename(daylength_hr = diap_50_est) %>%
  mutate(prop_repro = rep(-Inf,8))
CDL_est_22_daylength_hr$temperature <- factor(CDL_est_22_daylength_hr$temperature, levels = c("38", "28"))


## By population NEW VERSION
ggplot(oviposition_data2022, aes(x = daylength_hr, y = prop_repro, color = temperature, group = temperature)) +
  geom_vline(data = CDL_est_22_daylength_hr, aes(xintercept = daylength_hr, color = temperature), size = .8) +
  geom_rect(data = CDL_est_22_daylength_hr,
            aes(xmin=daylength_hr-SE, xmax=daylength_hr+SE, ymin=prop_repro, ymax=Inf, fill = temperature),
            alpha = 0.2, color = NA) +
  geom_point(size = 3, position = position_dodge(width = 0.05)) +
  geom_line(lwd = 1.2,  position = position_dodge(width = 0.05), aes(group = interaction(temperature, daylength_block))) +
  facet_wrap( ~ population, nrow = 4, strip.position = "right", labeller = labeller(population = pop_labels)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = 'none') +
    labs(y = "Proportion Reproductive", x = "Daylength (hr)", color = "Temperature\nRegime",
       shape = "Experiment Block") +
  theme_bw(base_size = 15) +
  theme(strip.background=element_rect(color = 'black', fill="white", size = 1))

#Formal modeling
model1_2022 <- glmmTMB(cbind(no_reproductive, no_diapause) ~ daylength_hr * temperature * population + (1|daylength_block), data = oviposition_data2022,
                  family = binomial)

check_collinearity(model1_2022)                         # check for collinearity (with no interactions in model)
simulateResiduals(model1_2022) %>% plot()               # QQ plot & resid vs predicted

summary(model1_2022)                                    # model results
Anova(model1_2022, type = 3)
emmeans(model1_2022, pairwise ~ temperature|daylength_hr|population, type = "response")


# CDL 2022 ----
#DE cool
DeCool_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                  data = oviposition_data2022 %>% filter(population == "De" & temperature == "28"),
                  family = binomial(link = 'logit'))
#DE warm
DeWarm_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                  data = oviposition_data2022 %>% filter(population == "De" & temperature == "38"),
                  family = binomial)
#SG cool
SgCool_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                  data = oviposition_data2022 %>% filter(population == "Sg" & temperature == "28"),
                  family = binomial)
#SG warm
SgWarm_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                  data = oviposition_data2022 %>% filter(population == "Sg" & temperature == "38"),
                  family = binomial)
#BI cool
BiCool_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                     data = oviposition_data2022 %>% filter(population == "Bi" & temperature == "28"),
                     family = binomial)
#BI warm
BiWarm_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                     data = oviposition_data2022 %>% filter(population == "Bi" & temperature == "38"),
                     family = binomial)
#CI cool
CiCool_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                  data = oviposition_data2022 %>% filter(population == "Ci" & temperature == "28"),
                  family = binomial)
#CI warm
CiWarm_glm_22 <- glm(cbind(no_reproductive, no_diapause) ~ daylength_hr,
                  data = oviposition_data2022 %>% filter(population == "Ci" & temperature == "38"),
                  family = binomial)
coef(DeCool_glm)

# # Calculate CDLs
DeCool_cdl_22 <- dose.p(DeCool_glm_22, cf = 1:2, p= 0.5)
DeCool_cdl_22_CI <- DeCool_cdl_22 + attr(DeCool_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

DeWarm_cdl_22 <- dose.p(DeWarm_glm_22, cf = c(1,2), p= 0.5)
DeWarm_cdl_22_CI <- DeWarm_cdl_22 + attr(DeWarm_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

SgCool_cdl_22 <- dose.p(SgCool_glm_22, cf = c(1,2), p= 0.5)
SgCool_cdl_22_CI <- SgCool_cdl_22 + attr(SgCool_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

SgWarm_cdl_22 <- dose.p(SgWarm_glm_22, cf = c(1,2), p= 0.5)
SgWarm_cdl_22_CI <- SgWarm_cdl_22 + attr(SgWarm_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

BiCool_cdl_22 <- dose.p(BiCool_glm_22, cf = c(1,2), p= 0.5)
BiCool_cdl_22_CI <- BiCool_cdl_22 + attr(BiCool_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

BiWarm_cdl_22 <- dose.p(BiWarm_glm_22, cf = c(1,2), p= 0.5)
BiWarm_cdl_22_CI <- BiWarm_cdl_22 + attr(BiWarm_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

CiCool_cdl_22 <- dose.p(CiCool_glm_22, cf = c(1,2), p= 0.5)
CiCool_cdl_22_CI <- CiCool_cdl_22 + attr(CiCool_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

CiWarm_cdl_22 <- dose.p(CiWarm_glm_22, cf = c(1,2), p= 0.5)
CiWarm_cdl_22_CI <- CiWarm_cdl_22 + attr(CiWarm_cdl_22, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)

# CDL summary
CDL_est_22 <- rbind(
  cbind(DeCool_cdl_22, SE = attributes(DeCool_cdl_22)$SE[,1], lower.CI = DeCool_cdl_22_CI[,1], upper.CI = DeCool_cdl_22_CI[,2], population = "De", temperature = "28"),
  cbind(DeWarm_cdl_22, SE = attributes(DeWarm_cdl_22)$SE[,1], lower.CI = DeWarm_cdl_22_CI[,1], upper.CI = DeWarm_cdl_22_CI[,2], population = "De", temperature = "38"),
  cbind(CiCool_cdl_22, SE = attributes(CiCool_cdl_22)$SE[,1], lower.CI = CiCool_cdl_22_CI[,1], upper.CI = CiCool_cdl_22_CI[,2], population = "Ci", temperature = "28"),
  cbind(CiWarm_cdl_22, SE = attributes(CiWarm_cdl_22)$SE[,1], lower.CI = CiWarm_cdl_22_CI[,1], upper.CI = CiWarm_cdl_22_CI[,2], population = "Ci", temperature = "38"),
  cbind(BiCool_cdl_22, SE = attributes(BiCool_cdl_22)$SE[,1], lower.CI = BiCool_cdl_22_CI[,1], upper.CI = BiCool_cdl_22_CI[,2], population = "Bi", temperature = "28"),
  cbind(BiWarm_cdl_22, SE = attributes(BiWarm_cdl_22)$SE[,1], lower.CI = BiWarm_cdl_22_CI[,1], upper.CI = BiWarm_cdl_22_CI[,2], population = "Bi", temperature = "38"),
  cbind(SgWarm_cdl_22, SE = attributes(SgWarm_cdl_22)$SE[,1], lower.CI = SgWarm_cdl_22_CI[,1], upper.CI = SgWarm_cdl_22_CI[,2], population = "Sg", temperature = "38"),
  cbind(SgCool_cdl_22, SE = attributes(SgCool_cdl_22)$SE[,1], lower.CI = SgCool_cdl_22_CI[,1], upper.CI = SgCool_cdl_22_CI[,2], population = "Sg", temperature = "28"))
CDL_est_22 <- as.data.frame(CDL_est_22)
rownames(CDL_est_22) <- NULL
colnames(CDL_est_22)[1] <- "diap_50_est"
CDL_est_22$diap_50_est <- as.numeric(as.character(CDL_est_22$diap_50_est))
CDL_est_22$SE <- as.numeric(as.character(CDL_est_22$SE))
CDL_est_22$lower.CI <- as.numeric(as.character(CDL_est_22$lower.CI))
CDL_est_22$upper.CI <- as.numeric(as.character(CDL_est_22$upper.CI))
CDL_est_22$population = factor(CDL_est_22$population, levels=c("De", "Sg", "Bi", "Ci"))
CDL_est_22

# CDL plasticity plot
ggplot(CDL_est_22, aes(x = temperature, y = diap_50_est, color = population)) +
  geom_point(aes(shape = population), size = 3) +
  geom_line(aes(group = population)) +
  geom_errorbar(aes(x = temperature, ymin = lower.CI, ymax = upper.CI), width = .2) +
  # scale_color_manual(values = c("#9e9ac8", "#756bb1", "#54278f")) +
  ylab("Estimated CDL") +
  theme_bw()
#
# Combine CDL estimates from both years
CDL_est_22 <- CDL_est_22 %>% mutate(year = rep(2022, 8))
CDL_est <- CDL_est %>% mutate(year = rep(2020, 6))
CDL_bothyears <- rbind(CDL_est, CDL_est_22)
CDL_bothyears$year <- as.factor(CDL_bothyears$year)
CDL_bothyears$population = factor(CDL_bothyears$population, levels=c("De", "Sg", "Bi", "Ci"), labels = c("Delta (39°N)", "St. George (37°N)", "Big Bend (35°N)", "Cibola (33°N)"))


# Plots with both years ----
#Combine data
ovi_bothyears <- bind_rows(oviposition_data2022 %>% rename(n = no_pairs,
                                reproductive = no_reproductive,
                                not_reproductive = no_diapause,
                                reproductive_prop = prop_repro) %>%
                                mutate(year = 2022),
  summary_ovi %>% mutate(year = 2020) %>% mutate(daylength_block = "year1"))
ovi_bothyears$daylength_block <- as.factor(ovi_bothyears$daylength_block)
ovi_bothyears$year <- as.factor(ovi_bothyears$year)


CDL_est_plasticity_bothyears <- bind_rows(CDL_est_22_daylength_hr %>%
                                            rename(reproductive_prop = prop_repro),
                                          CDL_est_daylength_hr)
write.csv(CDL_est_plasticity_bothyears, "CDL.estimates.csv")

##both years in plasticity plot
ggplot(ovi_bothyears, aes(x = daylength_hr, y = reproductive_prop, color = temperature, group = temperature)) +

  #cdl vertical lines and shaded areas
  geom_vline(data = CDL_est_plasticity_bothyears, aes(xintercept = daylength_hr, color = temperature), size = .8, lty = 2) +
  geom_rect(data = CDL_est_plasticity_bothyears,
            aes(xmin=lower.CI, xmax=upper.CI, ymin=reproductive_prop, ymax=Inf, fill = temperature),
            alpha = 0.2, color = NA) +

  #proportion data points and lines
  geom_point(size = 2.5, position = position_dodge(width = 0.05)) +
  # geom_line(lwd = 1.2,  position = position_dodge(width = 0.05), aes(group = temperature)) +

  #predicted line
  # geom_line(aes(x = daylength_hr, y = 1-prediction, color = temperature, group = temperature)) +
  geom_line(data = ovi_bothyears_predictnew, aes(x = daylength_hr, y = predicted_prob, color = temperature, group = temperature),
            lwd = .9) +
  geom_ribbon(data = ovi_bothyears_predictnew, aes(x = daylength_hr, ymin = ci_low, ymax = ci_high,
                                                   fill = temperature, group = temperature),
              inherit.aes = F, alpha = 0.2) +

  # formatting
  facet_grid(population ~ year, labeller = labeller(population = pop_labels)) +
  # scale_color_viridis_d(direction = -1, aesthetics = c('color', 'fill')) +
  scale_color_manual(values = c("#E69F00", "#0072B2"), aesthetics = c('color', 'fill'),
                     labels = c("Warm: 38/23°C", "Cool: 28/13°C"),
                     guide = guide_legend(reverse = T)) +
  scale_y_continuous(n.breaks = 3) +
  guides(fill = 'none') +
  labs(y = "Proportion Reproductive", x = "Daylength (hr)", color = "Temperature Regime") +
  theme_bw(base_size = 15) +
  theme(strip.background=element_rect(color = 'black', fill="white", size = 1),
        legend.position = 'bottom',
        panel.grid = element_blank()
  )
##Export: 654 x 530

# CDL plasticity plot both years

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

ggplot(CDL_bothyears, aes(x = temperature, y = diap_50_est, group = interaction(population, year), color = population, shape = year)) +
  geom_errorbar(aes(x = temperature, ymin = lower.CI, ymax = upper.CI), width = .55, position = position_dodge(0.1)) +
  geom_line(position = position_dodge(0.1), aes(color = population)) +
  geom_point(size = 3.5, position = position_dodge(0.1), aes(color = population)) +
  scale_color_viridis_d() +
  scale_shape_manual(values = c(19, 17)) +
  scale_x_discrete(labels = c('38' = "Warm: 38/23°C", '28' = "Cool: 28/13°C")) +
  labs(y = "Critical Daylength (hrs)", x = "Temperature Regime", color = "Collection Site",
       shape = "Experiment Year", linetype = "Experiment Year") +
  annotate('text', x = 0.45, y = 15.8, label = "Earlier\ndiapause", hjust = 0) +
  annotate('text', x = 0.45, y = 6.5, label = "Later\ndiapause", hjust = 0) +
  # geom_text(x = 0, y = 15.8, hjust = 1, color = 'black', label = "Earlier\ndiapause")+
  # coord_cartesian(clip = 'off') +
  guides(shape = guide_legend(order = 2), linetype = guide_legend(order = 2), color = guide_legend(order = 1)) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())
#
#Export 600 x 425 px


##Formal modeling of CDL both years ----

cdl_mod <- lm(diap_50_est ~ population * temperature, data = CDL_bothyears)
summary(cdl_mod)
Anova(cdl_mod, type = 3)
emmeans(cdl_mod, pairwise ~ temperature|population)

## Daylength at sites
library(geosphere)
day_of_year <- seq(1:365)
delta_daylength <- daylength(lat = 39.10004, doy = day_of_year)
max(delta_daylength)
min(delta_daylength)

stgeorge_daylength <- daylength(lat = 37.086714, doy = day_of_year)

bigbend_daylength <- daylength(lat = 35.060760, doy = day_of_year)
cibola_daylength <- daylength(lat = 33.303369, doy = day_of_year)
cbind(day_of_year,delta_daylength,stgeorge_daylength,bigbend_daylength,cibola_daylength) %>% View()

as.Date(c(219,210,202,211,224,216,230,257,270,279,281), origin = "2022-01-01")

# Formal modeling of oviposition data ----

ovi_bothyears


str(ovi_bothyears)
model_full <- glm(cbind(reproductive, not_reproductive) ~ daylength_hr * temperature * population + year,
                  data = ovi_bothyears,
                  family = binomial)

check_collinearity(model_full)                         # check for collinearity (with no interactions in model)
simulateResiduals(model_full) %>% plot()               # QQ plot & resid vs predicted

summary(model_full)                                    # model results
Anova(model_full, type = 3)
exp(model_full$coefficients)
exp(confint(model_full))
model_full_null <- glm(cbind(reproductive, not_reproductive) ~ 1,
                  data = ovi_bothyears,
                  family = binomial)
1-logLik(model_full)/logLik(model_full_null)
lrtest(model_full, model_full_null)
# plot(not_reproductive ~ daylength_hr, data = ovi_bothyears)


emmeans(model_full, pairwise ~ temperature|population, type = "response")
emtrends(model_full, specs = pairwise ~ temperature|population, var = "daylength_hr")

#prediction
daylength_hr_new <- seq(6, 16, 0.1)
# predict_data <- data.frame(population = rep(c("De", "Sg", "Bi", "ci"), each = 2), temperature = c('38', '28'), year = c('2020', '2022'))

ovi_bothyears_predictnew <- ovi_bothyears %>% dplyr::select(population, temperature, year) %>% distinct() %>%
  slice(rep(1:n(), each = 101)) %>%
  mutate(daylength_hr = rep(daylength_hr_new, 14))


# predict_data <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("population", "temperatue", "year")))) %>%
#   mutate(population = rep(c("De", "Sg", "Bi", "ci"), each = 2), temperature = rep(c('38', '28'), 2))
inverse_logit = function(x){
  exp(x)/(1+exp(x))
}

# ovi_bothyears_predictnew$prediction <- predict(model_full, newdata = ovi_bothyears_predictnew, type = 'response', se = F)
predicted <- as.data.frame(predict(model_full, newdata = ovi_bothyears_predictnew, type = 'link', se.fit = T))
ci_high = inverse_logit(predicted$fit + (predicted$se.fit*1.96))
ci_low = inverse_logit(predicted$fit - (predicted$se.fit*1.96))
predicted_prob = inverse_logit(predicted$fit)
ovi_bothyears_predictnew <- cbind(ovi_bothyears_predictnew, predicted_prob, ci_high, ci_low)



