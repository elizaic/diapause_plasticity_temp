library(dplyr)
library(stringr)
library(car)
library(emmeans)

oviposition.dat <- read.csv("diapause_datasheet_final_data.csv",
                            header = T, strip.white = T, na.strings = "NA")
oviposition.dat <- oviposition.dat[-which(oviposition.dat$eggs_present_date == ""),]
oviposition.dat <- oviposition.dat[,-which(colnames(oviposition.dat) == "Notes")]
oviposition.dat <- oviposition.dat[,-which(colnames(oviposition.dat) == "mortality")]
oviposition.dat <-  oviposition.dat[complete.cases(oviposition.dat),]
oviposition.dat$pair_ID <- oviposition.dat$?..pair_ID

View(oviposition.dat)

oviposition.dat$pop <-
  substr(oviposition.dat$pair_ID, 1, 4) %>%
  as.factor()

oviposition.dat$pop_treatment_combo <- paste0(oviposition.dat$pop, "_", oviposition.dat$treatment_combo)

oviposition.dat <- mutate(oviposition.dat, photoperiod = str_split(treatment_combo, "_", simplify = T)[,1], 
         temp = str_split(treatment_combo, "_", simplify = T)[,2]) 

summary_ovi <- oviposition.dat %>%
  group_by (pop_treatment_combo) %>% 
  summarise(diapause = sum(eggs_present_date=="-"),
            egglaid =  sum(eggs_present_date!="-"),
            reproductive_rate = egglaid/(diapause+egglaid)) %>%
  as.data.frame()

summary_ovi_meta <- data.frame(str_split(summary_ovi$pop_treatment_combo, "_", simplify = T))
colnames(summary_ovi_meta) <- c("pop", "photoperiod", "temp")
summary_ovi <- cbind(summary_ovi, summary_ovi_meta)
summary_ovi$pop_f = factor(summary_ovi$pop, levels=c("De20", "Sg20", "Ci20"))

ggplot(data = summary_ovi) +
  geom_point(aes(x = photoperiod, y = reproductive_rate, 
                 shape = temp, color = temp), alpha = .8, size = 4) +
  geom_path(aes(x = photoperiod, y = reproductive_rate, group = temp, color = temp), size = 1) +
  # geom_smooth(aes(x = photoperiod, y = reproductive_rate, group = temp), 
  #             method = "glm", 
  #             method.args = list(family = "binomial"), 
  #             se = FALSE) +
  scale_color_manual(values = c("sky blue", "orange")) +
  theme_bw() +
  # facet_grid(~pop_f)
facet_grid(rows = vars(pop_f))


ggplot(data = summary_ovi) +
  geom_point(aes(x = temp, y = reproductive_rate, 
                 shape = pop, color = pop), alpha = .8, size = 4) +
  geom_path(aes(x = temp, y = reproductive_rate, group = pop, color = pop), size = 1) +
  # geom_smooth(aes(x = photoperiod, y = reproductive_rate, group = temp), 
  #             method = "glm", 
  #             method.args = list(family = "binomial"), 
  #             se = FALSE) +
  # scale_color_manual(values = c("sky blue", "orange")) +
  theme_bw() +
  # scale_x_discrete(expand = c(0, 0.5)) +
  # facet_grid(~pop_f)
  facet_grid(rows = vars(photoperiod))

cool <- summary_ovi %>%
  dplyr::filter(temp == "28/13")%>%
  select(pop, photoperiod, reproductive_rate) %>%
  rename(cool_reproductive_rate=reproductive_rate)

warm <- summary_ovi %>%
  dplyr::filter(temp == "38/23") %>%
  select(pop, photoperiod, reproductive_rate) %>%
  rename(warm_reproductive_rate=reproductive_rate)

plasticity <- merge(cool, warm) %>%
  mutate(temp_diff = warm_reproductive_rate - cool_reproductive_rate)


ggplot(plasticity) +
  geom_col(aes(x = photoperiod, y = temp_diff, 
                 fill = pop), position = "dodge", alpha = .8, size = 4) +
  ylab("proportion reproproductive in warm - cool treamtements") +
  theme_bw() 

ggplot(plasticity) +
  geom_hline(yintercept =0) +
  geom_point(aes(x = photoperiod, y = temp_diff, 
               color = pop, shape = pop), alpha = .8, size = 4) +
  geom_path(aes(x = photoperiod, y = temp_diff, group = pop, color = pop), size = 1) +
  ylab("proportion reproproductive in warm - cool treamtements") +
  theme_bw() 



glm(formula = summary_ovi$egglaid ~ summary_ovi$photoperiod * 
      summary_ovi$temp + summ)

model1 <- glm(cbind(egglaid, diapause) ~ pop * temp + photoperiod, family = binomial, data = summary_ovi)
summary(model1)
pchisq(model1$deviance, df = 25, lower.tail = F)
plot(residuals(model1, "pearson")~ fitted(model1))
Anova(model1, type = 3)
emmeans(model1, pairwise ~temp|pop, at = list(photoperiod = "11:25"), type = "response")
emmeans(model1, ~temp|pop, at = list(photoperiod = "14:10"), type = "response")


library(MASS)
dose.p(model1, cf = c(1,3), p= 0.5)
