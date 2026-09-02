library(dplyr)
library(stringr)

oviposition.dat <- read.csv("C://Users//eliza//Google Drive//GRAD SCHOOL//RESEARCH//D. carinulata//temperature//analyses from amanda//diapause_datasheet_final_data.csv",
                            header = T, strip.white = T, na.strings = "NA")
oviposition.dat <- oviposition.dat[-which(oviposition.dat$eggs_present_date == ""),]
oviposition.dat <- oviposition.dat[,-which(colnames(oviposition.dat) == "Notes")]
oviposition.dat <- oviposition.dat[,-which(colnames(oviposition.dat) == "mortality")]
oviposition.dat <-  oviposition.dat[complete.cases(oviposition.dat),]

# View(oviposition.dat)

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

library(ggplot2)
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



