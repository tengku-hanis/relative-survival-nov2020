#======================================================================================#
# Title: An introduction to relative survival analysis in R
# Author: Tengku Hanis 
# Date: 05-11-2020
#======================================================================================#
##----- First part: Expanding abridged life table

## Package
library(MortalityLaws)
library(tidyverse)
library(relsurv)

## Abridged life table (from github)
slovenia_male <- read.table("https://raw.githubusercontent.com/tengku-hanis/relative-survival-nov2020/main/slovenia_males.txt", 
                            header = TRUE, quote="\"")
slovenia_female <- read.table("https://raw.githubusercontent.com/tengku-hanis/relative-survival-nov2020/main/slovenia_females.txt", 
                              header=TRUE, quote="\"")

## Expand male
age_int <- c(0,1,seq(5, 110, by=5)) # age interval in abridged life table
age_range <- 0:110 # range of age to be expanded

# filter 1994-2005
slovenia_male$Year <- as.factor(slovenia_male$Year)
by_yearM <- slovenia_male %>% filter(Year %in% 1994:2005) %>% group_by(Year) %>% nest()

# Separate mx in list
mx_male <- vector("list", 0)
for (i in seq_along(by_yearM$data)) {
  mx_male[[i]] <- by_yearM$data[[i]]$mx
}

# Estimate coefficient
models_male <- vector("list", 0)
for (i in seq_along(mx_male)) {
  models_male[[i]] <- MortalityLaw(age_int, mx = mx_male[[i]], law = "siler", opt.method = "LF2") 
}

# Expand life table
male_1994_2005 <- vector("list", 0)
for (i in seq_along(models_male)) {
  male_1994_2005[[i]] <- LawTable(age_range, par = models_male[[i]]$coefficients, 
                                  law = "siler", sex = "male")
}

# Combine life table into data frame
male_list <- vector("list", 0)
for (i in seq_along(male_1994_2005)) {
  male_list[[i]] <- data.frame(male_1994_2005[[i]]$lt)
}
male_lt <- male_list %>% enframe() %>% unnest(cols = value)

## Expand female
# age interval and range of age as male

# filter 1994-2005
slovenia_female$Year <- as.factor(slovenia_female$Year)
by_yearF <- slovenia_female %>% filter(Year %in% 1994:2005) %>% group_by(Year) %>% nest()

# Separate mx in list
mx_female <- vector("list", 0)
for (i in seq_along(by_yearM$data)) {
  mx_female[[i]] <- by_yearM$data[[i]]$mx
}

# Estimate coeffcient
models_female <- vector("list", 0)
for (i in seq_along(mx_female)) {
  models_female[[i]] <- MortalityLaw(age_int, mx = mx_female[[i]], law = "siler", opt.method = "LF2") 
}

# Expand life table
female_1994_2005 <- vector("list", 0)
for (i in seq_along(models_female)) {
  female_1994_2005[[i]] <- LawTable(age_range, par = models_female[[i]]$coefficients, 
                                    law = "siler", sex = "female")
}

# Combine life table into data frame
female_list <- vector("list", 0)
for (i in seq_along(female_1994_2005)) {
  female_list[[i]] <- data.frame(female_1994_2005[[i]]$lt)
}
female_lt <- female_list %>% enframe() %>% unnest(cols = value)

## Make a rate table 
# select age(x) and probability of dying(qx)
# px (survival probability) = 1-qx
pop_m <- male_lt %>% mutate(year = rep(1994:2005, each = 111), px = 1-qx) %>% 
  select(x, year, px)
pop_f <- female_lt %>% mutate(year = rep(1994:2005, each = 111), px = 1-qx) %>% 
  select(x, year, px)

# long -> wide
# use px (survival probability)
pop_m.w <- pivot_wider(pop_m, names_from = year, values_from = px)
pop_f.w <- pivot_wider(pop_f, names_from = year, values_from = px)
str(pop_m.w)
str(pop_f.w)

# delete age column (x)
pop_m.w$x <- NULL
pop_f.w$x <- NULL

# as matrix
pop_fwm <- as.matrix(pop_f.w)
pop_mwm <- as.matrix(pop_m.w)
str(pop_fwm)
str(pop_mwm)

# ratetable
pop_rate <- transrate(men = pop_mwm, women = pop_fwm, yearlim = c(1994, 2005), int.length = 1)
is.ratetable(pop_rate)
str(pop_rate)

##----- Second part: Non-parametric relative survival

## Packages
library(relsurv)
library(survminer)

## Data
?colrec

# limit follow-up time to 5 years 
colrec$end <- colrec$diag + colrec$time
  # recode end time
colrec$end2 <- ifelse(colrec$end > as.date("31Dec2005", order = "dmy"), 
                      as.date("31Dec2005", order = "dmy"), colrec$end)
  # recode the status
colrec$stat2 <- ifelse(colrec$end > as.date("31Dec2005", order = "dmy"), 0, colrec$stat)
  # edit follow-up time
colrec$time2 <- colrec$end2 - colrec$diag

## Pohar-perme
rs_PP <- rs.surv(Surv(time2, stat2) ~ 1, 
                 rmap = list(age = age, sex = sex, year = diag), 
                 method = "pohar-perme",
                 ratetable = slopop,
                 data = colrec)
summary(rs_PP, times = 1:12 * 365.241)

# Plot 
plot(rs_PP, xlab = "Time(years)", ylab = "Relative survival", xscale = 365.24, 
     main="Pohar-perme method")

# Plot using survminer
ggsurvplot(rs_PP, conf.int = T, xscale="d_y", break.x.by= 365.24, xlab= "Time(years)", 
           title= "Pohar-perme method", censor=F, legend = "none")

## Other estimators
# Edere 1
rs_e1 <- rs.surv(Surv(time2, stat2) ~ 1, 
                 rmap = list(age = age, sex = sex, year = diag), 
                 method = "ederer1",
                 ratetable = slopop,
                 data = colrec)
# Ederer2
rs_e2 <- rs.surv(Surv(time2, stat2) ~ 1, 
                 rmap = list(age = age, sex = sex, year = diag), 
                 method = "ederer2",
                 ratetable = slopop,
                 data = colrec)
# Hakulinen
rs_h <- rs.surv(Surv(time2, stat2) ~ 1, 
                rmap = list(age = age, sex = sex, year = diag), 
                method = "hakulinen",
                ratetable = slopop,
                data = colrec)

## Compare plot of all estimators
rs_list <- list("Pohar-perme" = rs_PP, "Ederer1" = rs_e1, "Ederer2" = rs_e2, 
                "Hakulinen" = rs_h)
ggsurvplot(rs_list, data = colrec, conf.int = T, combine = T,  xscale="d_y", 
           break.x.by= 365.24, censor = F,
           xlab= "Time(years)", ylab="Relative survival",  
           title= "Relative survival", 
           legend="top", legend.title= "Estimators", 
           legend.labs=c("Pohar-perme", "Ederer1", "Ederer2", "Hakulinen"))

## Compare expanded life table and non-expanded life table
# Pohar-perme with expanded life table
rs_PP_expand <- rs.surv(Surv(time2, stat2) ~ 1, 
                        rmap = list(age = age, sex = sex, year = diag), 
                        method = "pohar-perme",
                        ratetable = pop_rate,
                        data = colrec)
summary(rs_PP_expand, times = 1:12 * 365.241)


plot(rs_PP_expand, xlab = "Time(years)", ylab = "Relative survival", xscale = 365.24, 
     main="Pohar-perme method")

# Compare plot
compare_lt <- list("Non-expanded" = rs_PP, "Expanded" = rs_PP_expand)
ggsurvplot(compare_lt, data = colrec, conf.int = T, combine = T,  xscale="d_y", 
           break.x.by= 365.24, censor = F,
           xlab= "Time(years)", ylab="Relative survival",  
           title= "Relative survival using Pohar-Perme", 
           legend="top", legend.title= "Life table:", 
           legend.labs=c("Non-expanded", "Expanded"))

# Data of the difference between survival estimate
diff_surv <-rs_PP$surv - rs_PP_expand$surv
dat_diff <- data.frame(day = rs_PP_expand$time, diff = diff_surv) 
dat_diff$year <- dat_diff$day/365.241

# Plot the difference
ggplot(dat_diff, aes(year, diff)) + geom_point() + 
  scale_x_continuous(breaks = 0:12) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Relative survival different between rate table", 
       subtitle = "(Non-expand - Expand)",
       x = "Years", y = "Relative survival different") +
  theme_bw()

## Compare relative survival between sites
diff_site <- rs.diff(Surv(time2, stat2) ~ site, 
                     rmap = list(age = age, sex = sex, year = diag), 
                     data = colrec, ratetable = slopop)
diff_site 

# Plot
diff_rs <- rs.surv(Surv(time2, stat2) ~ site, 
                   rmap = list(age = age, sex = sex, year = diag), 
                   data = colrec, ratetable = slopop)
ggsurvplot(diff_rs, data = colrec, conf.int = T, combine = T,  xscale="d_y", 
           break.x.by= 365.24, censor = F,
           xlab= "Time(years)", ylab="Relative survival",  
           title= "Relative survival using Pohar-Perme", 
           legend="top", legend.title= "Site:", 
           legend.labs=c("Colon", "Rectum"))

## Crude probability of death and year lost
cpdeath <- cmp.rel(Surv(time2, stat2) ~ site, 
                   rmap = list(age = age, sex = sex, year = diag), 
                   ratetable = slopop, data = colrec, tau = 3652.41)
  # tau = 3652.41, value after 10 years is censored
summary(cpdeath, times = c(1, 5, 10), scale = 365.24) #scale default 1:day
  # within 10 years after diagnosis 58.4% of patients have died dt 
  # colorectal cancer at rectum
  # and 54.5% of patients have died dt colorectal cancer at colon
cpdeath
  # area=year lost dt disease
  # patients with colorectal cancer at rectum lost 4.8 years
  # patients with colorectal cancer at colon lost 4.8 years

plot(cpdeath, xscale = 365.24 ,col = 1:4, conf.int = c(3, 1),
     xlab = "Time(years)", main = "Crude probabily of death")


