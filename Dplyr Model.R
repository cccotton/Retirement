cat("\014")
#Dplyr approach - Cary

#Functions
summaryFailureYear <- function (x) {
  balance <- portfolio
  
  runningbalance <- vector(mode = "numeric", length = 0)
  
  for (k in 1:50) {
    balance <- (balance - spending) * x[k]
    runningbalance[k] <- balance
    if (balance < 0) break
  }
  if(k == 50) k <- 99
  return(k)
}

simrun <- function (seed.init, index) {
  set.seed(seed.init)
  simulation <- tbl_df(data.frame(cohort.member = as.numeric(seq(1, 1000)))) %>%
    bind_cols(sample_n(cleanedlifeSpanJoint, 1000, replace = TRUE)) %>%
    mutate(series = sample(1:10000, 1000, replace=T)) %>%
    left_join(solvedMarketScenarios, by = "series") %>%
    mutate(event.year = pmin(death.years, failure.year),
           event.type = ifelse(failure.year < death.years, 1, 0)) %>%
    select(series, event.year, event.type)
  
  simtimes <- tbl_df(data.frame(time = as.numeric(seq(65, 115))))
  
  sfit <- survfit(Surv(as.numeric(simulation$event.year), simulation$event.type) ~ 1)
  sfit.tdf <- tbl_df(data.frame(time = as.numeric(sfit$time), surv = sfit$surv))
  
  crtp <- timepoints(cuminc(as.numeric(simulation$event.year), simulation$event.type, cencode = 2), seq(65, 115))$est
  
  crfit.tdf <- tbl_df(data.frame(index, time = seq(65, 115), cideath = crtp[2,], cifailure = crtp[1,]))
  
  simresults <- simtimes %>%
    left_join(sfit.tdf, by = "time") %>%
    left_join(crfit.tdf, by = "time") %>%
    mutate(surv = ifelse(time == 65, 1, surv),
           surv = ifelse(is.na(surv), 99, surv),
           cideath = cummax(cideath),
           cifailure = cummax(cifailure),
           surv = cummin(surv),
           time = as.numeric(time))
  return(simresults)
}

##Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(survival)
library(cmprsk)

##Files and formatting
marketReturns <- read.csv("Random Market Returns Log-N 10Kx50.csv",header = FALSE)
lifeSpanMale <- read.csv("Male Random Lifetimes from Age 65.csv",header=FALSE)
lifeSpanFemale <- read.csv("Female Random Lifetimes from Age 65.csv",header=FALSE)
lifeSpanJoint <- read.csv("Joint Random Lifetimes from Age 65.csv",header=FALSE)

#marketReturns
marketReturns.tdf <- tbl_df(marketReturns) %>%
  mutate(series = row_number())
cleanedMarketReturns <- marketReturns.tdf %>%
  gather(year, return.annual, -series, convert = T) %>%
  mutate(year = as.numeric(sub("V", "", year))) %>%
  arrange(series, year) %>%
  group_by(series)
cleanedMarketReturns

#lifeSpanMale
lifeSpanMale.tdf <- tbl_df(lifeSpanMale)

cleanedlifeSpanMale <- lifeSpanMale.tdf %>%
  rename(death.years = V1)
cleanedlifeSpanMale

#lifeSpanFemale
lifeSpanFemale.tdf <- tbl_df(lifeSpanFemale)

cleanedlifeSpanFemale <- lifeSpanFemale.tdf %>%
  rename(death.years = V1)
cleanedlifeSpanFemale

#lifeSpanJoint
lifeSpanJoint.tdf <- tbl_df(lifeSpanJoint)

cleanedlifeSpanJoint <- lifeSpanJoint.tdf %>%
  rename(death.years = V1)
cleanedlifeSpanJoint

##Plotting distributions
weight.inverse <- 1 / (max(cleanedMarketReturns$series) * max(cleanedMarketReturns$year))
ggplot(data = cleanedMarketReturns, aes(x = return.annual)) +
  geom_bar(aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Rate of annual returns") + ylab("Proportion of annual returns") +
  theme_classic()

weight.inverse <- 1 / 10000
ggplot(data = cleanedlifeSpanMale, aes(x = death.years)) + ggtitle("Male") +
  stat_bin(binwidth = 1, origin = 64.5, aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Age at death") + ylab("Proportion of population samples") +
  theme_classic()

ggplot(data = cleanedlifeSpanFemale, aes(x = death.years)) +
  stat_bin(binwidth = 1, origin = 64.5, aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Age at death") + ylab("Proportion of population samples") + ggtitle("Female") +
  theme_classic()

ggplot(data = cleanedlifeSpanJoint, aes(x = death.years)) +
  stat_bin(binwidth = 1, origin = 64.5, aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Age at death") + ylab("Proportion of population samples") + ggtitle("Joint") +
  theme_classic()

##Parameters
portfolio <- 1000000
wr = 0.03 # percentage of initial portfolio value to spend annually
spending <- wr * portfolio

##Failure years (conditional on above parameters)
solvedMarketScenarios <- cleanedMarketReturns %>%
  summarise(failure.year = summaryFailureYear(return.annual) + 65)

weight.inverse <- 1 / (max(solvedMarketScenarios$series))
ggplot(data = solvedMarketScenarios, aes(x = failure.year)) +
  geom_bar(aes(weight = weight.inverse)) +
  scale_y_continuous(labels = percent) +
  xlab("Year of Portfolio Failure") + ylab("Proportion of cases (%)") +
  theme_classic()

##Run a single interation of the simulation
simrun(1, 1)

##Run an indexed set of iterations of the simulation
sim.length <- length(seq(65, 115))
first.seed <- 114
n <- 100
index <- seq(1, n)
sims <- vector(mode = "list", length = n)

for (i in 1:n) {
  seed = 114 + i
  sims[[i]] <- simrun(seed.init = seed, index = i)
}

simulation.set <- bind_rows(sims) %>%
  select(-index) %>%
  group_by(time) %>%
  summarise_each(funs(mean, var))

##Graph the mean curves with simulation standard error 95% CL
ggplot(data = simulation.set, aes(x = time, y = surv_mean)) +
  geom_ribbon(aes(ymin = (surv_mean - 1.96 * sqrt(surv_var)), ymax = (surv_mean + 1.96 * sqrt(surv_var))), alpha = 0.05) +
  geom_path() +
  xlab("Age") + ylab("Proportion 1000-Member Cohort Without Portfolio Failure") + xlim(65, 100) + ylim(0, 1) +
  theme_classic()

crsim <- cbind(gather(select(simulation.set, time, cideath_mean, cifailure_mean), curve, mean, -time),
           select(gather(select(simulation.set, time, cideath_var, cifailure_var), curve, var, -time), var))

crsim$curve <- factor(crsim$curve, labels = c("Portfolio failure", "Death"))

ggplot(data = crsim, aes(x = time, y = mean)) +
  geom_path(aes(linetype = curve)) +
  geom_ribbon(aes(ymin = (mean - 1.96 * sqrt(var)),
                  ymax = (mean + 1.96 * sqrt(var)),
                  linetype = curve),
              alpha = 0.05) +
  xlab("Age") + ylab("Proportion 1000-Member Cohort Specified Outcome") + xlim(65, 100) +
  scale_linetype_discrete(name = "") +
  theme_classic()
