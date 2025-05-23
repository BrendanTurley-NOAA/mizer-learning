### MIZER online course
### https://mizer.course.sizespectrum.org/
### May 2025

### Part 1: Understand
### Single species spectra
## https://mizer.course.sizespectrum.org/understand/single-species-spectra.html

## Introduction

# remotes::install_github("sizespectrum/mizerExperimental")

library(mizer)
library(mizerExperimental)
library(tidyverse)


## Single-species model
params <- newSingleSpeciesParams(lambda = 2.05)


## Steady state spectrum
plotSpectra(params, power = 0)
plotSpectra(params, power = 1)
plotSpectra(params, power = 2)


### exercise #1
ex1_params <- newSingleSpeciesParams(lambda = 2.1)
plotSpectra(ex1_params, power = 0)
## end exercise #1


## Numbers
n <- initialN(params)
n
dimnames(n)
n[1, 61]

# It is important to realise that this is not the number of fish in the size class,
# but the average number density in the size class. To get the number of fish we have
# to multiply the number density by the width of the size class. Those widths can 
# be obtained with the dw() function. So the number of fish in each size class is 
# obtained with

numbers <- n * dw(params)
numbers[1, 61]

# You may be surprised by the small number if you interpret it as the number of fish
# between 1 gram and 1.12 gram in the entire ocean. However it looks more reasonable
# if it is the average number per square meter of sea. For more of a discussion 
# this issue of working with numbers per area, numbers per volume or numbers for 
# the entire system see https://sizespectrum.org/mizer/reference/setParams.html#units-in-mizer


### exercise #2
# Determine the total number of fish in the model with sizes between 10 grams and 
# 20 grams. You can use the sum() function to add together contributions from the 
# various size classes.
which(dimnames(n)$w=='10')
which(dimnames(n)$w=='20')

sum(n[1, which(dimnames(n)$w=='10'):which(dimnames(n)$w=='20')])
### end exercise #2

### Other respresentaions
plotSpectra(params)

# Now the green line representing the biomass density of the background has a slope
# of -1.05. The initial slope of the species biomass density is also negative, meaning
# that the biomass density in the species decreases with size.

# We can also plot the biomass density in log weight, i.e., the Sheldon spectrum, 
# supplying the argument power = 2 to plotSpectra().

plotSpectra(params, power = 2)

# This latest plot seems to indicate that most of the biomass of the species is 
# at larger sizes of around 30 grams, whereas the previous plot seemed to indicate
# that most of the biomass is at the smallest sizes. So which one is true? Please
# think about this question, because it really highlights the importance of not 
# confusing biomass density with biomass. Questions about the amount of biomass at
# a size do not make sense. Instead you have to ask about biomass in a size range.

# So for example, we might want to consider the prey biomass available to two different
# predators of our species, one small and one large. Assume that the smaller predator
# feeds on prey in the size range from 1g to 2g. The other predator, which we assume
# is 10 times larger, feeds on prey in the size range from 10g to 20g. These feeding
# intervals have the same width on the logarithmic weight axis. Therefore we should
# look at the plot of the biomass density in log weight to see that the larger predator
# has a lot more prey biomass from our species available to it than the smaller 
# This is in spite of the fact that the plot of biomass density in weight tells us
# that the biomass density is lower at 10g than at 1g.

### Biomass

# We already said above that we can obtain the biomass density in a size class from
# the number density by multiplying the number density by the weight of the individuals
# in the size class. To obtain the appropriate weights, we use the function w() that
# returns the weights at the start of each size class. So we calculate

biomass_density <- n * w(params)

# We obtain the total biomass in each size class by multiplying the biomass density
# in each size class by the width of each size class

biomass <- biomass_density * dw(params)

biomass[61]

# Let us briefly present yet another way to represent the size distribution. When 
# we talk about size spectra, we always have the representation in terms of densities in mind.
# We can similarly describe the size distribution of the biomass by a cumulative
# biomass distribution function, which gives the total biomass of all sizes up to a specific size.

# Initialise an array with the right dimensions
cumulative_biomass <- biomass
# Calculate the cumulative sum of all biomasses in previous bins
cumulative_biomass[] <- cumsum(biomass)
# Normalise this so that it is given as a percentage of the total biomass
cdf <- cumulative_biomass / cumulative_biomass[1, 101] * 100
# Melt the array to a data frame and then plot
p_biomass_cdf <- ggplot(melt(cdf), aes(x = w, y = value)) +
  geom_line() + 
  labs(x = "Weight [g]",
       y = "% of total biomass")
p_biomass_cdf

p_biomass_cdf +
  scale_x_log10()


### Allometric rates

# In a multi-species mizer model the mortality is an emergent property that depends
# on the abundance of predators. In this single species model the mortality rate is
# set to the one that would emerge if all the species in the fixed background community
# predated with the same ferocity as the target species. This leads to a rate that 
# scale allometrically with a power of n−1=3/4−1=−1/4. This means that the death 
# rate experienced by larger individuals is smaller than that of smaller individuals.

growth_rate <- getEGrowth(params)
growth_rate[1, 61]

# This is the instantaneous per-capita growth rate, measured in grams per year. Note
# that in mizer all rates are measured in units of 1/year, but for many people daily
# values are easier to understand. Since growth rate here is an instantaneous rate
# we can simply divide it by 365 to get a daily rate (although note that mizer does
# not simulate processes on daily time steps). This gives us a growth rate in grams
# per day for a 1g sized fish of

growth_rate[1, 61] / 365
growth_rate_frame <- melt(growth_rate)
names(growth_rate_frame)

p <- ggplot(growth_rate_frame) +
  geom_line(aes(x = w, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Growth rate [g/year]")
p

# We see that at least up to a size of a few grams the line is straight. Let’s isolate
# the growth rate for those smaller sizes

g_small_fish <- filter(growth_rate_frame, w <= 10)
lm(log(g_small_fish$value) ~ log(g_small_fish$w))


### Excercie #3
# Use the methods you have just seen to make a log-log plot of the mortality rate.
# You can get the mortality rate with the getMort() function. While adjusting 
# code to this new task, you need to take into account that the name of the size-
# dimension of the array returned by getMort() is "w_prey" instead of "w".

# Then fit a linear model to determine the slope and and intercept and thus the 
# allometric exponent r and the coefficient μ0 for the mortality rate

mort_rate <- getMort(params)
mort_rate_frame <- melt(mort_rate)
names(mort_rate_frame)

p <- ggplot(mort_rate_frame) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Mortality rate [1/year]")
p

mm <- lm(log(mort_rate_frame$value) ~ log(mort_rate_frame$w_prey))
mm

m0 <- exp(coef(mm)[[1]])
m0
### end of exercise #3

### Slope of juvenile spectrum
# We have seen that for juvenile fish the growth rate and the death rate are both
# power laws. By solving a differential equation we can derive that the juvenile 
# spectrum also follows a power law. Number density drops off faster with size 
# if the mortality rate coefficient μ0 is higher or if the growth rate coefficient 
# g0 is smaller, which is what we would expect.

# We can also check this claim numerically. Let’s look at the spectrum of individuals
# up to 10 grams. By now we know how to do this. We first convert the number density
# matrix n into a dataframe and then filter out all observations that do not have 
# w ≤ 10. The resulting data frame we pass to ggplot() and ask it to plot a line on log-log axes.

nf <- melt(n) %>% 
  filter(w <= 10)

ggplot(nf) +
  geom_line(aes(x = w, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Number density [1/g]")

lm(log(nf$value) ~ log(nf$w))

g0 <- 4.2 
-m0 / g0 - 3/4

# it is also amazing how we can calculate expected numbers of fish from basic assumptions
# and rules. Of course, natural ecosystems never look like that, but if we have theoretical
# expectations derived from clear assumptions (about growth and mortality rate and
# food availability), we can start asking questions about which processes in natural
# ecosystems deviate from these generic assumptions, why this happens and how it 
# should affect the observed size spectra.


### Shape of the adult spectrum
