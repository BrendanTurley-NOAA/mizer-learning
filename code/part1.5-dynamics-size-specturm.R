### MIZER online course
### https://mizer.course.sizespectrum.org/
### June 2025

### Part 1: Understand
### Dynamics of size spectra
## https://mizer.course.sizespectrum.org/understand/dynamics-of-spectra.html#size-spectrum-dynamics

### Size spectrum dynamics
# In previous tutorials we have concentrated on the steady state of the mizer model,
# where for each size class and each species, the rate at which individuals grow into
# the size class balances the rate at which individuals grow out of the size class or
# die, thus keeping the size spectrum constant. In this tutorial we explore the dynamic 
# that takes place when this balance is changed.

# Size spectrum dynamics is very intuitive: The rate at which the number of individuals 
# in a size class changes is the difference between the rate at which individuals 
# are entering the size class and the rate at which they are leaving the size class.
# Individuals can enter a size class by growing in or, in the case of the smallest 
# size class, by being born into it. They can leave by growing out or by dying.

library(mizer)
library(mizerExperimental)
library(tidyverse)

## Projections
# In the previous tutorial, in the section on trophic cascades, we already simulated
# the size-spectrum dynamics to find the new steady state. But we only looked at the
# final outcome once the dynamics had settled down to the new steady state. We reproduce
# the code here:

# Create trait-based model
mp <- newTraitParams() |> 
  # run to steady state with constant reproduction rate
  steady() |>
  # turn of reproduction and instead keep egg abundance constant
  setRateFunction("RDI", "constantEggRDI") |>
  setRateFunction("RDD", "noRDD")

# We make a copy of the model
mp_lessRes <- mp
# and set the resource interaction to 0.8 for species 8 to 11
given_species_params(mp_lessRes)$interaction_resource[8:11] <- 0.8

# We run the dynamics until we reach steady state
mp_lessRes_steady <- projectToSteady(mp_lessRes)

# We compare the steady states
plotSpectra2(mp_lessRes_steady, name1 = "less resource", 
             mp, name2 = "original",
             total = TRUE, power = 2,
             ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# 
# But we can also save and then display the spectra of all the species at intermediate 
# times. This is what the project() function does. It projects the current state forward
# in time and saves the result of that simulation in a larger object, a MizerSim object,
# which contains the resulting time series of size spectra. Let’s use it to project 
# forward by 24 years.

sim_lessRes <- project(mp_lessRes, t_max = 24)

animateSpectra(sim_lessRes, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# Note, for some species size spectra at the largest size class drop all the way to
# very small values (e.g. 10^-7) and for others they stop higher. This is just a
# discretisation artefact and is not important. Try to ignore it.


getTimes(sim_lessRes)

N(sim_lessRes)[6, 2, 1]

getBiomass(sim_lessRes)[6, 2]

plotBiomass(sim_lessRes)


## Reproduction dynamics
# The above simulation was run with constant abundance in the smallest size class 
# for each species. This of course is not realistic. The abundance of the smallest 
# individuals depends on the rate at which mature individuals spawn offspring, and 
# this in turn depends, among other things, on the abundance of mature individuals. 
# So if the abundance of mature individuals goes down drastically, as it did for 
# species 8 to 11 above, then the abundance of offsprings for those species will 
# go down as well.

# To see the effect we run the same code as above after deleting the two lines that
# turned off the reproduction dynamics. We also specify with t_save = 2 that we want
# to save the spectrum only every second year, which speeds up the display of the 
# animation.

# Create trait-based model and run to steady state
mp <- newTraitParams() |> steady()

# We make a copy of the model
mp_lessRes <- mp
# and set the resource interaction to 0.8 for species 8 to 11
species_params(mp_lessRes)$interaction_resource[8:11] <- 0.8

# We simulate the dynamics for 30 years, saving only every 2nd year
sim_lessRes <- project(mp_lessRes, t_max = 30, t_save = 2)

# We animate the result
animateSpectra(sim_lessRes, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))


# Note that now the fish species whose access to resources was decreased continue 
# to decrease in abundance as time goes on. Interestingly, species 1 appears to be 
# driven to extinction by the increased abundance of its predators that in turn is 
# due to the decrease in predation from the larger species.


## Energy invested into reproduction
# We already discussed the investment into reproduction in an earlier tutorial. As
# mature individuals grow, they invest an increasing proportion of their income into
# reproduction and at their asymptotic size they would be investing all income into
# reproduction. Summing up all these investments from mature individuals of a particular
# species gives the total rate ER at which that species invests energy into reproduction.
# This total rate of investment is multiplied by a reproduction efficiency factor 
# erepro, divided by a factor of 2 to take into account that only females reproduce,
# and then divided by the egg weight w_min to convert it into the rate at which eggs
# are produced. 

getRDI(mp)

# The erepro parameter or reproduction efficiency can vary between 0 and 1 (although
# 0 would be bad) and gives the proportion of energy invested into reproduction that
# is converted into viable eggs or larvae.


## Density-dependence in reproduction
# The stock-recruitment relationship is an emergent phenomenon in mizer, with several
# sources of density dependence. Firstly, the amount of energy invested into reproduction
# depends on the energy income of the spawners, which is density-dependent due to 
# competition for prey. Secondly, the proportion of larvae that grow up to recruitment
# size depends on the larval mortality, which depends on the density of predators, 
# and on larval growth rate, which depends on density of prey.

# However, there are other sources of density dependence that are not explicitly 
# modelled mechanistically in mizer. An example would be a limited carrying capacity
# of suitable spawning grounds and other spatial effects. So mizer has another species
# parameter Rmax that gives the maximum possible rate of recruitment. Imposing a finite 
# maximum reproduction rate leads to a non-linear relationship between energy invested 
# and eggs hatched. This density-dependent reproduction rate Rdd is given by a Beverton-Holt
# type function.

getRDD(mp)


### Reproduction level
# We have seen the two species parameters that determine how the energy invested 
# into reproduction is converted to the number of eggs produced: erepro and R_max.
# For neither of these is it obvious what value they should have. The choice of values
# influences two important properties of a model: the steady state abundances of the
# species and the density-dependence in reproduction. It is therefore useful to change
# to a new set of two parameters that reflect these two properties better. These are:
  
# The birth rate Rdd at steady state. This determines the abundance of a species.

# The ratio between Rdd and Rmax at steady state. This determines the degree of density dependence.

# The ratio Rdd/Rmax we denote as the reproduction level. This name may remind you 
# of the feeding level, which was the ratio between the actual feeding rate and the
# maximum feeding rate and described the level of density dependence coming from satiation.
# It takes a value between 0 and 1. It follows from our discussion in the previous section 
# that a species with a high reproduction level is more resilient to changes.

getReproductionLevel(mp)

# We see that by default newTraitParams() had given all species the same reproduction 
# level. We can change the reproduction level with the setBevertonHolt() function. 
# We can set different reproduction levels for each species, but here we will simply
# set it to 0.9 for all species:

mp9 <- setBevertonHolt(mp, reproduction_level = 0.9)

# Changing the reproduction level has no effect on the steady state, because that 
# only depends on the rate of egg production Rdd and that is kept fixed when changing
# the reproduction level. We can check that by running our new model to steady state 
# and plotting that steady state together with the original steady state.

mp9 <- projectToSteady(mp9)
plotSpectra2(mp, name1 = "reproduction_level = 0.25",
             mp9, name2 = "reproduction_level = 0.9",
             total = TRUE, power = 2, ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# They overlap perfectly. However the reproduction level does have an effect on 
# sensitive the system is to changes. As an example, let us look at the dynamics that
# is triggered by the reduction in interaction with the resource by species 8 through 11.

# We make a copy of the model
mp_lessRes9 <- mp9
# and set the resource interaction to 0.8 for species 8 to 11
species_params(mp_lessRes9)$interaction_resource[8:11] <- 0.8

sim_lessRes9 <- project(mp_lessRes9, t_max = 30, t_save = 2)

animateSpectra(sim_lessRes9, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# Notice how the species have settled down to a new steady state after 30 years without
# any extinctions and the impact on species 1 is much less extreme. As expected, the
# higher reproduction level has made the species more resilient to perturbations.

# The problem of course is that in practice the reproduction level is hardly ever 
# known. Instead, one will need to use any information one has about the sensitivity
# of the system from observed past perturbations to calibrate the reproduction levels.
# We’ll discuss this again towards the end of Part 2.

### Excercise #1
Go back to the example with fishing on individuals above 1kg from the section on 
fishing-induced cascades. Impose the same fishing, but now on the trait-based model
with reproduction dynamics left turned on and with a reproduction level of 0.5 for 
all species. Project the model for 20 years and animate the result.


mp5 <- setBevertonHolt(mp, reproduction_level = 0.5)
initial_effort(mp5) <- 1

mp_lessRes5 <- mp5
species_params(mp_lessRes5)$interaction_resource[8:11] <- 0.8

sim_lessRes5 <- project(mp_lessRes5, t_max = 20, t_save = 2)

animateSpectra(sim_lessRes5, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

# We make a copy of the trait-based model
mp_fishing <- mp
# Set the reproduction level set to 0.5 for all species
mp_fishing <- setBevertonHolt(mp, reproduction_level = 0.5)
# Set fishing effort to 1
initial_effort(mp_fishing) <- 1

# Simulate the dynamics for 20 years
sim_fishing <- project(mp_fishing, t_max = 20, t_save = 2)

# Animate the result
animateSpectra(sim_fishing, total = TRUE, power = 2, 
               ylim = c(1e-8, NA), wlim = c(1e-3, NA))

### end 


## Resource dynamics
# The resource spectrum is not described by size spectrum dynamics, because in reality
# it is typically not made up of individuals that grow over a large size range during
# their life time. In mizer, the resource number density in each size class is describe
# by semichemostat dynamics: the resource number density in each size class recovers
# after depletion, and this biomass growth or recovery rate will decrease as the number
# density gets close to a carrying capacity. If you want the mathematical details,
# you can find them in the mizer model description in the section on resource density.

# The effect of these dynamics is that if the number of fish consuming the resource
# in a certain size range increases, the resource abundance in that size range will
# decrease, if it cannot recover quickly enough (regeneration rate of the resource 
# is set by the user). So there is competition for the resource, which provides a
# stabilising influence on the fish abundances. We will be discussing this more in 
# later tutorials.


