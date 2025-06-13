### MIZER online course
### https://mizer.course.sizespectrum.org/
### May 2025

### Part 1: Understand
### Species interactions
## https://mizer.course.sizespectrum.org/understand/species-interactions.html

# In the previous tutorials we studied a single species interacting with a fixed 
# background community. In this tutorial we want to acknowledge that there is no 
# such thing as a fixed background community. Instead, all species form part of a 
# dynamical ecosystem in which changes to any species has knock-on effects on other
# species.

library(mizer)
library(mizerExperimental)
library(tidyverse)


## Trait-based model
# In this first part of the course we aim for understanding, not realism. So in this
# tutorial we investigate the tangled web of interactions in an idealised multi-species
# system. We choose a trait-based model in which the species making up the community
# differ from each other only in a single trait: their asymptotic body size (sometimes
# it is also called maximum body size).

# We use the newTraitParams() function to create our idealised trait-based multi-
# species model. The function has many parameters, but we will just keep the defaults.
# Unlike the newSingleSpeciesParams() function, the newTraitParams() function does not set
# the initial spectra to their steady state values. We thus need to run the result 
# through the steady() function. We assign the resulting MizerParams object to the 
# variable mp.

mp <- newTraitParams() |> steady()

plotSpectra(mp, power = 2, total = TRUE)

## Turn off reproduction dynamics
# As in previous tutorials, we want to concentrate on the shapes of the size spectra
# and we do not yet want to look at what determines the overall abundance of each species.
# Therefore we modify the model so that it keeps the abundances at egg size fixed 
# (i.e. numbers in the first size bin).

mp <- mp |>
  setRateFunction("RDI", "constantEggRDI") |>
  setRateFunction("RDD", "noRDD")


## Mortality from other species

plotlyDeath(mp, species = "8")

plotDeath(mp, species = "8", proportion = FALSE)


## Income from other species

plotDiet(mp, species = "8")

# The diet composition we see in the plot is shaped by two things: the predation 
# kernel (the size preference in the feeding of the predator) and the relative 
# abundances of prey at different sizes. First, a predator will only eat food that
# is within the predation kernel size range. But once in this size range the relative
# proportion of different species or resource consumed will simply depend on their 
# relative biomass. So if, for example, 80% of biomass in a specific prey size class
# consists of resource, 15% of species 1 and 5% of species 2, then the diet of the
# predator feeding in that size class will consist of 80% resource, 15% of species
# 1 and 5% of species 2.

# Of course, when we build a model for a real-world ecosystem we will have some 
# knowledge about the biology of the species and their food preferences. Perhaps 
# one species is actively selecting fish out of the resource, or predating on specific
# species only? This is where the interaction matrix comes in that we will discuss 
# in the next section.


### Excercise #1
# Now, check what the diets of other species look like by making a plot for each 
# species from 1 (smallest one) to 11 (the largest one). Hint: if you look at the
# documentation for plotDiet() you may find a convenient way to do this.

plotDiet(mp)

### end

## Interaction matrix
# Now we arrive to an interesting and challenging aspects of multi-species modelling
# - setting up parameters for species and resource interactions. By default, mizer 
# assumes that all species in the model can access all other species and resource 
# equally and the amount of different prey consumed just depends on their relative 
# abundance in the predator’s feeding size range. So the default interaction matrix
# of the species in our model looks very simple

interaction_matrix(mp)

# The matrix has all values set at 1 which means that all predators can access all
# prey species equally.

mp_modified <- mp
# We change row 11 (predator species 11) and column 2 (prey species 2) 
# to a smaller value
interaction_matrix(mp_modified)[11, 2] <- 0.2

plotDeath(mp, species = 2)
plotDeath(mp_modified, species = 2)

plotDiet(mp, species = 11)
plotDiet(mp_modified, species = 11)


## Resource interactions

species_params(mp)$interaction_resource


# Now we might want to reduce the availability of resource to some predators. Perhaps
# we know that certain species much prefer to feed on other fish rather than on similar
# sized plankton. Let us look at an example where species 8 through 11 have a 20% 
# reduction in their interaction with resource.


# We make a copy of the model
mp_lessRes <- mp
# and set the resource interaction to 0.8 for species 8 to 11
given_species_params(mp_lessRes)$interaction_resource[8:11] <- 0.8
# We print out the result to check
species_params(mp_lessRes)$interaction_resource

plotDiet(mp, species = 9)
plotDiet(mp_lessRes, species = 9)


# The change seems small enough. However, now that we changed the availability of 
# resources, which is so important for larval stages, these four species will experience
# a much reduced growth rate during their juvenile stage. We can see that effect by
# recalculating the single-species spectra with

mp_lessRes_sss <- steadySingleSpecies(mp_lessRes)
plotSpectra(mp_lessRes_sss, power = 2, total = T)


## Trophic cascades
# As we just discussed, the above picture does does not show a steady state of the
# ecosystem. Species now find themselves with a different abundance of predators and 
# prey and this will change their mortality and their growth and hence their size spectra.

# The easiest way to find the new steady state that the ecosystem will settle into 
# is to simulate the full multi-species dynamics forward in time. Mizer refers to 
# this simulation to find the future state of the ecosystem as “projecting”. We can
# use the function projectToSteady() to project forward in time far enough so the 
# system has settled down again close to the new steady state.

mp_lessRes_steady <- projectToSteady(mp_lessRes)
plotSpectra2(mp_lessRes_steady, name1 = "less resource", 
             mp, name2 = "original", 
             total = TRUE, power = 2, ylim = c(1e-8, NA), wlim = c(1e-3, NA))


## Fishing-induced cascades
# Let’s investigate these trophic cascades a bit more. This time we can look at 
# fishing large fish will affect the ecosystem.

# The model has been set up with a knife-edge fishing gear that selects all individuals
# above 1kg, irrespective of species. This is not a realistic gear and mizer can do 
# much better, as we will see in Part 3. But it serves our current purpose, because 
# it will impose a fishing mortality that only impacts the larger species that actually
# grow to sizes above 1kg. To use that gear we just have to set a non-zero fishing 
# effort. We create a new model mp_fishing with a fishing effort of 1:

mp_fishing <- mp
initial_effort(mp_fishing) <- 1

mp_fishing_sss <- steadySingleSpecies(mp_fishing)
plotSpectra(mp_fishing_sss, power = 2, total = T)


### Excerise #2
# Project the mp_fishing model to its steady state and then make a plot comparing 
# it to the steady state of the un-fished system. Do you see a trophic cascade?

mp_fishing_steady <- projectToSteady(mp_fishing)
plotSpectra2(mp_fishing_steady, name1 = "fishing", 
             mp, name2 = "original", 
             total = TRUE, power = 2, ylim = c(1e-6, NA), wlim = c(1e-2, NA))

