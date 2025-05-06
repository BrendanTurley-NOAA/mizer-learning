### MIZER online course
### https://mizer.course.sizespectrum.org/
### May 2025

### Part 1: Understand
### Observed size spectrum


### Introduction

# install.packages(c("mizer", "tidyverse", "plotly", "remotes", "usethis",
#                    "rmarkdown", "rstudioapi"))
library(tidyverse)


## The data

# download.file("https://github.com/gustavdelius/mizerCourseNew/raw/master/understand/size-data.rds",
              # destfile = "data/size-data.rds")

size_data <- readRDS("data/size-data.rds")

str(size_data)
unique(size_data$species)


## Histogram

hist(size_data$weight, xlab = 'weight [g]', ylab = 'number of fish', las = 1, nclass = 30)
plot(hist(size_data$weight)$mids, hist(size_data$weight)$counts, log="y", type='h', lwd=10, lend=2)

p <- ggplot(size_data) +
  geom_histogram(aes(weight), fill = "blue", colour = "black") +
  labs(x = "Weight [g]",
       y = "Number of fish")
p
p + scale_y_log10()


log_breaks <- seq(from = 0, to = 11, by = 1)
log_breaks

bin_breaks <- 2 ^ log_breaks
bin_breaks

p2 <- ggplot(size_data) +
  geom_histogram(aes(weight), fill = "blue", colour = "black",
                 breaks = bin_breaks) +
  labs(x = "Weight [g]",
       y = "Number of fish") +
  scale_y_log10()
p2


### exercise #1
log_breaks2 <- seq(from = 0, to = 11, by = .5)
log_breaks2

bin_breaks2 <- 2 ^ log_breaks2
bin_breaks2

p3 <- ggplot(size_data) +
  geom_histogram(aes(weight), fill = "blue", colour = "black",
                 breaks = bin_breaks2) +
  labs(x = "Weight [g]",
       y = "Number of fish") +
  scale_y_log10()
p3
## end exercise #1


p2 + scale_x_log10()


## Density in size
## Binning

data_with_bins <- size_data |>
  mutate(bin = cut(weight, breaks = bin_breaks, right = FALSE,
                   labels = FALSE))
head(data_with_bins)


## Number density

binned_numbers <- data_with_bins |> 
  group_by(bin) |> 
  summarise(Number = n())
binned_numbers

table(data_with_bins$bin)

binned_numbers <- binned_numbers |> 
  mutate(bin_start = bin_breaks[-length(bin_breaks)],
         bin_end = bin_breaks[-1],
         bin_width = bin_end - bin_start,
         Number_density = Number / bin_width)
binned_numbers



binned_numbers <- binned_numbers |>
  mutate(bin_midpoint = exp((log(bin_start) + log(bin_end)) / 2))

p_number_density <- ggplot(binned_numbers) +
  geom_line(aes(x = bin_midpoint, y = Number_density)) +
  labs(x = "Weight [g]", y = "Number density [1/g]")
p_number_density

p_number_density + 
  scale_x_log10() + 
  scale_y_log10()



### exercise #2
data_with_bins2 <- size_data |>
  mutate(bin = cut(weight, breaks = bin_breaks2, right = FALSE,
                   labels = FALSE))
head(data_with_bins2)

binned_numbers2 <- data_with_bins2 |> 
  group_by(bin) |> 
  summarise(Number = n())
binned_numbers2

binned_numbers2 <- binned_numbers2 |> 
  mutate(bin_start = bin_breaks2[-length(bin_breaks2)],
         bin_end = bin_breaks2[-1],
         bin_width = bin_end - bin_start,
         Number_density = Number / bin_width)
binned_numbers2



binned_numbers2 <- binned_numbers2 |>
  mutate(bin_midpoint = exp((log(bin_start) + log(bin_end)) / 2))

p_number_density <- ggplot(binned_numbers2) +
  geom_line(aes(x = bin_midpoint, y = Number_density)) +
  labs(x = "Weight [g]", y = "Number density [1/g]") + 
  scale_x_log10() + 
  scale_y_log10()
p_number_density
## end exercise #2


## Fitting a power law

model <- lm(log(Number_density) ~ log(bin_midpoint), data = binned_numbers)
model
summary(model)


ggplot(binned_numbers, aes(x = bin_midpoint, y = Number_density)) +
  geom_line() + 
  scale_x_log10(name = "Weight [g]") + 
  scale_y_log10(name = "Number density [1/g]") +
  geom_smooth(method = 'lm')


n <- nrow(size_data)
w_min <- min(size_data$weight)
lambda <- 1 + n / sum(log(size_data$weight / w_min))
lambda



## Biomass density

binned_biomass <- data_with_bins |> 
  group_by(bin) |> 
  summarise(Biomass = sum(weight)) |>
  mutate(bin_start = bin_breaks[-length(bin_breaks)],
         bin_end = bin_breaks[-1],
         bin_width = bin_end - bin_start,
         bin_midpoint = exp((log(bin_start) + log(bin_end)) / 2),
         Biomass_density = Biomass / bin_width)

ggplot(binned_biomass, aes(x = bin_midpoint, y = Biomass_density)) +
  geom_line() + 
  scale_x_log10(name = "Weight [g]") + 
  scale_y_log10(name = "Biomass density") +
  geom_smooth(method = 'lm')


lm(log(Biomass_density) ~ log(bin_midpoint), data = binned_biomass)


binned_data <- 
  left_join(binned_numbers, binned_biomass,
            by = c("bin", "bin_start", "bin_end", "bin_width", 
                   "bin_midpoint")) |>
  mutate(bin_width_log_w = log10(bin_end) - log10(bin_start),
         Number_density_log_w = Number / bin_width_log_w,
         Biomass_density_log_w = Biomass / bin_width_log_w)


long_data <- binned_data |>
  pivot_longer(cols = contains("density"),
               names_to = "Type", values_to = "Density")

ggplot(long_data, aes(x = bin_midpoint, y = Density, colour = Type)) +
  geom_line() + 
  scale_x_log10(name = "Weight [g]") + 
  scale_y_log10(name = "Density") +
  geom_smooth(method = 'lm', se = F)


### Excercise #3
names(binned_data)

lm(log(Biomass_density_log_w) ~ log(bin_midpoint), data = binned_data)
### end exercise #3


## Sheldon's observation
# Sheldon's density = biomass density in long weight
## Size spectra of individual species

p <- ggplot(size_data) +
  geom_density(aes(weight, stat(count), colour = species), adjust = 4) +
  geom_density(aes(weight, stat(count)), colour = "black", lwd = 1.2, adjust = 4) +
  scale_x_continuous(trans = "log10", name = "Weight [g]") +
  scale_y_continuous(trans = "log10", limits = c(1, NA), name = "Number density in log w")

plotly::ggplotly(p)


