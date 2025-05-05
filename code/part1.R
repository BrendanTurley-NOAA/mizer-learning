### Part 1

install.packages(c("mizer", "tidyverse", "plotly", "remotes", "usethis",
                   "rmarkdown", "rstudioapi"))
library(tidyverse)

download.file("https://github.com/gustavdelius/mizerCourseNew/raw/master/understand/size-data.rds",
              destfile = "data/size-data.rds")

size_data <- readRDS("data/size-data.rds")

str(size_data)
unique(size_data$species)

hist(size_data$weight, xlab = 'weight [g]', ylab = 'number of fish', las = 1)
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



data_with_bins <- size_data |>
  mutate(bin = cut(weight, breaks = bin_breaks, right = FALSE,
                   labels = FALSE))
head(data_with_bins)


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

