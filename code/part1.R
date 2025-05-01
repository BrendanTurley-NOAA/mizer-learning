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
