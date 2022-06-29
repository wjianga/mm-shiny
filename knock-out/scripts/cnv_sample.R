library(tidyverse)
library(colorspace)
library(glue)

cnv = read.csv("./data/CNV/cnv.csv")

plot_cnv <- function(sample) {
  
  cell_line <- cnv %>% 
    filter(cell_line == sample)
  
  reorder = unique(cell_line$chromosome)

  organised_data <- cell_line %>% 
    mutate(chromosome = fct_relevel(chromosome, reorder)) %>% 
    group_by(type, chromosome) %>% 
    summarise(
      count = n()
    )

  ggplot(organised_data) +
    aes(x= ifelse( type =="duplication",count, -count), y = chromosome, fill = type) +
    geom_col(position = "stack")+
    scale_x_continuous(
      limits = c(-180, 180), 
      breaks = seq(-180, 180, 20), 
      labels = abs(seq(-180, 180, 20)), 
      name = "Variation Count") +
    scale_fill_discrete(
      name = "Variation") +
    scale_y_discrete(
      name = "Chromosome") +
    theme_bw() +
    labs(title = sample)
}
  

