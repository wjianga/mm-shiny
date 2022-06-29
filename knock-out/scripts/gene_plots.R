library(tidyverse)

#load in neoantigen prediction table
neoantigen <- read.csv("./data/neoantigen/neoantigen.csv")


#select genes to plot
gene_list <- c("CREB3L2","TCF3","CTNNB1","XPO1","ARNT","GNAS","KRAS",
               "SMARCA4","TPR","NRAS","NUP98","BIRC3","DDX6","MYH9","FGFR3","JAK1",
               "MLH1","PRDM1","RB1","IDH2","BRAF","CDKN2C")


#gene plot function
gene_plot <- function(Gene_name){
  neoantigen %>% 
    filter(.data$Gene.Name %in% .env$Gene_name) %>% 
    filter(.data$Corresponding.WT.Score > 500) %>% 
    group_by(.data$Gene.Name) %>% 
    summarise(
      count = n()
    ) %>%
    ggplot()+
    aes(y = fct_reorder(.data$Gene.Name, count), x = count)+
    geom_col(fill = "steelblue")+
    scale_y_discrete(
      name = "Gene"
    )+
    scale_x_continuous(
      name = "Neoantigen Load", 
      expand = expansion(mult = c(0,0.06))
    )+
    theme( 
      axis.text = element_text(
        colour ="black",
        size = 13
      ),
      axis.title = element_text(
        size = 15
      ),
      axis.title.x = element_text(
        vjust = 0.1
      ),
      panel.grid.major.y = element_blank()
    )
  
  
  
}
#call function
gene_plot(gene_list)
