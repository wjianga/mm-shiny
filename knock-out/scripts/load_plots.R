library(tidyverse)

#load in neoantigen prediction table
neoantigen <- read.csv("./data/neoantigen/neoantigen.csv")


#select genes to plot
gene_list <- c("CREB3L2","TCF3","CTNNB1","XPO1","ARNT","GNAS","KRAS",
               "SMARCA4","TPR","NRAS","NUP98","BIRC3","DDX6","MYH9","FGFR3","JAK1",
               "MLH1","PRDM1","RB1","IDH2","BRAF","CDKN2C")

#select cell lines
cell_line = unique(neoantigen$cell_line)



`%!in%` <- Negate(`%in%`) #create the not operator for knocking genes

#function to plot neoantigen load for selected cell lines and knocking genes
neoantigen_load <- function(cell_line, knock_genes){
  neoantigen %>% 
    filter(Corresponding.WT.Score > 500) %>%  
    filter(.data$cell_line %in% .env$cell_line, .data$Gene.Name %!in% .env$knock_genes) %>% 
    group_by(.data$cell_line) %>% 
    summarise(
      count = n()
    ) %>%
    ggplot()+
    aes(y = fct_reorder(.data$cell_line, count), x = count)+
    geom_col(fill = "steelblue")+
    scale_y_discrete(
      name = "Cell Line"
    )+
    scale_x_continuous(
      name = "Neoantigen Load", 
      expand = expansion(mult = c(0,0.06))
    )+
    theme( 
      axis.text = element_text(
        colour ="black",
        size = 8
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
neoantigen_load(cell_line, gene_list)

