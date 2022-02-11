library(shiny)

library(arrow)
library(DT)
library(reactable)
library(data.table)

library(maftools)

# library(rsconnect) # connect to shiny.io

cols = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Refseq_mRNA_Id", "Tumor_Sample_Barcode", "dbSNP_RS", "Genome_Change", "cDNA_Change", "Codon_Change",
         "Protein_Change", "Description", "CGC_Name", "CGC_Cancer_Somatic_Mut", "CGC_Cancer_Germline_Mut")

#---------------------------------------- SNV

# a data.frame object that contains SNVs
gatk_df = read_feather("/Users/alan/Desktop/MM_Shiny/data/SNVs/gatk-funcotator-snakemake.feather")

# filter the SNVs and keep essential ones
gatk_df = subset(gatk_df, select = names(gatk_df) %in% cols)

# a maf object that contains all the SNVs
gatk_maf = read.maf(gatk_df)


#---------------------------------------- Expressed Mutation

# a data.frame that contains all the expressed mutations
exp_mut = data.frame(fread("/Users/alan/Desktop/MM_Shiny/data/exp_muts/all-samples.filtered.expressed-mutations.tsv", sep = "\t", header = TRUE, quote = "", drop = "Number"))

# filter the expressed mutations and keep the essential ones
exp_mut = exp_mut[, names(exp_mut) %in% cols]

# a maf object that contains all the expressed mutations
exp_maf <- read.maf(exp_mut)


################################################################################## UI

ui <- fluidPage(
    
  # Titles and Logo
  titlePanel(
    div(img(src = "logo.png", height = 90), 
        "Multiple Myeloma")
  ),
  
  # Create a tab bars
  tabsetPanel(
    id = "tabs",
    
############################################## Landing page  
    tabPanel(
      title = "Welcome",
      column(6, textAreaInput(
        inputId = "genes",
        label = h4("Enter genes you want to knock off"),
        width = "90%",
        rows = 10,
        placeholder = "sfpq msh2 fus msn ... (Please use space to separate individual Hugo_symbol)"
      ),
      
      actionButton(
        inputId = "goGene",
        label = "Go!",
        class = "btn-primary"
      )
  
  ),
  
  column(6, textAreaInput(
    inputId = "samples",
    label = h4("Enter samples you want to see"),
    width = "90%",
    rows = 10,
    placeholder = "almc1 xg7 amo1 ... (Please use space to separate individual cell line)"
  )
),

  ),

############################################## Gene Information 
  tabPanel(
    title = "Mutation Burden",
    
    dataTableOutput(
       outputId = "mb_diff"
    ),
    
    textOutput(
      outputId = "test"
    )
    
  )
  )

)

################################################################################## Server

server <- function(input, output, session) {
  
  genes <- reactive(unlist(strsplit(toupper(input$genes), " ")))
  
  samples <- reactive(unlist(strsplit(toupper(input$samples), " ")))

  
############################################## Handle the button "Go" on the landing page
  observeEvent(
    input$goGene, {
      updateTabsetPanel(
        inputId = "tabs",
        selected = "Mutation Burden" # page to switch to
      )
    }
  )

############################################## Display the difference mutation burden table
  
  exp_mb <- reactive(tmb(subsetMaf(maf = exp_maf, 
                                   genes = setdiff(exp_maf@data[["Hugo_Symbol"]], genes())))[, c(1, 3)])
  
  snv_mv <- reactive(tmb(subsetMaf(maf = gatk_maf, 
                                   genes = setdiff(gatk_maf@data[["Hugo_Symbol"]], genes())))[, c(1, 3)])

  output$mb_diff <- renderDataTable(
    datatable(
      merge(exp_mb(), snv_mv(), by = "Tumor_Sample_Barcode", all = T),
      colnames = c("Tumor_Sample_Barcode", "Expressed_Mutation", "SNV")
    )
  )
  
  output$test <- renderText(
    print(genes())
  )


}

shinyApp(ui, server)
