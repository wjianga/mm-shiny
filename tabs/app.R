#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    [ http://shiny.rstudio.com/ ]
#

library(shiny)
library(arrow)
library(DT)
library(data.table)
library(maftools)
# library(rsconnect)

cols = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification",
         "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Refseq_mRNA_Id", 
         "Tumor_Sample_Barcode", "dbSNP_RS", "Genome_Change", "cDNA_Change", "Codon_Change",
         "Protein_Change", "Description", "CGC_Name", "CGC_Cancer_Somatic_Mut", "CGC_Cancer_Germline_Mut")

# a data.frame that contains all the expressed mutations
exp_mut = data.frame(fread("/Users/alan/Desktop/MM_Shiny/data/exp_muts/all-samples.filtered.expressed-mutations.tsv", 
                           sep = "\t", header = TRUE, quote = "", drop = "Number"))

# filter the expressed mutations and keep the essential ones
exp_mut = exp_mut[, names(exp_mut) %in% cols]

# a maf object that contains all the expressed mutations
exp_maf <- read.maf(exp_mut)

#----------------------------------------

# a data.frame object that contains SNVs
gatk_df = read_feather("/Users/alan/Desktop/MM_Shiny/data/SNVs/gatk-funcotator-snakemake.feather")

# a maf object that contains all the SNVs
gatk_maf = read.maf(gatk_df)
#gatk_maf = load("./SNV.RData")

# filter the SNVs and keep essential ones
gatk_df = subset(gatk_df, select = names(gatk_df) %in% cols)

#----------------------------------------

# a data.frame object that contains SVs Union
#sv_union_df = data.frame(fread("./Data/SVs/SVsSharedByMoreThan3Tools-margin-100bps.csv"))

# a data.frame object that contains SVs Intersection
sv_table_using = data.frame(fread("/Users/alan/Desktop/MM_Shiny/data/SVs/SVsSharedBy3Tools-margin:100bps.csv"))

###---------------------------------------------------------------------------------------
# Define UI for application that draws a histogram
ui <- fluidPage(
    
        navbarPage(title = "70 Multiple Myeloma Cell Lines",
                   tabPanel(title = "Expressed Mutations",
                            
                            sidebarLayout(
                                sidebarPanel(
                                    selectInput(inputId = "cellLines",
                                                label = "Cell Lines:",
                                                choices = c("Click to see options", 
                                                            unique(exp_mut$Tumor_Sample_Barcode)),
                                                selected = "Click to see options"),
                                    
                                    textInput(inputId = "multipleCellLines",
                                              label = "Select Multiple Cell Lines
                                              (Please use space to separate. e.g. 'almc1 xg7 ...')"),
                                    
                                    selectInput(inputId = "chromosome",
                                                label = "Chromosome",
                                                choices = c("Click to see options", 
                                                            unique(exp_mut$Chromosome)),
                                                selected = "Click to see options"),
        
                                    textInput(inputId = "multipleChromosome",
                                              label = "Select Multiple Chromosomes
                                              (Please use space to separate. e.g. 'chr1 chr2 ...')"),
                                    
                                    selectInput(inputId = "hugo_gene",
                                                label = "Genes",
                                                choices = c("Click to see options", 
                                                            unique(exp_mut$Hugo_Symbol)),
                                                selected = "Click to see options"),
                                    
                                    textInput(inputId = "multipleGenes",
                                              label = "Select Multiple Genes
                                              (Please use space to separate. e.g. 'sfpq msh2 ...')"),
                                    
                                    textInput(inputId = "lolliSample",
                                              label = "Enter gene to visualize in Lollipop Plot
                                              (One gene at a time please)")
                                    ),
                                mainPanel(
                                    DTOutput(outputId = "table"),
                                    plotOutput(outputId = "exp_mut_summary"),
                                    plotOutput(outputId = "exp_oncoplot"),
                                    plotOutput(outputId = "exp_TiTv"),
                                    plotOutput(outputId = "exp_lolipop")
                                )
                            )
                    ),
                   
                   #----------------------------------------
                   
                   tabPanel(title = "Single Nucleotides Variants (SNV)",
                            
                            sidebarLayout(
                                sidebarPanel(
                                    selectInput(inputId = "gatk_cellLines",
                                                label = "Cell Lines:",
                                                choices = c("Click to see options", 
                                                            unlist(strsplit(toString(unique(gatk_df[,Tumor_Sample_Barcode])), 
                                                            split = ", "))),
                                                selected = "Click to see options"),

                                    textInput(inputId = "gatk_multipleCellLines",
                                              label = "Select Multiple Cell Lines
                                              (Please use space to separate. e.g. 'almc1 xg7 ...')"),

                                    selectInput(inputId = "gatk_chromosome",
                                                label = "Chromosome",
                                                choices = c("Click to see options", unique(gatk_df$Chromosome)),
                                                selected = "Click to see options"),

                                    textInput(inputId = "gatk_multipleChromosome",
                                              label = "Select Multiple Chromosomes
                                              (Please use space to separate. e.g. 'chr1 chr2 ...')"),

                                    textInput(inputId = "gatk_multipleGenes",
                                              label = "Select Multiple Genes
                                              (Please use space to separate. e.g. 'sfpq msh2 ...')"),
                                    
                                    textInput(inputId = "sample_to_rain",
                                              label = "Enter sample to visualize in Rainfall Plot
                                              (One sample at a time please)"
                                    )
                                ),

                                mainPanel(
                                    DTOutput(outputId = "SNV_table"),
                                    plotOutput(outputId = "rainfallPlot"),
                                    plotOutput(outputId = "snv_mut_summary"),
                                    plotOutput(outputId = "snv_oncoplot"),
                                    plotOutput(outputId = "snv_TiTv")
                                )
                            )
                            
                    ),
                   
                   #----------------------------------------
                   
                   tabPanel(title = "Structural Variants (SV)",
                            sidebarLayout(
                                sidebarPanel(
                                    # selectInput(inputId = "sv_table_mode",
                                    #             label = "Please choose 'union' or 'intersection' across 3 tools:",
                                    #             choices = c("Union", "Intersection"),
                                    #             selected = "Union"),
                                    
                                    selectInput(inputId = "sv_cellLines",
                                                label = "Cell Lines:",
                                                choices = c("Click to see options", 
                                                            unlist(strsplit(toString(unique(sv_table_using[,"CellLineName"])), split = ", "))),
                                                selected = "Click to see options"),
                                    
                                    textInput(inputId = "sv_multipleCellLines",
                                              label = "Select Multiple Cell Lines
                                              (Please use space to separate. e.g. 'almc1 xg7 ...')"),
                                    
                                    selectInput(inputId = "sv_chr1",
                                                label = "Chr1 options:",
                                                choices = c("Click to see options",
                                                            unlist(strsplit(toString(unique(sv_table_using[,"CHROM1"])), split = ", "))),
                                                selected = "Click to see options"),
                                    
                                    textInput(inputId = "sv_multipleChr1",
                                              label = "Select Multiple Chr1
                                              (Please use space to separate. e.g. 'chr1 ch2 ...'"),
                                    
                                    selectInput(inputId = "sv_chr2",
                                                label = "Chr2 options:",
                                                choices = c("Click to see options",
                                                            unlist(strsplit(toString(unique(sv_table_using[,"CHROM2"])), split = ", "))),
                                                selected = "Click to see options"),
                                    
                                    textInput(inputId = "sv_multipleChr2",
                                              label = "Select Multiple Chr2
                                              (Please use space to separate. e.g. 'chr1 ch2 ...'")
                                ),
                                
                                mainPanel(
                                    DTOutput(outputId = "SV_table"),
                                    htmlOutput(outputId = "sv_html")
                                )
                            )
                            
                   ),
                   
                   #----------------------------------------
                   
                   tabPanel(title = "Copy Number Variants (CNV)"
                            
                            
                   )
        )
)

###---------------------------------------------------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output) {
    # Reactive function that filter expressed mutation dataset
    df_subset = reactive({
        if (input$multipleCellLines != ""){
            exp_mut = exp_mut[tolower(exp_mut$Tumor_Sample_Barcode) %in% unlist(strsplit(tolower(input$multipleCellLines), " ")), ]
        } else {
            exp_mut = exp_mut
        }

        if (input$multipleChromosome != ""){
            exp_mut = exp_mut[tolower(exp_mut$Chromosome) %in% unlist(strsplit(tolower(input$multipleChromosome), split = " ")), ]
        } else {
            exp_mut = exp_mut
        }
        
        if (input$multipleGenes != ""){
            exp_mut = exp_mut[tolower(exp_mut$Hugo_Symbol) %in% unlist(strsplit(tolower(input$multipleGenes), split = " ")), ]
        } else {
            exp_mut = exp_mut
        }
        
        return(exp_mut)
        })
    
    output$table = renderDT(
        expr = df_subset(),
        options = list(scrollX = TRUE)
        )
    
    output$exp_mut_summary = renderPlot(
        plotmafSummary(exp_maf)
    )
    
    output$exp_oncoplot = renderPlot(
        oncoplot(exp_maf)
    )
    
    exp.titv = titv(maf = exp_maf, plot = FALSE, useSyn = TRUE)
    
    output$exp_TiTv = renderPlot(
        plotTiTv(res = exp.titv)
    )
    
    output$exp_lolipop = renderPlot(
        lollipopPlot(exp_maf, gene = toupper(input$lolliSample), AACol = "Protein_Change", showMutationRate = T, tsb = "INA6")
    )
    
    #lollipopPlot(subsetMaf(exp_maf, tsb = c("INA6", "ALMC1")), gene = "SET")
   
#---------------------------------
    
    ## Output for "SNVs"
    # gatk_maf_subset = read.maf(gatk_df)

    gatk_snv_subset = reactive({
        if (input$gatk_multipleCellLines != "" | input$gatk_multipleChromosome != "" | input$gatk_multipleGenes != ""){
            if (input$gatk_multipleCellLines != ""){
                gatk_df = gatk_df[tolower(gatk_df$Tumor_Sample_Barcode) %in% unlist(strsplit(tolower(input$gatk_multipleCellLines), " ")), ]
            }

            if (input$gatk_multipleChromosome != ""){
                gatk_df = gatk_df[tolower(gatk_df$Chromosome) %in% unlist(strsplit(tolower(input$gatk_multipleChromosome), split = " ")), ]
            }

            if (input$gatk_multipleGenes != ""){
                gatk_df = gatk_df[tolower(gatk_df$Hugo_Symbol) %in% unlist(strsplit(tolower(input$gatk_multipleGenes), split = " ")), ]
            }
        } else {
            gatk_df = NULL
        }
        
        # gatk_maf_subset = read.maf(gatk_df)

        return(gatk_df)
    })
    
    # gatk_maf_subset = read.maf(gatk_snv_subset())
    
    # if (! is.null(gatk_df)) {
    #     gatk_maf_subset = read.maf(gatk_maf_subset)
    #     
    #     output$snv_mut_summary = renderPlot(
    #         plotmafSummary(gatk_maf_subset)
    #     )
    #     
    #     output$snv_oncoplot = renderPlot(
    #         oncoplot(gatk_maf_subset)
    #     )
    #     
    #     snv.titv = titv(maf = gatk_maf_subset, plot = FALSE, useSyn = TRUE)
    #     
    #     output$snv_TiTv = renderPlot(
    #         plotTiTv(res = snv.titv)
    #     )
    # }
    
    output$rainfallPlot = renderPlot(
        rainfallPlot(gatk_maf, tsb = toupper(toString(input$sample_to_rain)))
    )
    
    output$SNV_table = renderDT(
        expr = gatk_snv_subset(),
        options = list(scrollX = TRUE)
    )
    

    # output$snv_mut_summary = renderPlot(
    #     plotmafSummary(gatk_maf_subset)
    # )

    # output$snv_oncoplot = renderPlot(
    #     oncoplot(gatk_maf_subset)
    # )
    # 
    # snv.titv = titv(maf = gatk_maf_subset, plot = FALSE, useSyn = TRUE)
    # 
    # output$snv_TiTv = renderPlot(
    #     plotTiTv(res = snv.titv)
    # )
    
#-----------------------------------
    
    ## Output for "SVs"
    # Reactive function that filter expressed mutation dataset
    sv_subset = reactive({
        # if (input$sv_table_mode == "Union"){
        #     sv_table_using = sv_union_df
        # } else {
        #     sv_table_using = sv_intersection_df
        # }
        
        if (input$sv_multipleCellLines != ""){
            sv_table_using = sv_table_using[tolower(sv_table_using$CellLineName) %in% unlist(strsplit(tolower(input$sv_multipleCellLines), " ")), ]
        } else {
            sv_table_using = sv_table_using
        }
        
        if (input$sv_multipleChr1 != ""){
            sv_table_using = sv_table_using[tolower(sv_table_using$CHROM1) %in% unlist(strsplit(tolower(input$sv_multipleChr1), split = " ")), ]
        } else {
            sv_table_using = sv_table_using
        }

        if (input$sv_multipleChr2 != ""){
            sv_table_using = sv_table_using[tolower(sv_table_using$CHROM2) %in% unlist(strsplit(tolower(input$sv_multipleChr2), split = " ")), ]
        } else {
            sv_table_using = sv_table_using
        }
        
        return(sv_table_using)
    })
    
    output$SV_table = renderDT(
        expr = sv_subset(),
        options = list(scrollX = TRUE)
    )
    
    # getPage = function(){
    #     return(includeHTML('./Data/SVs/XG7_PLB.html'))
    # }
    # 
    # output$sv_html = renderUI({getPage()})
}

# Run the application 
shinyApp(ui = ui, server = server)
