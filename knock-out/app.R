library(shiny)

# library(rsconnect) # connect to shiny.io

# library(maftools) # use when publishing
# library(arrow) # use when publishing

library(arrow, lib.loc = "/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/packages") # read in feather file
library(maftools, lib.loc = "/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/packages") # snv manipulation

library(DT) # data table display
library(tidyverse)
library(plotly)
library(data.table)

source("./scripts/circosplot.R") # script to plot circos plot
source("./scripts/load_plots.R") # script to plot the neoantigen load
source("./scripts/gene_plots.R") # plot the driver genes
source("./scripts/cnv_sample.R") # plot cnv sample
source("./scripts/calculateFGA.R") # calculate FGA

# source("/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/scripts/circosplot.R") # script to plot circos plot
# source("/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/scripts/load_plots.R") # script to plot the neoantigen load
# source("/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/scripts/gene_plots.R") # plot the driver genes
# source("/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/scripts/cnv_sample.R") # plot cnv sample
# source("/stor/work/FRI-BigDataBio/wj3972/mm_shiny/knock-out/scripts/calculateFGA.R") # calculate FGA

cols = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand", 
         "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", 
         "Tumor_Seq_Allele2", "Refseq_mRNA_Id", "Tumor_Sample_Barcode", "dbSNP_RS", 
         "Genome_Change", "cDNA_Change", "Codon_Change", "Protein_Change", "Description", 
         "CGC_Name", "CGC_Cancer_Somatic_Mut", "CGC_Cancer_Germline_Mut")

`%!in%` = Negate(`%in%`)

#---------------------------------------- SNV - read in file

# a data.frame object that contains SNVs
gatk_df = read_feather("./data/SNVs/snakemake-67samp-March22-noFiltered.feather") %>% 
    # keep only columns in <cols>
    select(all_of(cols))

gatk_maf = read.maf(gatk_df) # a maf object that contains all the SNVs


#---------------------------------------- Expressed Mutation - read in file

# a data.frame that contains all the expressed mutations
exp_mut = data.frame(fread("./data/exp_muts/mutect2-March22-67samp-expressed-mutations.tsv", 
                           sep = "\t", header = TRUE, quote = "", drop = "Number")) %>% 
            select(all_of(cols)) # keep only columns in <cols>

exp_maf <- read.maf(exp_mut) # a maf object that contains all the expressed mutations


#---------------------------------------- SV - read in file

sv = read.csv("./data/SVs/SVsSharedBy3Tools-margin:100bps-withGeneColumn-filteredForCircosPlots.csv")


#---------------------------------------- CNV - read in file

cnv = read.csv("./data/CNV/cnv.csv")

#---------------------------------------- neoantigen - read in file

neo_csv = read.csv("./data/neoantigen/neoantigen.csv")

# neo_csv <- neo_csv %>% 
#     select(Best.MT.Score, Best.MT.Score.Method, Tumor.DNA.Depth, Tumor.DNA.VAF, 
#            Tumor.RNA.Depth, Tumor.RNA.VAF, Normal.Depth, Normal.VAF, Gene.Expression,
#            Transcript.Expression, Corresponding.WT.Score, Gene.Name)


############################################################################################################################# UI

ui <- fluidPage(
    
    # Titles and Logo
    titlePanel(
        div(img(src = "logo.png", height = '90'), # FRI logo
            "Characterization of 70 Multiple Myeloma Cell Lines", # title
            img(src = "livestrong.png", height = '95', align = "right")) # livestrong logo
    ),
    
    # multiple tabs starts here (e.g. Welcome, SNV, ...)
    # the backbone of UI
    tabsetPanel(
        id = "tabs",
        
        ########################################################################################################## Landing page - UI
        tabPanel(
            title = "Welcome", # tab name
            
            fluidRow(
                # first text box that takes in genes
                column(6, textAreaInput(
                    inputId = "genes",
                    label = h4("Enter genes you want to investigate"),
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
                
                # second text box that takes in samples used to filter table
                column(6, textAreaInput(
                    inputId = "samples",
                    label = h4("Enter cell lines you want to visualize"),
                    width = "90%",
                    rows = 10,
                    placeholder = "almc1 xg7 amo1 ... (Please use space to separate individual cell line) \nWhen empty, the entire dataset will be used."
                )
                )
            ),
            
            # place where instruction on the first page is displayed
            column(6, 
                htmlOutput(outputId = "introduction"),
                # align text to left
                align = "left",
                # move the text block a bit to the right
                offset = 3
            )
            
        ),
        
        ########################################################################################################## Mutation Burden - UI
        tabPanel(
            title = "Mutation Burden",
            
            # display the user entered genes list
            titlePanel(
                h3(textOutput(
                    outputId = "gene_entered_burden_tab"
                ))),
            
            # download tab
            downloadButton(outputId = "mutation_burden_download", 
                           label = "Download CSV"),
            
            # table of mutation burden
            dataTableOutput(
                outputId = "mb_diff"
            ),
            
            uiOutput(
                outputId = "proportion_stacked_barplot"
            )
            
        ),
        
        
        ########################################################################################################## SNV - UI
        
        tabPanel(
            title = "Single Nucleotide Variants",
            
            navlistPanel(
                
                tabPanel(
                    title = "Overview Visuals",
                    plotOutput(outputId = "snv_summary"),
                    plotOutput(outputId = "snv_oncoplot")
                ),
                
                tabPanel(
                    title = "Selected Samples Visuals",
                    h3(textOutput(
                        outputId = "sample_entered_snv_tab"
                    )),
                    uiOutput(outputId = "rainfall_snv")
                ),
                
                tabPanel(
                    title = "SNV Table",
                    
                    downloadButton(outputId = "snv_download", 
                                   label = "Download CSV"),
                    
                    dataTableOutput(outputId = "snv_table")
                ),
                
                # ratio of (navbar : main panel)
                widths = c(3, 9)
            )
        ),
        
        ########################################################################################################## Expressed Mutation - UI
        
        tabPanel(
            title = "Expressed Mutation",
            
            navlistPanel(
                
                tabPanel(
                    title = "Overview Visuals",
                    
                    plotOutput(outputId = "exp_summary"),
                    plotOutput(outputId = "exp_oncoplot")
                ),
                
                tabPanel(
                    title = "Selected Samples Visuals",
                    h3(textOutput(
                        outputId = "sample_entered_em_tab"
                    )),
                    uiOutput("rainfall_em")
                ),
                
                tabPanel(
                    title = "Expressed Mutation Table",
                    
                    downloadButton(outputId = "em_download",
                                   label = "Download CSV"),
                    
                    dataTableOutput(outputId = "em_table")
                ),
                
                # ratio of (navbar : main panel)
                widths = c(3, 9)
            )
            
        ),
        
        ########################################################################################################## SV - UI
        
        tabPanel(
            title = "Structural Variants",
            
            navlistPanel(
                
                tabPanel(
                    title = "Overview Visuals",
                    htmlOutput(outputId = "overview_circos")
                ),
                
                tabPanel(
                    title = "Selected Samples Visuals",
                    h3(textOutput(
                        outputId = "sample_entered_sv_tab"
                    )),
                    htmlOutput(outputId = "selected_circos")
                ),
                
                tabPanel(
                    title = "SV Table",
                    
                    downloadButton(outputId = "sv_download",
                                   label = "Download CSV"),
                    
                    dataTableOutput(outputId = "sv_table")
                ),
                
                # ratio of (navbar : main panel)
                widths = c(3, 9)
            )
        ),
        
        ########################################################################################################## CNV - UI
        
        tabPanel(
            title = "Copy Number Variants",
            
            navlistPanel(
                
                tabPanel(
                    title = "Overview Visuals",
                    plotOutput(outputId = "cnv_overview_plot")
                ),
                
                tabPanel(
                    title = "Selected Samples Visuals",
                    h3(textOutput(
                        outputId = "sample_entered_cnv_tab"
                    )),
                    uiOutput(outputId = "cnv_sample_plot")
                ),
                
                tabPanel(
                    title = "CNV Table",
                    
                    downloadButton(outputId = "cnv_download",
                                   label = "Download CSV"),
                    
                    dataTableOutput(outputId = "cnv_table")
                ),
                
                # ratio of (navbar : main panel)
                widths = c(3, 9)
            )
        ),
        
        
        ########################################################################################################## Neoantigen - UI
        
        tabPanel(
            title = "Neoantigen",
            
            navlistPanel(
                
                tabPanel(
                    title = "Overview Visuals",
                    plotOutput(outputId = "neo_count"),
                    # plotOutput(outputId = "neo_burden")
                    h4(textOutput(
                        outputId = "driver_genes_text"
                    )),
                    plotOutput(outputId = "driver_plot")
                ),
                
                tabPanel(
                    title = "Selected Samples Visuals",
                    h3(textOutput(
                        outputId = "sample_entered_neo_tab"
                    )),
                    plotOutput(outputId = "neo_selected")
                ),
                
                tabPanel(
                    title = "Neoantigen Table",
                    
                    downloadButton(outputId = "neo_download",
                                   label = "Download CSV"),
                    
                    dataTableOutput(outputId = "neo_table")
                ),
                
                # ratio of (navbar : main panel)
                widths = c(3, 9)
            )
            
        )
    )
    
)

########################################################################################################## Server

server <- function(input, output, session) {
    
    # store the user input in "genes" box into a vector
    genes <- reactive(unlist(strsplit(toupper(trimws(input$genes)), " ")))
    
    # store the user input in "samples" box into a vector
    samples <- reactive(unlist(strsplit(toupper(trimws(input$samples)), " ")))
    
    output$introduction = renderText(
        paste("<h4>Instructions to use this dashboard:</h4><br>
        
        The 'Mutation Burden' displays the mutation burden, and the following five tabs each displays information for a particular variant type.<br>
        <br>
        The left text box takes in genes that you want to remove, and mutation burden change will be shown under 'Mutation Burden' tab.<br>
        <br>
        The right text box takes in cell lines that you want to investigate, and visualizations associated with selected cell lines and truncated table are displayed under each variant type's tab. If no cell line is inputted, then the whole dataset will be displayed")
    )
    
    snv_filtered_by_samples <- reactive(
        if (length(samples()) == 0){
            gatk_df
        } else {
            gatk_df %>% filter(Tumor_Sample_Barcode %in% samples())
        }
    )
    
    em_filtered_by_samples <- reactive(
        if (length(samples()) == 0){
            exp_mut
        } else {
            exp_mut %>% filter(Tumor_Sample_Barcode %in% samples())
        }
    )
    
    sv_filtered_by_samples <- reactive(
        if (length(samples()) == 0){
            sv
        } else {
            sv %>% 
                mutate_at("CellLineName", str_remove, "_.*") %>% 
                filter(CellLineName %in% samples())
        }
    )
    
    neo_filtered_by_samples <- reactive(
        if (length(samples()) == 0){
            neo_csv
        } else {
            neo_csv %>% filter(cell_line %in% samples())
        }
    )
    
    cnv_filtered_by_samples <- reactive(
        if (length(samples()) == 0){
            cnv
        } else {
            cnv %>% filter(cell_line %in% samples())
        }
    )
    
    
    ########################################################################################################## Landing Page - Server
    observeEvent(
        input$goGene, {
            updateTabsetPanel(
                inputId = "tabs",
                selected = "Mutation Burden" # page to switch to
            )
        }
    )
    
    ########################################################################################################## Mutation Burden - Server
    
    # generate expressed mutation mutation burden
    # including mutation burden changes and original maf object burden
    exp_mb <- reactive(tmb(subsetMaf(maf = exp_maf, 
                                     genes = setdiff(exp_maf@data[["Hugo_Symbol"]], genes()))) %>% 
                           select(1, 3) %>%
                           rename(Expressed_Mutation_Burden = total_perMB) %>% 
                           full_join(tmb(exp_maf), by = "Tumor_Sample_Barcode") %>%
                           select(1, 2, 4) %>% 
                           mutate(EM_Burden_Change = Expressed_Mutation_Burden - total_perMB) %>% 
                           select(1, 2, 4))
    
    # generate snv burden
    # including mutation burden change and original maf object burden
    snv_mb <- reactive(tmb(subsetMaf(maf = gatk_maf, 
                                     genes = setdiff(gatk_maf@data[["Hugo_Symbol"]], genes()))) %>% 
                           select(1, 3) %>%
                           rename(SNV_Burden = total_perMB) %>% 
                           full_join(tmb(gatk_maf), by = "Tumor_Sample_Barcode") %>%
                           select(1, 2, 4) %>% 
                           mutate(SNV_Burden_Change = SNV_Burden - total_perMB) %>% 
                           select(1, 2, 4))
    
    
    sv_mb <- reactive(
        # sv %>% 
        #   select(1, 6) %>% 
        #   full_join(sv[(!grepl(paste(genes(), collapse = "|"), sv[[6]])),], by = "CellLineName") %>% 
        #   group_by(CellLineName) %>% 
        #   summarize(count = n()) %>%
        #   mutate(SV_Burden = count / 41.6) %>%
        #   mutate_at("CellLineName", str_remove, "_.*") %>%
        #   select(1, 3) %>%
        #   rename("Tumor_Sample_Barcode" = "CellLineName")
        
        sv %>% 
            # count the number of rows for each cell line
            group_by(CellLineName) %>% 
            summarize(count = n()) %>%
            # calculate the mutation burden for each cell line and add a new column
            mutate(SV_Burden = count / 41.6) %>%
            mutate_at("CellLineName", str_remove, "_.*") %>%
            select(1, 3) %>%
            rename("Tumor_Sample_Barcode" = "CellLineName") %>% 
            full_join((sv[(!grepl(paste(genes(), collapse = "|"), sv[[6]])),] %>% 
                           group_by(CellLineName) %>% 
                           summarize(count = n()) %>% 
                           mutate(SV_Burden_Off = count / 41.6) %>%
                           mutate_at("CellLineName", str_remove, "_.*") %>%
                           select(1, 3) %>%
                           rename("Tumor_Sample_Barcode" = "CellLineName")), by = "Tumor_Sample_Barcode") %>% 
            mutate(SV_Burden_Change = SV_Burden_Off - SV_Burden) %>% 
            select(1, 2, 4)
    )
    
    neo_mb <- reactive(
        neo_csv %>% 
            filter(Corresponding.WT.Score > 500) %>% 
            # keep two columns: 1. hugo symbol 2. cell line name
            select(Gene.Name, cell_line) %>% 
            group_by(cell_line) %>% 
            summarize(Neoantigen_Burden = n() / 41.6) %>% 
            full_join(
                (
                    neo_csv[(!grepl(paste(genes(), collapse = "|"), neo_csv[["Gene.Name"]])),] %>% 
                        filter(Corresponding.WT.Score > 500) %>% 
                        select(Gene.Name, cell_line) %>% 
                        group_by(cell_line) %>% 
                        summarize(count_after = n() / 41.6)
                ), by = "cell_line") %>% 
            mutate(Neoantigen_Change = count_after - Neoantigen_Burden) %>% 
            select(c("cell_line", "Neoantigen_Burden", "Neoantigen_Change"))
    )
    
    cnv_mb <- reactive(
        cnv %>% 
            # this line excludes chromosomes Y and M (chromosome X shouldn't be in the dataset)
            filter(chromosome != "chrY" & chromosome != "chrM") %>% 
            group_by(cell_line) %>% 
            summarize(CNV_burden = sum(end - start) / 3099734149) %>% 
            full_join(
                cnv[(!grepl(paste(genes(), collapse = "|"), cnv[["gene"]])), ] %>% 
                    # this line excludes chromosomes Y and M (chromosome X shouldn't be in the dataset)
                    filter(chromosome != "chrY" & chromosome != "chrM") %>% 
                    group_by(cell_line) %>% 
                    summarize(CNV_burden_after = sum(end - start) / 3099734149),
                by = "cell_line"
            ) %>% 
            mutate(CNV_Change = CNV_burden_after - CNV_burden) %>% 
            select(!CNV_burden_after)
    )
    
    # join all the variants into a big table to display
    to_display <- reactive(
        if (length(samples()) == 0){
            snv_mb() %>% 
                full_join(exp_mb(), by = "Tumor_Sample_Barcode") %>% 
                full_join(sv_mb(), by = "Tumor_Sample_Barcode") %>% 
                full_join(cnv_mb(), by = c("Tumor_Sample_Barcode" = "cell_line")) %>% 
                full_join(neo_mb(), by = c("Tumor_Sample_Barcode" = "cell_line"))
        } else {
            snv_mb() %>% 
                full_join(exp_mb(), by = "Tumor_Sample_Barcode") %>% 
                full_join(sv_mb(), by = "Tumor_Sample_Barcode") %>% 
                full_join(cnv_mb(), by = c("Tumor_Sample_Barcode" = "cell_line")) %>% 
                full_join(neo_mb(), by = c("Tumor_Sample_Barcode" = "cell_line")) %>% 
                filter(Tumor_Sample_Barcode %in% samples())
        }
        
    )
    
    output$gene_entered_burden_tab <- renderText(
        paste(c("Genes that have been removed:", genes(), collapse = " "))
    )
    
    output$mutation_burden_download <- downloadHandler(
        filename = "mutation_burden.csv",
        content = function(file) {
            write.csv(to_display(), file)
        }
    )
    
    # render to_display() into a data table
    output$mb_diff <- renderDataTable(
        datatable(
            to_display(),
            rownames = F,
            options = list(scrollX = TRUE)
        )
    )
    
    output$proportion_stacked_barplot <- renderUI(
        ggplotly(to_display() %>%
                     filter(Tumor_Sample_Barcode != "KMS11") %>% 
            select(!matches(".*_change")) %>%
            mutate_at(2:5, function(x)
                (x / rowSums(select_if(., is.numeric), na.rm = T))) %>%
            pivot_longer(cols = !Tumor_Sample_Barcode,
                         names_to = "type",
                         values_to = "burden") %>%
            ggplot(aes(x = Tumor_Sample_Barcode,
                       y = burden,
                       fill = type)) +
            geom_bar(position = "stack", stat = "identity") +
            scale_y_continuous(expand = c(0, 0)) +
            labs(fill = "Variant Type", y = "Proportion of mutation burden explained") +
            theme(axis.text.x = element_text(angle = 45)))
    )
    
    # output$proportion_stacked_barplot <- renderPlot(
    #     to_display() %>%
    #         select(!matches(".*_change")) %>%
    #         mutate_at(2:5, function(x)
    #             (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = T))) %>%
    #         mutate_at(2:5, function(x)
    #             (x / rowSums(select_if(., is.numeric), na.rm = T))) %>%
    #         pivot_longer(cols = !Tumor_Sample_Barcode,
    #                      names_to = "type",
    #                      values_to = "burden") %>%
    #         ggplot(aes(y = Tumor_Sample_Barcode,
    #                    x = burden,
    #                    fill = type)) +
    #         geom_bar(position = "stack", stat = "identity") +
    #         scale_x_continuous(expand = c(0, 0)) +
    #         labs(fill = "Variant Type", x = "Proportion of mutation burden explained (Scaled)")
    # )
    
    ########################################################################################################## SNV tab - Server
    
    output$snv_summary = renderPlot(
        plotmafSummary(gatk_maf)
    )
    
    output$snv_oncoplot = renderPlot(
        oncoplot(gatk_maf)
    )
    
    # output$rainfallPlot_snv = renderPlot(
    #     rainfallPlot(gatk_maf, tsb = toupper(toString(samples())))
    # )
    
    output$sample_entered_snv_tab = renderText(
        paste(c("Cell Lines that have been selected: ", samples(), collapse = " "))
    )
    
    output$rainfall_snv = renderUI({
        plot_output_list_snv <- lapply(1:10, function(i) {
            plotname <- paste("plot_snv", i, sep="")
            plotOutput(plotname)
        })
        
        tagList(plot_output_list_snv)
    })
    
    sapply(1:10, function(i) {
        local({
            plotname <- paste("plot_snv", i, sep="")
            output[[plotname]] <- renderPlot({
                rainfallPlot(gatk_maf, tsb = samples()[i])
            })
        })
    })
    
    output$snv_table = renderDataTable(
        data.table(snv_filtered_by_samples()), options = list(scrollX = TRUE)
    )
    
    output$snv_download = downloadHandler(
        filename = "snv.csv",
        content = function(file) {
            write.csv(snv_filtered_by_samples(), file)
        }
    )
    
    ########################################################################################################## Expressed Mutation tab - Server
    
    output$exp_summary = renderPlot(
        plotmafSummary(exp_maf)
    )
    
    output$exp_oncoplot = renderPlot(
        oncoplot(exp_maf)
    )
    
    # output$rainfallPlot_exp = renderPlot(
    #     rainfallPlot(exp_maf, tsb = toupper(toString(samples())))
    # )
    
    output$sample_entered_em_tab = renderText(
        paste(c("Cell Lines that have been selected: ", samples(), collapse = " "))
    )
    
    output$em_table = renderDataTable(
        em_filtered_by_samples(),
        options = list(scrollX = TRUE)
    )
    
    output$em_download = downloadHandler(
        filename = "expressed_mutation.csv",
        content = function(file) {
            write.csv(em_filtered_by_samples(), file)
        }
    )
    
    output$rainfall_em = renderUI({
        plot_output_list_em <- lapply(1:10, function(i) {
            plotname <- paste("plot_em", i, sep="")
            plotOutput(plotname)
        })
        
        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        # do.call(tagList, plot_output_list)
        tagList(plot_output_list_em)
    })
    
    sapply(1:10, function(i) {
        local({
            plotname <- paste("plot_em", i, sep="")
            
            output[[plotname]] <- renderPlot({
                rainfallPlot(exp_maf, tsb = samples()[i])
            })
        })
    })
    
    ########################################################################################################## SV tab - Server
    
    output$circo_plot = renderImage(
        list(src = "./www/4Tools-0bps.png"),
        deleteFile=FALSE
    )
    
    output$sample_entered_sv_tab = renderText(
        paste(c("Cell Lines that have been selected: ", samples(), collapse = " "))
    )
    
    output$circo_plot_ina6 = renderImage(
        list(src = "./www/ina6.png"),
        deleteFile=FALSE
    )
    
    output$sv_table = renderDataTable(
        sv_filtered_by_samples(),
        options = list(scrollX = TRUE)
    )
    
    output$sv_download = downloadHandler(
        filename = "structural_variants.csv",
        content = function(file) {
            write.csv(sv_filtered_by_samples(), file)
        }
    )
    
    output$overview_circos = renderUI(
        produceCircosPlot(sv)
    )
    
    output$selected_circos = renderUI(
        produceCircosPlot(sv_filtered_by_samples())
    )
    
    ########################################################################################################## CNV - Server
    
    output$cnv_table = renderDataTable(
        cnv_filtered_by_samples(),
        options = list(scrollX = TRUE)
    )
    
    output$cnv_download = downloadHandler(
        filename = "cnv.csv",
        content = function(file) {
            write.csv(cnv_filtered_by_samples(), file)
        }
    )
    
    # output$cnv_overview_plot = renderPlot(
    #     plot_cnv(unique(cnv$cell_line))
    # )
    
    output$sample_entered_cnv_tab = renderText(
        paste(c("Cell Lines that have been selected: ", samples(), collapse = " "))
    )
    
    output$cnv_sample_plot = renderUI({
        plot_output_list_cnv <- lapply(1:10, function(i) {
            plotname <- paste("plot_cnv", i, sep="")
            plotOutput(plotname)
        })
        
        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        # do.call(tagList, plot_output_list)
        tagList(plot_output_list_cnv)
    })
    
    sapply(1:10, function(i) {
        local({
            plotname <- paste("plot_cnv", i, sep="")
            
            output[[plotname]] <- renderPlot({
                plot_cnv(samples()[i])
            })
        })
    })
    
    ########################################################################################################## Neoantigen - Server
    
    output$neo_table = renderDataTable(
        neo_filtered_by_samples(),
        options = list(scrollX = TRUE)
    )
    
    output$neo_download = downloadHandler(
        filename = "neoantigens.csv",
        content = function(file) {
            write.csv(neo_filtered_by_samples(), file)
        }
    )
    
    output$neo_count = renderPlot(
        neo_csv %>%
            filter(Corresponding.WT.Score > 500) %>%
            group_by(cell_line) %>%
            summarize(count = n()) %>%
            ggplot(aes(x = reorder(cell_line, count), y = count)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  panel.background = element_blank()) +
            labs(title = "Number of Neoantigen for 42 cell lines",
                 x = "Cell Line",
                 y = "Count")
        
        # neoantigen_load(neo_csv$cell_line, NULL)
    )
    
    # output$neo_burden = renderPlot(
    #   neo_csv %>% 
    #     group_by(cell_line) %>% 
    #     summarize(burden = n() / 41.6) %>% 
    #     ggplot(aes(x = reorder(cell_line, burden), y = burden)) +
    #     geom_bar(stat = "identity", fill = "seagreen3") +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1),
    #           panel.background = element_blank()) +
    #     labs(title = "Mutation Burden for 42 cell lines",
    #          x = "Cell Line",
    #          y = "Mutation Burden")
    # )
    
    output$sample_entered_neo_tab = renderText(
        paste(c("Cell Lines that have been selected: ", samples(), collapse = " "))
    )
    
    output$neo_selected = renderPlot(
        neoantigen_load(samples(), genes())
    )

    output$driver_genes_text = renderText(
        paste('Here is a list of cancer driver genes searched: "CREB3L2","TCF3","CTNNB1","XPO1","ARNT","GNAS","KRAS",
               "SMARCA4","TPR","NRAS","NUP98","BIRC3","DDX6","MYH9","FGFR3","JAK1",
               "MLH1","PRDM1","RB1","IDH2","BRAF","CDKN2C". \n
              The graph below is what we find')
    )
    
    output$driver_plot = renderPlot(
        gene_plot(gene_list)
    )

}

shinyApp(ui, server)
