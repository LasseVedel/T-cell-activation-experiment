library(shiny)
library(ggplot2)
library(tidyverse)
library(matrixStats) 


server <- function(input, output) {
  
  plotInput <- reactive({
    #Read data file and convert to a matrix 
    normalized_tibble_with_symbols <- read_csv("C:/Users/Lasse/Desktop/Lasse Vedel Jørgensen - Master Project/Bulk RNA-seq/T cell activation experiment/Data/normalized_tibble_with_symbols.csv")
    matsymbol <- as.matrix(normalized_tibble_with_symbols[, 2:21])
    row.names(matsymbol) <-  normalized_tibble_with_symbols$...1
    
    
    ## The gene is found in the matsymbol to extract all normalized read counts 
    ## A matrix of dim=5x4 is formed 
    gene_counts <- t(matrix(matsymbol[input$gene,], nrow=4))
    
    ## Define rownames for the matrix 
    rownames(gene_counts) <- c("Non-stimulated", 
                               "Stimulated, 24h", 
                               "Stimulated, 48h", 
                               "Stimulated, 48h + LV", 
                               "Stimulated, 72h + LV")
    
    ## Calculate rowMeans and rowSDs for each row in the matrix 
    row_means_gene <- rowMeans(gene_counts)
    row_sds_gene <- matrixStats::rowSds(gene_counts)
    
    ## Collect to a dataframe which can be used for ggplot  
    df_gene <- as.data.frame(cbind(row_means_gene, row_sds_gene))
    
    ## Plot the expression using ggplot 
    p_gene <- ggplot(df_gene, aes(x=rownames(df_gene), y=row_means_gene, fill = rownames(df_gene))) + 
      geom_bar(stat="identity", color="grey", position=position_dodge(), width = 0.7) + 
      geom_errorbar(aes(ymin=row_means_gene-row_sds_gene, ymax=row_means_gene+row_sds_gene), width=0.2,
                    position=position_dodge(.9), color = "#404040") + 
      scale_fill_manual("Condition", values = c("Non-stimulated" = "blue", 
                                                "Stimulated, 24h" = "red", 
                                                "Stimulated, 48h" = "green",
                                                "Stimulated, 48h + LV" = "yellow",
                                                "Stimulated, 72h + LV" = "black")) + 
      labs(x="Condition", y = "Normalized expression (read counts) +/- s.d.") +
      ggtitle(label = paste(input$gene, "expression")) + 
      theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 30, hjust=1))
    
    print(p_gene)
  })
  
  
  output$plot <- renderPlot({
    print(plotInput())
    ggsave("plot.png", plotInput(), dpi = 300)
  })
  
  output$dndPlot <- downloadHandler(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      file.copy("plot.png", file, overwrite=TRUE)
    }
  ) 
}
