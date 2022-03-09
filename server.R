library(shiny)
library(ggplot2)
library(tidyverse)
library(matrixStats) 


server <- function(input, output) {
  
  DataInput <- reactive({
    #Read data file and convert to a matrix 
    normalized_tibble_with_symbols <- read_csv("normalized_tibble_with_symbols.csv")
    matsymbol <- as.matrix(normalized_tibble_with_symbols[, 2:21])
    row.names(matsymbol) <-  normalized_tibble_with_symbols$...1
    return(matsymbol)
  }) 
  
    plotInput1 <- reactive({
    ## The gene is found in the matsymbol to extract all normalized read counts 
    ## A matrix of dim=5x4 is formed 
    gene_counts <- t(matrix(DataInput()[input$gene,], nrow=4))
    
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
  
  
  plotInput2 <- reactive({
    ## The gene is found in the matsymbol to extract all normalized read counts 
    ## A matrix of dim=5x4 is formed 
    gene_counts <- t(matrix(DataInput()[input$gene,], nrow=4))
    values <- as.vector(t(gene_counts)) 
    names <- c(rep("Non-stimulated", 4),rep("Stimulated, 24h", 4),rep("Stimulated, 48h", 4),rep("Stimulated, 48h + LV", 4),rep("Stimulated, 72h + LV", 4))
    count_df <- data.frame(value = values, name = names)
    
    
    b_gene <- ggplot(count_df, aes(x=names, y=values, fill = names)) +
      geom_boxplot(color="#4d4d4d", lwd = 0.2) + 
      geom_jitter(alpha=0.6, width=0.00, color = "#8c8c8c") + 
      scale_fill_manual("Condition", 
                        values = c("Non-stimulated" = "blue", 
                                   "Stimulated, 24h" = "red", 
                                   "Stimulated, 48h" = "green",
                                   "Stimulated, 48h + LV" = "yellow",
                                   "Stimulated, 72h + LV" = "black")) + 
      labs(x="Condition", y = "Normalized expression (read counts)") +
      ggtitle(label = paste(input$gene, "expression")) + 
      theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 30, hjust=1))
    
    print(b_gene)
    
  })
  
  
  output$plot1 <- renderPlot({
    print(plotInput1())
    ggsave("plot1.png", plotInput1(), dpi = 300)
  })
  
  output$plot2 <- renderPlot({
    print(plotInput2())
    ggsave("plot2.png", plotInput2(), dpi = 300)
  })
  
  output$dndPlot1 <- downloadHandler(
    filename = function() {
      "plot1.png"
    },
    content = function(file) {
      file.copy("plot1.png", file, overwrite=TRUE)
    }
  )
  
  output$dndPlot2 <- downloadHandler(
    filename = function() {
      "plot2.png"
    },
    content = function(file) {
      file.copy("plot2.png", file, overwrite=TRUE)
    }
  )
}
