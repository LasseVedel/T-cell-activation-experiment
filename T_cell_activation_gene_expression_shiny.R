library(shiny)
library(dplyr)
library(ggplot2)


normalized_tibble_with_symbols <- read_csv("normalized_tibble_with_symbols.csv")
matsymbol <- as.matrix(normalized_tibble_with_symbols[, 2:21])
row.names(matsymbol) <-  normalized_tibble_with_symbols$...1
gene_names <- sort(row.names(matsymbol))


ui <- fluidPage(
  selectInput("gene", "Please enter a gene of interest (Examples: IL2, CD28, LDLR):", gene_names),
  plotOutput("plot")
)

server <- function(input, output) {
  
  data1 <- reactive({
    #Read data file and convert to a matrix 
    normalized_tibble_with_symbols <- read_csv("normalized_tibble_with_symbols.csv")
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
    row_sds_gene <- rowSds(gene_counts)
    
    ## Collect to a dataframe which can be used for ggplot  
    df_gene <- as.data.frame(cbind(row_means_gene, row_sds_gene))
    return(df_gene)
    
  })
  
  
  output$plot <- renderPlot({
    
    req(data1())
    
    ## Plot the expression using ggplot 
    p_gene <- ggplot(data1(), aes(x=rownames(data1()), y=row_means_gene, fill = rownames(data1()))) + 
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
    p_gene
  })
}
shinyApp(ui, server)

