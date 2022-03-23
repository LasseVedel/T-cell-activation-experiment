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
  
  DataInput2 <- reactive({
    
    log2FCs_df <- read_csv("log2FCs_df.csv")
    log2FCs_mat <- as.matrix(log2FCs_df[,1:7])
    row.names(log2FCs_mat) <- log2FCs_df$names
    return(log2FCs_mat)
    
  })
  
  DataInput3 <- reactive({
    
    padjust_df <- read_csv("padjust_df.csv")
    padjust_mat <- as.matrix(padjust_df[,1:7])
    row.names(padjust_mat) <- padjust_df$names
    return(padjust_mat)
    
  })
  
  
  plotInput1 <- reactive({
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
      theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 30, hjust=1))
    
    print(b_gene)
    
  })
  
  
  plotInput2 <- reactive({
    
    gene <- input$gene
    
    #Extracting values for the gene of interest
    gene_vals <- as.numeric(DataInput2()[gene,])
    p_vals <- as.numeric(DataInput3()[gene,])
    p_vals <- format(p_vals, digits=3,scientific = T)
    
    #Creating data frame that contains information on specific gene 
    df_to_plot <- as.data.frame(cbind("Conditions" = colnames(DataInput2()),
                                      "Log2FCs" = as.numeric(gene_vals)))
    #defining levels
    levels <- c("Stimulated, 24h vs. Non-stimulated",
                "Stimulated, 48h vs. Non-stimulated",
                "Stimulated, 48h + LV vs. Non-stimulated",
                "Stimulated, 72h + LV vs. Non-stimulated",
                "Stimulated, 48h vs. Stimulated, 24h",
                "Stimulated, 48h + LV vs. Stimulated, 48h",
                "Stimulated, 72h + LV vs. Stimulated, 48h + LV")
    
    df_to_plot$Conditions <- factor(df_to_plot$Conditions,                                    
                                    levels = levels)
    
    ## Plot the Log2FCs
    L2FC_p <- ggplot(df_to_plot, aes(x=Conditions, y=as.numeric(Log2FCs))) +
      geom_point(size=3, fill=alpha("black", 1.0), alpha=0.7, shape=21, stroke=2) +
      geom_segment( aes(x=Conditions, xend=Conditions, y=0, yend=as.numeric(Log2FCs))) + 
      geom_hline(yintercept=0, linetype="dashed", color = "red", size=2) + 
      labs(x="Tested conditions", y = "Log2 fold changes (Log2FC)",
           subtitle = "Shrunken Log2FCs, p-values obtained by Wald testing (Benjamini-Hochberg adjusted)") +
      ggtitle(label = paste("Log2 fold changes for", gene)) + 
      theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 30, hjust=1)) + 
      theme(
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 10, margin=margin(b = 10, unit = "pt")),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(margin=margin(b = 20, unit = "pt"))) +
      annotate("label", 
             x = levels, 
             y = c((max(gene_vals)-min(gene_vals))/10+gene_vals), 
             label = c(paste("p =", p_vals[1]),
                       paste("p =", p_vals[2]),
                       paste("p =", p_vals[3]),
                       paste("p =", p_vals[4]),
                       paste("p =", p_vals[5]),
                       paste("p =", p_vals[6]),
                       paste("p =", p_vals[7])))
    
    print(L2FC_p)
    
  })
  
  
  
  #Generating the table displaying log2FCs and p-values
  tableInput1 <- reactive({
    
    levels <- c("Stimulated, 24h vs. Non-stimulated",
                "Stimulated, 48h vs. Non-stimulated",
                "Stimulated, 48h + LV vs. Non-stimulated",
                "Stimulated, 72h + LV vs. Non-stimulated",
                "Stimulated, 48h vs. Stimulated, 24h",
                "Stimulated, 48h + LV vs. Stimulated, 48h",
                "Stimulated, 72h + LV vs. Stimulated, 48h + LV")
    gene <- input$gene
    
    #Extracting values for the gene of interest
    gene_vals <- as.numeric(DataInput2()[gene,])
    gene_vals_digits <- format(gene_vals, digits=2,scientific = F)
    p_vals <- as.numeric(DataInput3()[gene,])
    p_vals <- format(p_vals, digits=3,scientific = T)
    
    #Creating the table 
    table_with_values <- as.data.frame(cbind("Shrunken log2FCs" = gene_vals_digits,
                                             "Adjusted p-values" = p_vals))
    
    row.names(table_with_values) <- levels
    table_with_values <- t(table_with_values)
    return(table_with_values)
    
  })
  
  
  
  
  #
  #All of the outputs generated from this script 
  #
  
  #Plot for the normalized expression 
  output$plot1 <- renderPlot({
    print(plotInput1())
    ggsave("plot1.png", plotInput1(), dpi = 300)
  })
  
  #Plot for the log2FC plot 
  output$plot2 <- renderPlot({
    print(plotInput2())
    ggsave("plot2.png", plotInput2(), dpi = 300)
  })
  
  #Download handler for normalized expression plot 
  output$dndPlot1 <- downloadHandler(
    filename = function() {
      "plot1.png"
    },
    content = function(file) {
      file.copy("plot1.png", file, overwrite=TRUE)
    }
  )
  
  #Download handler for log2FC plot 
  output$dndPlot2 <- downloadHandler(
    filename = function() {
      "plot2.png"
    },
    content = function(file) {
      file.copy("plot2.png", file, overwrite=TRUE)
    }
  )
  
  #Output for log2FC/p-value table 
  output$table1 <- renderTable({
    tableInput1()
  },
  width = "900px",
  rownames = T)
  
}
