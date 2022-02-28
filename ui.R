library(shiny)
library(tidyverse)
library(ggplot2)
library(matrixStats) 

normalized_tibble_with_symbols <- read_csv("normalized_tibble_with_symbols.csv")
matsymbol <- as.matrix(normalized_tibble_with_symbols[, 2:21])
row.names(matsymbol) <-  normalized_tibble_with_symbols$...1
gene_names <- sort(row.names(matsymbol))



ui <- fluidPage(
  
  titlePanel(title = span(img(src = "citco.png", height = 35), "T cell activation experiment"), windowTitle = "T cell activation experiment"),
  
  sidebarLayout(
    
    sidebarPanel(
      selectInput("gene", "Please enter a gene of interest (Examples: IL2, CD28, LDLR):", gene_names),
      p("Use this app to visualize gene expression data for activated T cells and T cells exposed to VSV-G pseudotyped lentivirus (LV)"),
      p(""),
      strong("Experimental design:"),
      p("- 4 donors"),
      p("- T cells stimulated with CD3/CD28 beads at cell:bead ratio of 1:1"),
      p("- T cells kept at 1E6 cells/mL"),
      p("- Lenti viral vector (VSV-G pseudotyped) added at MOI 10"),
      p("- Total RNA extraction using RNeasy Mini Kit (Qiagen)"),
      width = 4
    ),
    
    mainPanel(
      plotOutput("plot", width = "700px"),
      downloadButton("dndPlot", "Download Graph")
    )
  )
)
