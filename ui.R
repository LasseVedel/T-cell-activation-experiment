library(shiny)
library(tidyverse)

normalized_tibble_with_symbols <- read_csv("normalized_tibble_with_symbols.csv")
matsymbol <- as.matrix(normalized_tibble_with_symbols[, 2:21])
row.names(matsymbol) <-  normalized_tibble_with_symbols$...1
gene_names <- sort(row.names(matsymbol))



ui <- fluidPage(
  selectInput("gene", "Please enter a gene of interest (Examples: IL2, CD28, LDLR):", gene_names),
  plotOutput("plot")
)
