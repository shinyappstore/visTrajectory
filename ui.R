#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(plotly)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    useShinyjs(),
    # Application title
    titlePanel("Discover Ordering of the Data Points"),
      
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(position = "right", width = 15,
                   fluidRow(
                     column(4, selectInput(inputId = "data_type", 
                                           label =  h5(strong("Sequencing data type:")), 
                                           choices = list("microbiome" = "microbiome", 
                                                       "gene expression" = "gene"), 
                                           selected = "gene")),
                     column(4, textInput("covariate", label = h5(strong("Covariate")), 
                                         value = "hpf")),
                     column(4, numericInput("nCenters", value = 50,
                                            label = h5(strong("Number of points to highlight:"))))
                   ), 
                   fluidRow(
                     column(4, fileInput(inputId = "file_countTable", 
                               label = h5(strong("Count table file input (.csv)")),
                               accept = ".csv")),
                     column(4, selectInput(inputId = "method", 
                                 label =  h5(strong("Visualization method:")), 
                                 choices = list("PCoA" = "PCoA", "tSNE" = "tSNE"), 
                                 selected = "PCoA")),
                     column(2, numericInput("K", value = NULL,
                                            label = h5(strong("k-nearest neighbors:")))),
                     column(2, numericInput("nPaths", value = 50,
                                            label = h5(strong("Number of paths to plot:"))))
                   ), 
                   fluidRow(
                     column(4, fileInput(inputId = "file_sampleData", 
                                         label = h5(strong("Sample data file input (.csv)")),
                                         accept = ".csv")),
                     column(4, selectInput(inputId = "transform_distances", 
                                           label =  h5(strong("Apply dissimilarity transformation:")), 
                                           choices = list("yes" = TRUE, "no" = FALSE), 
                                           selected = TRUE)),
                     column(2, textInput("sample_idx", label = h5(strong("Chosen samples")), 
                                         value = "15, 30, 50")),
                     column(2, textInput("feat_idx", label = h5(strong("Chosen features")), 
                                         value = ""))
                   )
      ),
      # Show a plot of the generated distribution
      mainPanel(width = 20,
                fluidRow(
                  column(3, offset = 1, actionButton("rerunButton", "Rerun!")),
                  column(6, h4(id = "text"))
                ),
                fluidRow(
                  column(4, plotOutput("plot_rank_tau")),
                  column(4, plotOutput("plot_data_vs_tau")),
                  column(4, plotOutput("plot_X"))
                ),
                fluidRow(
                  column(3, offset = 1, downloadButton(outputId = "down_rank_tau", label = "Download plot")),
                  column(3, offset = 1, downloadButton(outputId = "down_data_tau", label = "Download plot")),
                  column(3, offset = 1, downloadButton(outputId = "down_plot_X", label = "Download plot"))
                ),
                fluidRow(
                  column(4, plotOutput("plot2D")),
                  column(4, plotOutput("plot_density")),
                  column(4, plotOutput("plot_contours"))
                ),
                fluidRow(
                  column(3, offset = 1, downloadButton(outputId = "down_buds2D", label = "Download plot")),
                  column(3, offset = 1, downloadButton(outputId = "down_density", label = "Download plot")),
                  column(3, offset = 1, downloadButton(outputId = "down_contours", label = "Download plot"))
                ),
                fluidRow(
                  column(6, plotlyOutput("plot3D")),
                  column(6, plotOutput("plot_features"))
                ),
                fluidRow(
                  column(3, offset = 7, downloadButton(outputId = "down_features", label = "Download plot"))
                )
      )
    )
  )
)
