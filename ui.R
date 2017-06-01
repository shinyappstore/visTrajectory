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
                     column(4, selectInput(inputId = "dist_method", 
                                           label =  h5(strong("Distance metric:")), 
                                           choices = list("Jaccard" = "jaccard", 
                                                       "correlation" = "correlation",
                                                       "Manhattan" = "manhattan",
                                                       "Euclidean" = "euclidean",
                                                       "Exp Manhattan" = "exp manhattan",
                                                       "Exp Euclidean" = "exp euclidean",
                                                       "Maximum" = "maximum",
                                                       "Canberra" = "canberra",
                                                       "Binary" = "binary",
                                                       "Minkowski"= "minkowski"), 
                                           selected = "correlation")),
                     column(4,  uiOutput("covariate")), 
                     column(2, selectInput(inputId = "init", 
                                           label =  h5(strong("Initialize from:")), 
                                           choices = list("Principal Curve" = "principal_curve", 
                                                          "Random" = "random"), 
                                           selected = "prin_curve")),
                     column(2, numericInput("K", value = NA,
                                            label = h5(strong("k-nearest neighbors:"))))
                   ), 
                   fluidRow(
                     column(4, fileInput(inputId = "file_countTable", 
                               label = h5(strong("Count table file input (.csv)")),
                               accept = ".csv")),
                     column(4, selectInput(inputId = "method", 
                                 label =  h5(strong("Visualization method:")), 
                                 choices = list("PCoA" = "PCoA", "tSNE" = "tSNE"), 
                                 selected = "PCoA")),
                     column(2, numericInput("nCenters", value = 50,
                                            label = h5(strong("No. highlighted of points:")))),
                     column(2, numericInput("nPaths", value = 50,
                                            label = h5(strong("No. of plotted paths:"))))
                   ), 
                   fluidRow(
                     column(4, fileInput(inputId = "file_sampleData", 
                                         label = h5(strong("Sample data file input (.csv)")),
                                         accept = ".csv")),
                     column(2, selectInput(inputId = "log_transform_data", 
                                           label =  h5(strong("Log transform data:")), 
                                           choices = list("yes" = TRUE, "no" = FALSE), 
                                           selected = TRUE)),
                     column(2, selectInput(inputId = "transform_distances", 
                                           label =  h5(strong("Transform dissimilarities:")), 
                                           choices = list("yes" = TRUE, "no" = FALSE), 
                                           selected = FALSE)),
                     column(2, textInput("sample_idx", label = h5(strong("Chosen samples")), 
                                         value = "15, 30, 50")),
                     column(2, textInput("feat_idx", label = h5(strong("Chosen features")), 
                                         value = ""))
                   )
      ),
      # Show a plot of the generated distribution
      mainPanel(width = 20,
                fluidRow(
                  column(2, offset = 1, actionButton("loadDefault", "Load default data!")),
                  column(2, actionButton("runButton", "Run BUDS!")),
                  column(2, actionButton("updateButton", "Update colors!")),
                  column(5, h4(id = "text"))
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
