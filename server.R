#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
# Import helper functions
source("R/distcomps.R")
source("R/get_data_to_plot.R")
source("R/plot_utils.R")


# Packages
sapply(c("buds", "coda", "DistatisR", "dplyr", "ggplot2", "MCMCglmm", 
         "plotly", "plyr", "princurve", "reshape2", "rstan", "shiny",
         "viridis"), require, character.only = TRUE)

# Options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(shiny.maxRequestSize=30*1024^2) 
theme_set(theme_bw())

# Parameters
min_row_sum <- 100
min_row_prevalence <- 5
B <- 100
init_buds_from <- "principal_curve"

# Default data files
countTable_default_file <- "data/frog_processed_counts.csv"
sampleData_default_file <- "data/frog_sample_data.csv"
frog_D_resfile <- "data/frog_buds_res.rds"
frog_D0_resfile <- "data/frog_buds_res_untrans.rds"

shinyServer(function(input, output) {
  shinyjs::html("text", "")
  
  # Count table
  X <- reactive({
    if (!is.null(input$file_countTable)) {
      message("Loading count table...")
      shinyjs::html("text", "Loading count table...")
      inFile <- input$file_countTable
      X <- read.csv(inFile$datapath, row.names = 1)
      return(X)
    }
    if (input$loadDefault) {
      message("Loading default count table...")
      shinyjs::html("text", "Loading default count table...")
      X <- read.csv(countTable_default_file, row.names = 1)
      return(X)
    }
  })
  
  # External attributes of the sample
  sampleData <- reactive({
    if (!is.null(input$file_sampleData)) {
      message("Loading sample data...")
      shinyjs::html("text", "Loading sample data...")
      inFile <- input$file_sampleData
      sampleData <- read.csv(inFile$datapath, row.names = 1)
      return(sampleData)
    }
    if (input$loadDefault) {
      message("Loading default sample data...")
      shinyjs::html("text", "Loading default sample data...")
      sampleData <- read.csv(sampleData_default_file, row.names = 1)
      return(sampleData)
    }
  })
  
  # Drop-down selection box for which data set
  output$covariate <- renderUI({
    selectInput("covariate", "Sample covariate:", colnames(sampleData()))
  })
  
  covariate_name <- eventReactive(input$updateCovariate, {
    if(input$covariate %in% colnames(sampleData())) {
      return(input$covariate)
    } else if (input$loadDefault & input$covariate == "") {
      return("hpf")
    } else {
      return(NA)
    }
  })
  
  # Covariate used for coloring points
  sample_covariate <- eventReactive(covariate_name(), {
    if(!(covariate_name() %in% colnames(sampleData()))) {
      warning("Selected covariate is not in the sample data.")
      shinyjs::html("text", "ERROR: Selected covariate is not in the sample data.")
      return(rep(1, ncol(X())))
    }
    sample_covariate <- sampleData()[, covariate_name()]
    return(sample_covariate)
  })
  
  # Chosen features indexes
  chosen_feats <- reactive({
    if(is.null(input$feat_idx) || input$feat_idx == "") {
      return(NA)
    }
    features <- strsplit(input$feat_idx, ",")[[1]]
    if (all(features %in% rownames(X()))) {
      return(features)
    }
    idx <- as.numeric(features)
    idx <- idx[idx %in% 1:nrow(X())]
    return(idx)
  })
  
  # Chosen sample indexes
  chosen_samples <- reactive({
    samples <- strsplit(input$sample_idx, ",")[[1]]
    if (all(samples %in% colnames(X()))) {
      return(samples)
    }
    idx <- as.numeric(samples)
    idx <- idx[idx %in% 1:ncol(X())]
    return(idx)
  })
  
  # Constant (K) for selecting number of kNN
  K <- reactive({
    if(is.null(input$K)) {
      K <- floor(ncol(X())/10)
    } else {
      K <- input$K
    }
    return(K)
  })
  
  # Dissimilarity matrix
  D0 <- eventReactive(c(input$runButton, input$loadDefault), {
    req(input$runButton | input$loadDefault)
    message("Computing dissimilarities...")
    shinyjs::html("text", "Computing dissimilarities...")
    typeData <- input$data_type
    if (typeData == "microbiome") {
      D <- jacc_dist(X(), min_row_sum, min_row_prevalence)
    } else {
      D <- cor_dist(X(), log_trans = TRUE, scale = TRUE, base = 1e6)
    }
    return(D)
  })
  
  
  # Transformed dissimilarity matrix
  D <- reactive({
    req(D0())
    if (input$transform_distances) {
      message("Transforming dissimilarities...")
      shinyjs::html("text", "Transforming dissimilarities...")
      D <- transform_dist(D0())
    }
    return(D)
  }) 
  
  # Fit 1D latent coordinates, tau, with BUDS
  budsFit <- eventReactive(c(input$runButton, input$loadDefault), {
    if(input$loadDefault & is.null(isolate(input$file_countTable))){
      res_file <- ifelse(input$transform_distances,
                         frog_D_resfile, frog_D0_resfile)
      buds_fit <- readRDS(res_file)
      return(buds_fit)
    }
    req(D())
    message("Fitting BUDS model...")
    shinyjs::html("text", "Fitting BUDS model...")
    fit <- buds::fit_buds(D(), K = K(), method = "vb",
                          init_from = init_buds_from)
    buds_fit <- fit$fit_buds
    return(buds_fit)
  })
  
  # Gather tau samples drawn from the posterior
  tau_df <- reactive({
    req(budsFit())
    return(get_tau_df(budsFit(), prob = 0.95))
  })
  
  # Choose a subset of samples to display trajectory for
  idxBigger <- reactive({
    DF <- tau_df()
    rownames(DF) <- 1:nrow(DF)
    ordTau_df <- DF[order(DF$tau), ]
    idx <- seq(1, nrow(DF), length.out = input$nCenters)
    idx <- as.numeric(rownames(ordTau_df)[idx])
    return(idx) 
  })
  
  Y2D <- reactive({
    return(low_dim_vis(D0(), method = input$method, dims = 2))
  })
  
  Y3D <- reactive({
    return(low_dim_vis(D0(), method = input$method, dims = 3))
  })
  
  distatisData <- reactive({
    req(budsFit())
    fitParams <- rstan::extract(budsFit())
    nDraws <- nrow(fitParams$tau)
    B <- min(B, nDraws) 
    boot <- get_D_copies(D(), budsFit(), B)
    distatis_input <- get_input_for_distatis(D = D(), 
                                             D.lst = boot$D.lst,
                                             tau_mode =  tau_df()$tau,
                                             tau.lst = boot$booData.lst, 
                                             sample_data = isolate(sampleData()))
    message("Computing input data for DiSTATIS...")
    shinyjs::html("text", "Computing input data for DiSTATIS...")
    return(distatis_input) 
  })
  
  distatis_res <- reactive({
    res <- run_distatis(bootD = distatisData()$bootD, dims =2,
                        booData.lst = distatisData()$booData.lst, 
                        modeData = distatisData()$modeData)
    message("Running DiSTATIS...")
    shinyjs::html("text", "Running DiSTATIS...")
    return(res)
  })
  
  
  ####################### PLOTTING #########################
  
  # # Plot average distances to kNN along tau interval
  # plot_kNN_dist <- reactive({
  #   req(covariate_name())
  #   plt <- plot_kNN_mean_dist(D(), tau_df()$tau, color = sample_covariate(),
  #                             color_label = covariate_name(), 
  #                             K = input$K, error_bars = FALSE)
  #   return(plt)
  # })
  
  # Plot tau estimates vs rankof tau
  plot_rank_tau <- eventReactive(c(tau_df(), sample_covariate()), {
    plt <- plot_buds_1D(tau_df(), covariate = NULL,
                        color = sample_covariate(), 
                        color_label = covariate_name(), 
                        idxBigger = idxBigger()) 
    return(plt)
  })
  
  # Plot tau against chosen covariate
  plot_data_vs_tau <- eventReactive(c(tau_df(), sample_covariate()), {
    plt <- plot_buds_1D(tau_df(), covariate = sample_covariate(),
                        color = sample_covariate(), 
                        color_label = covariate_name(), 
                        idxBigger = idxBigger()) 
    return(plt)
  })
  
  # Plot 2D visualization of the data and trajectories
  plot2D <- eventReactive(c(Y2D(), covariate_name()), {
    plt <- plot_buds_trajectory(budsFit(), Y2D()$Y, Y2D()$eigs,
                                sample_data = sampleData(), 
                                covariate_name = covariate_name(), 
                                path_col = "#2171B5", 
                                nPaths = input$nPaths, 
                                nCenters = input$nCenters)
    return(plt)
  })
  
  # Plot 3D visualization of the data and trajectory on idxBigger
  plot3D <- eventReactive(c(Y3D(), covariate_name()), {
    plt <- plot_buds_trajectory(budsFit(), Y3D()$Y, Y3D()$eigs,
                                sample_data = sampleData(), 
                                covariate_name = covariate_name(), 
                                path_col = "#2171B5", 
                                nPaths = input$nPaths, 
                                nCenters = input$nCenters)
    return(plt)
  })
  
  plot_X <- eventReactive(c(input$runButton, input$loadDefault), {
    req(ncol(X()) == nrow(sampleData()))
    plt <- plot_ordered_matrix(X(), tau_df()$tau, log_trans = TRUE,
                               keep_fatures = NULL, 
                               nfeatures = max(500, 2*ncol(X())),
                               byMean = TRUE, window = NULL)
    return(plt)
  })
  
  plot_density <- eventReactive(c(distatis_res(), covariate_name()), {
    distatis_df <- distatis_res()$partial
    consensus_df <- distatis_res()$consensus
    plt <- plot_distatis(distatis_df, consensus_df, 
                         color_label = covariate_name()) 
    return(plt)
  })
  
  plot_contours <- eventReactive(c(distatis_res(), covariate_name()), {
    distatis_df <- distatis_res()$partial
    consensus_df <- distatis_res()$consensus
    plt <- plot_point_contours(distatis_df, consensus_df, 
                               idx_list = chosen_samples(), 
                               color_label = covariate_name())  
    return(plt)
  })
  
  plot_features <- eventReactive(c(input$runButton, input$loadDefault), {
    plt <- plot_features_curves(X(), tau_df()$tau, 
                                feat_idx = chosen_feats(), 
                                log_trans = TRUE)
    return(plt)
  })
  
  
  output$plot_features <- renderPlot({
    plot_features()
  })  
  
  output$plot_density <- renderPlot({
    plot_density()
  })  
  
  output$plot_contours <- renderPlot({
    plot_contours()
  })  
  
  output$plot_kNN_dist <- renderPlot({
    plot_kNN_dist()
  })  
  
  output$plot_rank_tau <- renderPlot({
    plot_rank_tau()
  })
  
  output$plot_data_vs_tau <- renderPlot({
    plot_data_vs_tau()
  })
  
  output$plot2D <- renderPlot({
    plot2D()
  })
  
  output$plot3D <- renderPlotly({
    plot3D()
  })
  
  output$plot_X <- renderPlot({
    plot_X()
  })  
  
  output$down_rank_tau <- downloadHandler(
    filename = "rank_tau.png",
    content = function(file) {
      ggsave(file, plot_rank_tau())
    }
  )
  
  output$down_data_tau <- downloadHandler(
    filename = "tau_vs_covariate.png",
    content = function(file) {
      ggsave(file, plot_data_vs_tau())
    }
  )
  
  output$down_plot_X <- downloadHandler(
    filename = "heatmap.png",
    content = function(file) {
      ggsave(file, plot_X())
    }
  )
  
  output$down_buds2D <- downloadHandler(
    filename = "trajectory2D.png",
    content = function(file) {
      ggsave(file, plot2D())
    }
  )
  
  output$down_density <- downloadHandler(
    filename = "data_density.png",
    content = function(file) {
      ggsave(file, plot_density())
    }
  )
  
  output$down_contours <- downloadHandler(
    filename = "contours.png",
    content = function(file) {
      ggsave(file, plot_contours())
    }
  )
  
  output$down_features <- downloadHandler(
    filename = "features.png",
    content = function(file) {
      ggsave(file, plot_features())
    }
  )
  
})



