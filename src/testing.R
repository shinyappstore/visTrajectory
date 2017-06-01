setwd("~/MEGA/BIOSTATS GROUP/Projects/visBUDS/visTrajectory/")
rm(list = ls())
source("R/distcomps.R")
source("R/get_data_to_plot.R")
source("R/plot_utils.R")
library(shiny)
library(shinyjs)
library(plotly)


# Packages
sapply(c("buds", "coda", "DistatisR", "dplyr", "ggplot2", "MCMCglmm", 
         "plotly", "plyr", "princurve", "reshape2", "rstan", "shiny",
         "viridis"), require, character.only = TRUE)

# Options
source("./R/distcomps.R")
source("./R/get_data_to_plot.R")
source("./R/plot_utils.R")


# Packages
sapply(c("buds", "coda", "DistatisR", "dplyr", "ggplot2", "MCMCglmm", 
         "plotly", "plyr", "princurve", "reshape2", "rstan", "shiny",
         "viridis"), require, character.only = TRUE)

# Options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(shiny.maxRequestSize=30*1024^2) 
theme_set(theme_classic())
theme_update(text=element_text(size=20))

# Parameters
min_row_sum <- 100
min_row_prevalence <- 5
B <- 100
min_sigma <- 0.05
hparams <- list(
  "gamma_tau"= 2.5,
  "gamma_epsilon" = 2.5,
  "gamma_bias" = 2.5,
  "gamma_rho" = 2.5,
  "min_sigma" = min_sigma
)

# Default data files
countTable_default_file <- "data/114_US_Senate_binVotes.csv"
sampleData_default_file <- "data/114_US_Senate_legisData.csv"

X <- read.csv(countTable_default_file, row.names = 1)
sampleData <- read.csv(sampleData_default_file, row.names = 1)
covariate_name <- "party"
sample_covariate <- sampleData[, covariate_name]

K <- floor(ncol(X)/10)

#D0 <- cor_dist(X, log_trans = TRUE, scale = FALSE, base = 1e6)
#D0 <- generic_dist(X, method = "jaccard")
D0 <- generic_dist(X, method = "exp manhattan", log_trans = FALSE)
D <- D0
#D <- transform_dist(D0, threshold = F)

budsFit <- buds::fit_buds(D, K = 10, method = "vb", hyperparams = hparams,
                          init_from = "random")
budsParams <- (rstan::extract(budsFit$fit_buds))
tau_df <- get_tau_df(budsParams, prob = 0.95)


theme_set(theme_bw())
party_cols <- c("D" = "#1f78b4", "R" = "#e31a1c", "Indep" = "#ff7f00")
sampleData$Legis <- rownames(sampleData)
sampleData$Legis <- factor(sampleData$Legis, 
                           levels = sampleData$Legis[order(-tau_df$tau)])
(ptau1 <- tau_df %>%
    ggplot(aes(x = sampleData$Legis, y = tau, col = sampleData$party)) +
    geom_errorbar(aes(ymin = tau_lower, ymax = tau_upper), width = 1, lwd = 0.9) +
    geom_point(pch = 21, aes(fill = sampleData$party), color = "white", size = 3) + 
    scale_fill_manual(name = "Party", values = party_cols) +
    scale_color_manual(name = "Party", values = party_cols) + 
    ylab("BUDS tau") + xlab("Senator") + 
    theme(text = element_text(size = 20), legend.position = c(0.9, 0.8), 
          axis.text.x = element_text(angle = 90, size = 8, face = "bold")))

theme_set(theme_classic())
theme_update(text = element_text(size = 20))
plt <- plot_buds_1D(tau_df, covariate = NULL,
                    color = sample_covariate, 
                    color_label = covariate_name, 
                    idxBigger = NULL) 
plt + geom_errorbar(aes(ymin = tau_lower, ymax = tau_upper), lwd = 0.9, width = 2) +
  geom_point(aes(fill = color), color = "white", pch = 21, size = 3) +
  scale_color_manual(name = "Party", values = party_cols) +
  scale_fill_manual(name = "Party", values = party_cols) 


################################################################################


theme_set(theme_bw())
(diab_ptau1 <- tau_df %>%
    ggplot(aes(x = sampleData$Age_at_Collection, y = tau, col = sampleData$Age_at_Collection)) +
    geom_errorbar(aes(ymin = tau_lower, ymax = tau_upper), width = 1, lwd = 0.9) +
    geom_point(pch = 21, aes(fill = sampleData$Age_at_Collection), color = "white", size = 3) + 
    scale_fill_viridis(name = "Age [days]") +
    scale_color_viridis(name = "Age [days]") + 
    ylab("BUDS tau") +
    theme(text = element_text(size = 20)))


(diab_ptau2 <- tau_df %>%
    ggplot(aes(x = rank(tau), y = tau, col = sampleData$Age_at_Collection)) +
    geom_errorbar(aes(ymin = tau_lower, ymax = tau_upper), width = 5, lwd = 0.7) +
    geom_point(pch = 21, aes(fill = sampleData$Age_at_Collection), color = "white", size = 3) + 
    scale_fill_viridis(name = "Age [days]") +
    scale_color_viridis(name = "Age [days]") + 
    theme(text = element_text(size = 20)))

theme_set(theme_classic())
theme_update(text = element_text(size = 20))
plt <- plot_buds_1D(tau_df, covariate = NULL,
                    color = sample_covariate, 
                    color_label = covariate_name, 
                    idxBigger = NULL) 
plt + geom_errorbar(aes(ymin = tau_lower, ymax = tau_upper), lwd = 0.7, width = 7) +
  geom_point(aes(fill = color), color = "grey77", pch = 21, size = 2) 
  

boot <- get_D_copies(D, budsFit, B)
distatis_input <- get_input_for_distatis(D = D, 
                                         D.lst = boot$D.lst,
                                         tau_mode =  tau_df$tau,
                                         tau.lst = boot$booData.lst, 
                                         sample_data = sampleData)




res <- run_distatis(bootD = distatisData$bootD, dims =2,
                    booData.lst = distatisData()$booData.lst, 
                    modeData = distatisData()$modeData)



