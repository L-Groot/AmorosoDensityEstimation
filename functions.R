######################

# FUNCTIONS

######################

# dAmoroso()
# --> Amoroso density function

# LL_Amoroso()
# --> Likelihood Function for Amoroso

# BIC_Amoroso()
# --> Calculates the BIC for a given set of Amoroso parameters and data

# estimateAmoroso_MLE_MDE()
# --> Estimates the Amoroso to certain data vector using MLE and MDE methods
# --> Returns a dataframe with parameter estimates and BIC for all fits
# --> Plots Amoroso fits on histogram and together with R Kernel density estimate

#----------------------

library(AmoRosoDistrib)
library(palmerpenguins)
library(tidyverse)
library(ggplot2)
library(gridExtra)

#-----------------------------------------------
## Testing function
# source("estimate_amoroso_MLE_MDE.R")
# 
# library(AmoRosoDistrib)
# library(palmerpenguins)
# library(multimode)
# 
# dat <- palmerpenguins::penguins
# dat <- rnorm(100, 40, 5)
# 
# modetest(dat)
# modetest(dat, method = "HH")
# 
# hist(dat)
# 
# estimate_amoroso_MLE_MDE(vec = dat)
#-----------------------------------------------


#--------------------------------------------
# Define Amoroso probability density function
#--------------------------------------------
# a: scale
# lambda: shape
# c: shape
# mu: location

dAmoroso <- function(x, a, lambda, c, mu) {
  c1 <- 1/(gamma(lambda))
  c2 <- abs(c/a)
  c3 <- ((x-mu)/a)^(lambda*c-1)
  c4 <- exp(-((x-mu)/a)^c)
  return(c1*c2*c3*c4)
}

#--------------------------------------------------------
# Define log-likelihood function for Amoroso distribution
#--------------------------------------------------------
LL_Amoroso <- function(data, a, lambda, c, mu) {
  # Calculate the log-likelihood for each data point
  log_likelihood_values <- log(dAmoroso(data, a, lambda, c, mu))
  
  # Sum the log-likelihood values
  log_likelihood <- sum(log_likelihood_values)
  
  return(log_likelihood)
}


#--------------------------------
# Define BIC function for Amoroso
#--------------------------------
BIC_Amoroso <- function(data, params) {
  # Extract parameters
  a <- params[1]
  lambda <- params[2]
  c <- params[3]
  mu <- params[4]
  
  # Calculate log-likelihood
  log_likelihood <- LL_Amoroso(data, a, lambda, c, mu)
  
  # Number of parameters in the model
  num_params <- length(params)
  
  # Number of data points
  n <- length(data)
  
  # Calculate BIC
  bic <- -2 * log_likelihood + num_params * log(n)
  
  return(bic)
}


#-------------------------------------------------------------------
# Define short info texts about plots for estimate_amoroso_MLE_MDE()
#-------------------------------------------------------------------
one_plot_info <- "\nABOUT THE PLOT:\n
  The plot contains the histogram of the variable (grey bars),
  the nonparametric R Kernel density estimator fit (dark grey line) and the
  Amoroso fits from the initial parameter estimation (as describes by Combes
  et al. (2022). Each coloured line corresponds to one Amoroso fit per method.
  Since each method has two sets of parameter estimates (+ve and -ve space), the
  plot shows only the one that fits the data better according to the BIC. For
  some methods this may be the parameter set in negative space and for other
  methods this may be the parameter set in positive space.\n"
grid_plots_info <- "\nABOUT THE PLOTS:\n
  Each of the 9 plots shows the two Amoroso fits for each method: the fit from
  the parameter estimates in positive space (green line) and the fit from the
  parameter estimates in negative space (red line). The dark-grey line is the
  fit of the nonparametric R Kernel density.\n"

#-----------------------------------------------------------------
# Define function to fit Amoroso with MLE and MDE and plot results
#-----------------------------------------------------------------
estimateAmoroso_MLE_MDE <- function(vec=NULL,
                               dataframe=NULL, variable=NULL,
                               plotpermethod = FALSE,
                               breaks=20,
                               varname=NULL) {
  
  #----------
  # Define y
  #----------
  if (!is.null(vec) && is.vector(vec)) {
    # If vector is provided directly, use it as y
    y <- na.omit(vec)
    cat("Using the provided vector", deparse(substitute(vec)), "\n")
    # Check if any NAs were removed
    num_nas_removed <- length(vec) - length(na.omit(vec))
    if(num_nas_removed > 0) {
      message(c(as.character(num_nas_removed), " NAs were removed from the data."))
    }
  } else if (!is.null(dataframe) && !is.null(variable) && is.data.frame(dataframe) && is.character(variable)) {
    # If dataframe and variable are provided, extract variable from df and assign to y
    y <- na.omit(eval(substitute(dataframe[[variable]])))
    cat("Using the", variable, "variable from the", deparse(substitute(dataframe)), "dataframe. \n")
    # Check if any NAs were removed
    num_nas_removed <- length(eval(substitute(dataframe[[variable]]))) - length(na.omit(eval(substitute(dataframe[[variable]]))))
    if(num_nas_removed > 0) {
      message(num_nas_removed, " NAs were removed from the data.")
    }
  } else {
    stop("Invalid arguments. \n. Provide either \n (1) a vector object directly to the
         `vec` argument or \n (2) a dataframe object to the 'dataframe argument
         AND the name of the variable's column name as a string to the 'variable'
         argument. \n Otherwise check whether the objects, columns etc. you wish to use
         indeed exist and that spelling is correct. \n")
  }
  
  #-----------------------------------------
  # Initialize dataframe to store parameters
  #-----------------------------------------
  df <- data.frame(
    method = rep(c("Initial","MLE","Kullback-Leibler","Jensen-Shanon","Hellinger",
                   "Wasserstein","Squared","Hellinger (CDF)","Wasserstein (CDF)",
                   "Squared (CDF)"), each = 2),
    code = rep(c("init","mle","kull","jens","hell","wass","sq","hell_cdf",
                 "wass_cdf","sq_cdf"), each = 2),
    par_space = rep(c("+","-"), length.out = 20),
    a = NA,
    l = NA,
    c = NA,
    mu = NA
  )
  
  #-----------------------------
  # Find initializing parameters
  #-----------------------------
  # For a > 0
  init.pos = init.theta(data = y, -20, 20, length = 1000, a.pos = TRUE)[1:4]
  df[which(df$code == "init" & df$par_space == "+"), c("a","l","c","mu")] <-
    init.pos
  # For a < 0
  init.neg = init.theta(data = y, -20, 20, length = 1000, a.pos = FALSE)[1:4]
  df[which(df$code == "init" & df$par_space == "-"), c("a","l","c","mu")] <-
    init.neg
  
  #-----------------------------------------------------------------
  # Estimate parameters for MLE and MDE methods in +ve and -ve space
  #-----------------------------------------------------------------
  # MLE
  df[which(df$code == "mle" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mle(y, init.pos)$par
  df[which(df$code == "mle" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mle(y, init.neg)$par
  # Kullback-Leibler divergence (based on PDF)
  df[which(df$code == "kull" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mkle(y, init.pos)$par
  df[which(df$code == "kull" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mkle(y, init.neg)$par
  # Jensen-Shanon divergence (based on PDF)
  df[which(df$code == "jens" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mjse(y, init.pos)$par
  df[which(df$code == "jens" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mjse(y, init.neg)$par
  # Hellinger distance (based on PDF)
  df[which(df$code == "hell" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mhe(y, init.pos)$par
  df[which(df$code == "hell" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mhe(y, init.neg)$par
  # Wasserstein distance (based on PDF)
  df[which(df$code == "wass" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mwe(y, init.pos)$par
  df[which(df$code == "wass" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mwe(y, init.neg)$par
  # Squared distance (based on PDF)
  df[which(df$code == "sq" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.msqe(y, init.pos)$par
  df[which(df$code == "sq" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.msqe(y, init.neg)$par
  # Hellinger distance (based on CDF)
  df[which(df$code == "hell_cdf" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mhdfe(y, init.pos)$par
  df[which(df$code == "hell_cdf" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mhdfe(y, init.neg)$par
  # Wasserstein distance (based on CDF)
  df[which(df$code == "wass_cdf" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.mwdfe(y, init.pos)$par
  df[which(df$code == "wass_cdf" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.mwdfe(y, init.neg)$par
  # Squared distance (based on CDF)
  df[which(df$code == "sq_cdf" & df$par_space == "+"), c("a","l","c","mu")] <-
    fit.msqdfe(y, init.pos)$par
  df[which(df$code == "sq_cdf" & df$par_space == "-"), c("a","l","c","mu")] <-
    fit.msqdfe(y, init.neg)$par
  
  
  #-------------------------
  # Model selection with BIC
  #-------------------------
  df <- df %>%
    rowwise() %>%
    # Calculate BIC for each method and parameter space
    mutate(BIC = BIC_Amoroso(data = y, params = c(a, l, c, mu))) %>%
    ungroup() %>%
    # Per method, identify better parameter space set (+ or -)
    group_by(method) %>%
    mutate(win_pair = ifelse(BIC == min(BIC), 1, 0)) %>%
    ungroup() %>%
    # Identify model with lowest BIC of all
    mutate(win = ifelse(BIC == min(BIC), 1, 0))
  
  #-------
  # PLOTS
  #-------
  # Define x
  xx <- seq(min(y), max(y), length = 10000)
  # Calculate kernel density estimate
  dens <- density(y)
  
  if (plotpermethod == FALSE) {
    #--------------------------------------
    # One plot with best Amoroso per method
    #--------------------------------------
    
    # Print plot info to console
    cat(one_plot_info,"\n")
    
    # Create plot title
    if (!is.null(varname)) {
      title <- paste("Histogram of", varname)
    } else if (is.null(variable) && !is.null(vec)) {
      title <- paste("Histogram of", deparse(substitute(vec)))
    } else {
      title <- paste("Histogram of", variable)
    }
    
    # Subset df for only the winning models (but keep both initial)
    df_win <- df %>% filter(code == "init" | win_pair == 1)
    
    # Define unique code to separate init+ from init-
    df_win$code_unique <- paste0(df_win$code,df_win$par_space)
    
    
    ### Add the lower BIC model of each method to histogram
    
    # Create empty list to store densities in
    dens_list <- vector("list", length(df_win$code_unique))
    
    names(dens_list) <- df_win$code_unique
    
    # Extract the 4 parameters
    for (i in 1:nrow(df_win)) {
      parvec <- c(
        (df_win[i,"a"] %>% pull()),
        (df_win[i,"l"] %>% pull()),
        (df_win[i,"c"] %>% pull()),
        (df_win[i,"mu"]) %>% pull())
      
      # Calculate density for current method
      dens_obj <- dgg4(xx,
                       a = parvec[1],
                       l = parvec[2],
                       c = parvec[3],
                       mu = parvec[4])
      
      # Save density for method in a new object and assign to list
      index_method <- df_win$code_unique[i]
      dens_list[[index_method]] <- dens_obj
    }
    
    # Put densities of all methods in a dataframe
    dens_df <- data.frame(matrix(NA, ncol = length(df_win$code_unique), nrow = 10000))
    colnames(dens_df) <- df_win$code_unique
    for (i in 1:ncol(dens_df)) {
      dens_df[,i] <- dens_list[[i]]
    }
    
    # Add x values
    dens_df$x <- xx
    
    # Transform to long format for plotting
    dens_long <- dens_df %>%
      pivot_longer(cols = -x, names_to = "Method", values_to = "Density")
    
    # Make histogram with Amoroso fits
    theplot <- ggplot() +
      # Create histogram
      geom_histogram(aes(x = y, y = ..density..), bins = 30, fill = "grey77", color = "grey40") +
      # Add standard nonparametric kernel density estimate
      geom_density(aes(x=y, colour = "lightgrey")) +
      # Add Amoroso fits
      geom_line(data = dens_long, aes(x = x, y = Density, color = Method), size = 1.2) +
      scale_color_manual(values = c("init+" = "black",
                                    "init-" = "dimgrey",
                                    "hell_cdf+" = "firebrick1",
                                    "hell_cdf-" = "firebrick1",
                                    "hell+" = "darkorange",
                                    "hell-" = "darkorange",
                                    "jens+" = "gold",
                                    "jens-" = "gold",
                                    "kull+" = "darkolivegreen2",
                                    "kull-" = "darkolivegreen2",
                                    "mle+" = "aquamarine",
                                    "mle-" = "aquamarine",
                                    "sq_cdf+" = "cyan3",
                                    "sq_cdf-" = "cyan3",
                                    "sq+" = "blue2",
                                    "sq-" = "blue2",
                                    "wass_cdf+" = "darkorchid1",
                                    "wass_cdf-" = "darkorchid1",
                                    "wass+" = "deeppink",
                                    "wass-" = "deeppink"
      )) +
      labs(title = title) +  # Add title
      theme(plot.title = element_text(hjust = 0.5)) # Center title
    
    # Show plot
    print(theplot)
    
  } else {
    #--------------------------------------------------------
    # 3x3 grid with one plot per method with -ve and +ve fits
    #--------------------------------------------------------
    
    # Print plot info in console
    cat(grid_plots_info,"\n")
    
    # Remove the initial estimate models
    df_all <- df[3:nrow(df),]
    # Define unique code to separate + from - models
    df_all$code_unique <- paste0(df_all$code,df_all$par_space)
    
    # Create empty list to store densities in
    dens_list <- vector("list", nrow(df_all))
    
    names(dens_list) <- df_all$code_unique
    
    # Extract the 4 parameters
    for (i in 1:nrow(df_all)) {
      parvec <- c(
        (df_all[i,"a"] %>% pull()),
        (df_all[i,"l"] %>% pull()),
        (df_all[i,"c"] %>% pull()),
        (df_all[i,"mu"]) %>% pull())
      
      # Calculate density for current method
      dens_obj <- dgg4(xx,
                       a = parvec[1],
                       l = parvec[2],
                       c = parvec[3],
                       mu = parvec[4])
      
      # Save density for method in a new object and assign to list
      index_method <- df_all$code_unique[i]
      dens_list[[index_method]] <- dens_obj
    }
    
    # Put densities of all methods in a dataframe
    dens_df <- data.frame(matrix(NA, ncol = length(df_all$code_unique), nrow = 10000))
    colnames(dens_df) <- df_all$code_unique
    for (i in 1:ncol(dens_df)) {
      dens_df[,i] <- dens_list[[i]]
    }
    
    # Add x values
    dens_df$x <- xx
    
    # Transform to long format for plotting
    dens_long <- dens_df %>%
      pivot_longer(cols = -x, names_to = "Method", values_to = "Density")
    
    # Create empty list to store plots in
    plots <- list()
    
    # Create string that filters out both par spaces for each method
    methnames <- unique(df$code)[-1]
    # Titles
    methodtitles <- unique(df$method)[-1]
    
    # Create 9 plots (1 per method)
    for (i in 1:length(methnames)) {
      # Create codes for both par spaces
      codes <- c(paste0(methnames[i],"+"),paste0(methnames[i],"-"))
      # Filter data for current method
      methdat <- dens_long %>% filter(Method == codes)
      # Make one plot for each method
      p <- ggplot() +
        # Create histogram
        geom_histogram(aes(x = y, y = ..density..), bins = 30, fill = "grey77", color = "grey40") +
        # Add standard nonparametric kernel density estimate
        geom_density(aes(x=y, colour = "R NP Kernel Density")) +
        # Add Amoroso fit
        geom_line(data = methdat, aes(x = x, y = Density, color = Method), size = 1.2) +
        # Add title
        labs(title = methodtitles[i]) +  # Add title
        theme(plot.title = element_text(hjust = 0.5)) # Center title
      
      # Add the plot to the list
      plots[[length(plots) + 1]] <- p
    }
    
    # Combine plots into a grid
    grid_plot <- do.call("grid.arrange", c(plots, ncol = 3))
    
    # Show plots
    print(grid_plot)
    
  }
  
  # Print best model
  win_mod <- df %>% filter(win == 1)
  win_mod_name <- paste0(win_mod$method, win_mod$par_space)
  cat("Model with the lowest BIC: \n", win_mod_name, "\n")
  cat(" with BIC:", win_mod$BIC, "\n\n")
  
  options(scipen = 999)
  # Return dataframe with parameter estimates and BIC for each model
  return(df)
  
}


#--------------------------------------------------------------------------

#-------------------
# Test the function
#-------------------
#estimateAmoroso_MLE_MDE(dataframe = palmerpenguins::penguins, variable = "bill_length_mm")
#estimateAmoroso_MLE_MDE(dataframe = palmerpenguins::penguins, variable = "bill_length_mm",)
#estimateAmoroso_MLE_MDE(vec = penguins$flipper_length_mm)
#estimateAmoroso_MLE_MDE(vec = penguins$bill_depth_mm, varname = "Penguin Bill Depth", plotpermethod =  F)

#-------------------
# Save the function
#-------------------
