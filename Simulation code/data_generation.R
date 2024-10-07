#####################File info#####################
# Code that stores the logic for data generation 

# The functions stored is:
#     1. generate_univariate_mixture_gaussian_data
#     2. generate_multivariate_mixture_gaussian_data
#     3. plot_univariate_gaussian_data
#     4. TODO: plot_bivariate_gaussian_data

#libraries
library(docstring)

generate_univariate_mixture_gaussian_data <- function(means, stds, probs, n){
    #' Univariate mixture of Gaussians data
    #' 
    #' @description
    #' The function generates mixtures of Gaussians (univariate) data. 
    #' 
    #'@param means vector. A vector of means, a mean for each component.
    #'@param sds vector. A vector of standard deviations, a sd for each component.
    #'@param probs vector. A vector of component mixing probabilities. 
    #'
    #'@return data vector. A vector storing the data (matrix if multivariate).
    
    ## Ensure different parameter lengths match
    if(length(means) != length(stds) || length(means) != length(probs)){
        stop("Number of means, standard deviations and components probabilities 
             must be equal.")
    }
    data <- matrix(nrow = n)
    for(i in 1:n){
        component <- which(rmultinom(1, 1, probs) == 1)
        data[i] <- rnorm(1, mean = means[component], sd = stds[component])
    }
    return(data)
}

plot_univariate_gaussian_data <- function(data, n_components, params, 
                                          overlay = T, freq = F){
    #' Plot univariate mixture of Gaussians data
    #' 
    #' @description
    #' The function plots univariate mixture of Gaussians data.   
    #' 
    #'@param data vector (if univariate). Vector storing data points.
    #'@param n_components integer. The number of components contained in data.
    #'@param params list. List of model parameters (mean, sd, mixing_probs).
    #'@param overlay bool (Default = True). Overlay component densities.
    #'@param freq bool (Default = False). If TRUE then histogram shows frequencies.
    
    ## Extracting parameters
    means <- params$means
    stds <- params$stds
    probs <- params$mixing_probs
    ## Histogram
    hist(data, breaks = 50, col = "lightblue", xlab = "Value", 
         freq = freq)
    ## Overlayed line plot
    if (overlay) {
        colors <- rainbow(n_components)
        for (i in 1:n_components){
            curve(dnorm(x, mean = means[i], sd = stds[i])*probs[i], add = TRUE, 
                  col = colors[i], lwd = 2)
        }
    }
}

library(ggplot2)

# Combined function to calculate prediction and plot GMM over histogram
plot_gmm_predictions <- function(data_nodes, params, bins = 50, hist_fill = "lightblue", line_color = "red") {
  
  # Function to calculate the predicted density
  calculate_prediction <- function(data, params) {
    n_components <- length(params$means)
    prediction <- rep(0, length(data))
    
    for (i in 1:n_components) {
      prediction <- prediction + (params$mixing_probs[i] * dnorm(data, 
                                                                 mean = params$means[i], 
                                                                 sd = params$stds[i]))
    }
    return(prediction)
  }
  
  # Calculate the predicted densities
  densities <- calculate_prediction(unlist(data_nodes), params = params)
  
  # Create a data frame for ggplot
  data_frame <- data.frame(
    x = unlist(data_nodes),
    density = densities
  )
  
  # Plot histogram with density estimation line
  ggplot(data_frame, aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = bins, fill = hist_fill, color = "black", alpha = 0.7) +
    geom_line(aes(y = density), color = line_color, size = 1) +
    labs(
      title = "Histogram and Density Estimation of GMM",
      x = "Value",
      y = "Density"
    ) +
    theme_minimal()
}

# Example usage
# plot_gmm_predictions(data_nodes, results$parameters)

