# Compare models 
file_path <- '/Users/bronwynmccall/Documents/GitHub/masters_research_code/Simulation code/'
source(paste(file_path, 'IEM_model.R', sep = ""))
source(paste(file_path, 'parallel_model.R', sep = ""))
source(paste(file_path, 'centralised_model.R', sep = ""))
library(psych)


# parameters
# SET MODEL HYPERPARAMETERS AND SETTINGS
K <- 3 # number of components
n <- 1000 # sample size
m <- 4 # number of nodes
n_m <- n/m # sample size per node

# TRUE PARAMETERS
pi_1 <- 0.3
pi_2 <- 0.2
pi_3 <- 0.5
mu_1 <- 0
mu_2 <- 5
mu_3 <- 2
si_1 <- 1
si_2 <- 0.3
si_3 <- 3
probs <- c(pi_1, pi_2, pi_3)
means <- c(mu_1, mu_2, mu_3)
stds <- c(si_1, si_2, si_3)
params <- list("mixing_probs" = probs, "means" = means, "stds" = stds)
params
pi_1_init <- 0.4
pi_2_init <- 0.2
pi_3_init <- 0.4
mu_1_init <- 0.2
mu_2_init <- 4
mu_3_init <- 1
si_1_init <- 1
si_2_init <- 0.5
si_3_init <- 2
probs_init <- c(pi_1_init, pi_2_init, pi_3_init)
means_init <- c(mu_1_init, mu_2_init, mu_3_init)
stds_init <- c(si_1_init, si_2_init, si_3_init)
params_init <- list("mixing_probs" = probs_init,
                    "means" = means_init,
                    "stds" = stds_init)

data_chunks <- lapply(1:m, function(x) {
    generate_univariate_mixture_gaussian_data(means = means,
                                              stds = stds,
                                              probs = probs,
                                              n = n_m)
})
params
plot_univariate_gaussian_data(unlist(data_chunks), 3, params)
plot_gmm_predictions(data_chunks, params)
run_and_measure_IEM <- function(data_chunks, params_init){
    # Get initial global and local statistics
    initial_local_statistics <- lapply(1:m, function(x) calculate_node_local_sufficient_statistics(
        data_chunks[[x]], params = params_init
    ))
    initial_global_statistics <- list('global_suff_stat_one' = 0,
                                 'global_suff_stat_two' = 0,
                                 'global_suff_stat_three' = 0)
    
    for(local_sufficient_stat in initial_local_statistics){
        initial_global_statistics$global_suff_stat_one <- initial_global_statistics$global_suff_stat_one + local_sufficient_stat$suff_stat_one
        initial_global_statistics$global_suff_stat_two <- initial_global_statistics$global_suff_stat_two + local_sufficient_stat$suff_stat_two
        initial_global_statistics$global_suff_stat_three <- initial_global_statistics$global_suff_stat_three + local_sufficient_stat$suff_stat_three
    }
    start_time <- Sys.time()
    iem_distributed_model_results <- iem_distributed_model(data_chunks,
                                                           initial_local_statistics,
                                                           initial_global_statistics,
                                                           params_init)
    total_time <- Sys.time() - start_time
    return(list('model_results' = iem_distributed_model_results,
                'time_taken' = total_time))
}
incremental_test_results <- run_and_measure_IEM(data_chunks,
                                                params_init)
incremental_test_results
run_and_measure_parallel <- function(data_chunks, 
                                     params_init){
    start_time <- Sys.time()
    parallel_distributed_model_results <- GMM_using_parallel(data_chunks,
                                                             params_init)
    total_time <- Sys.time() - start_time
    return(list('model_results' = parallel_distributed_model_results,
                'time_taken' = total_time))
}
parallel_test_results <- run_and_measure_parallel(data_chunks,
                                                  params_init)
parallel_test_results
run_and_measure_centralised <- function(data_chunks,
                                        params_init){
    start_time <- Sys.time()
    data <- NULL
    for (i in 1:length(data_chunks)){
      data <- rbind(data, data_chunks[[i]])
    }
    centralised_model_results <- centralised_model(unlist(data_chunks),
                                                   params_init)
    total_time <- Sys.time() - start_time
    return(list('model_results' = centralised_model_results,
                'time_taken' = total_time))
}
centralised_test_results <- run_and_measure_centralised(data_chunks,
                                                        params_init)
centralised_test_results

# Initialize storage for the results of each model
incremental_means <- c()
incremental_stds <- c()
incremental_mixing_probs <- c()
incremental_likelihood <- c()
incremental_time <- c()

parallel_means <- c()
parallel_stds <- c()
parallel_mixing_probs <- c()
parallel_likelihood <- c()
parallel_time <- c()

centralised_means <- c()
centralised_stds <- c()
centralised_mixing_probs <- c()
centralised_likelihood <- c()
centralised_time <- c()

# Number of iterations
iterations <- 1000

# Main loop for 1000 iterations
for (i in 1:iterations) {
    # Generate data chunks for each iteration
    data_chunks <- lapply(1:m, function(x) {
        generate_univariate_mixture_gaussian_data(means = means,
                                                  stds = stds,
                                                  probs = probs,
                                                  n = n_m)
    })
    
    # Incremental EM
    IEM_result <- run_and_measure_IEM(data_chunks, params_init)
    incremental_means <- rbind(incremental_means, IEM_result$model_results$params$means)
    incremental_stds <- rbind(incremental_stds, IEM_result$model_results$params$stds)
    incremental_mixing_probs <- rbind(incremental_mixing_probs, IEM_result$model_results$params$mixing_probs)
    incremental_likelihood <- rbind(incremental_likelihood,
                                    calculate_log_likelihood_with_full_data(unlist(data_chunks), IEM_result$model_results$params))
    incremental_time <- rbind(incremental_time, as.numeric(IEM_result$time_taken))
    
    # Parallel EM
    parallel_result <- run_and_measure_parallel(data_chunks, params_init)
    parallel_means <- rbind(parallel_means, parallel_result$model_results$params$means)
    parallel_stds <- rbind(parallel_stds, parallel_result$model_results$params$stds)
    parallel_mixing_probs <- rbind(parallel_mixing_probs, parallel_result$model_results$params$mixing_probs)
    parallel_likelihood <- rbind(parallel_likelihood,
                                    calculate_log_likelihood_with_full_data(unlist(data_chunks), parallel_result$model_results$params))
    parallel_time <- rbind(parallel_time, as.numeric(parallel_result$time_taken))
    
    # Centralised EM
    centralised_result <- run_and_measure_centralised(unlist(data_chunks), params_init)
    centralised_means <- rbind(centralised_means, centralised_result$model_results$params$means)
    centralised_stds <- rbind(centralised_stds, centralised_result$model_results$params$stds)
    centralised_mixing_probs <- rbind(centralised_mixing_probs, centralised_result$model_results$params$mixing_probs)
    centralised_likelihood <- rbind(centralised_likelihood,
                                    calculate_log_likelihood_with_full_data(unlist(data_chunks), centralised_result$model_results$params))
    centralised_time <- rbind(centralised_time, as.numeric(centralised_result$time_taken))
}

# Plotting histograms for means, stds, and mixing probabilities for all models
par(mar = c(4, 4, 2, 1))  # Set smaller margins
# Plot histograms for means
lapply(1:K, function(i) {
    hist(incremental_means[, i], main = paste("Incremental Means (Component", i, ")"), xlab = "Means")
    hist(parallel_means[, i], main = paste("Parallel Means (Component", i, ")"), xlab = "Means")
    hist(centralised_means[, i], main = paste("Centralised Means (Component", i, ")"), xlab = "Means")
})

# Set up a 3x1 grid for subplots (3 rows, 1 column for each component)
lapply(1:K, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns
  
  # Plot incremental means
  hist(incremental_means[, i], 
       main = paste("Incremental Means (Component", i, ")"), 
       xlab = "Means", 
       col = "lightblue")
  
  # Plot parallel means
  hist(parallel_means[, i], 
       main = paste("Parallel Means (Component", i, ")"), 
       xlab = "Means", 
       col = "lightgreen")
  
  # Plot centralised means
  hist(centralised_means[, i], 
       main = paste("Centralised Means (Component", i, ")"), 
       xlab = "Means", 
       col = "lightpink")
})

# Reset plot layout to default (1 plot per window)
par(mfrow = c(1, 1))

# Set up a 3x1 grid for subplots (3 rows, 1 column for each component)
lapply(1:K, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns
  
  # Plot incremental means
  hist(incremental_stds[, i], 
       main = paste("Incremental stds (Component", i, ")"), 
       xlab = "Means", 
       col = "lightblue")
  
  # Plot parallel means
  hist(parallel_stds[, i], 
       main = paste("Parallel stds (Component", i, ")"), 
       xlab = "Means", 
       col = "lightgreen")
  
  # Plot centralised means
  hist(centralised_stds[, i], 
       main = paste("Centralised stds (Component", i, ")"), 
       xlab = "Means", 
       col = "lightpink")
})

# Reset plot layout to default (1 plot per window)
par(mfrow = c(1, 1))

# Set up a 3x1 grid for subplots (3 rows, 1 column for each component)
lapply(1:K, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns
  
  # Plot incremental means
  hist(incremental_mixing_probs[, i], 
       main = paste("Incremental mixing probabilities (Component", i, ")"), 
       xlab = "Means", 
       col = "lightblue")
  
  # Plot parallel means
  hist(parallel_mixing_probs[, i], 
       main = paste("Parallel mixing probabilities (Component", i, ")"), 
       xlab = "Means", 
       col = "lightgreen")
  
  # Plot centralised means
  hist(centralised_mixing_probs[, i], 
       main = paste("Centralised mixing probabilities (Component", i, ")"), 
       xlab = "Means", 
       col = "lightpink")
})

# Reset plot layout to default (1 plot per window)
par(mfrow = c(1, 1))

# Plot histograms for standard deviations
lapply(1:K, function(i) {
    hist(incremental_stds[, i], main = paste("Incremental STDs (Component", i, ")"), xlab = "Standard Deviations")
    hist(parallel_stds[, i], main = paste("Parallel STDs (Component", i, ")"), xlab = "Standard Deviations")
    hist(centralised_stds[, i], main = paste("Centralised STDs (Component", i, ")"), xlab = "Standard Deviations")
})

# Plot histograms for mixing probabilities
lapply(1:K, function(i) {
    hist(incremental_mixing_probs[, i], main = paste("Incremental Mixing Probs (Component", i, ")"), xlab = "Mixing Probabilities")
    hist(parallel_mixing_probs[, i], main = paste("Parallel Mixing Probs (Component", i, ")"), xlab = "Mixing Probabilities")
    hist(centralised_mixing_probs[, i], main = paste("Centralised Mixing Probs (Component", i, ")"), xlab = "Mixing Probabilities")
})

# Display time taken for each model
hist(incremental_time)
hist(parallel_time)
hist(centralised_time)

describe(incremental_means)
describe(parallel_means)
describe(centralised_means)

describe(incremental_stds)
describe(parallel_stds)
describe(centralised_stds)

describe(incremental_mixing_probs)
describe(parallel_mixing_probs)
describe(centralised_mixing_probs)


