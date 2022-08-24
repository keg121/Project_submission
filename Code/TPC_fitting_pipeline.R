#!/usr/bin/env Rscript

library(minpack.lm)
library(nls.multstart)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)

rm(list=ls())
graphics.off()

#####################
# F U N C T I O N S #
#####################

fit_Sharpe_Schoolfield <- function(dataset)
{
  
  # Set the Boltzmann constant.
  k <- 8.617 * 10^-5
  
  # Set the minimum trait measurement at the rise of the TPC as the 
  # starting value for B_0.
  B_0_start <- min(dataset$trait_value[dataset$temp == min(dataset$temp)])
  
  # Set the starting value of T_pk to the temperature at which the 
  # trait reaches its maximum value.
  T_pk_start <- max(dataset$temp[dataset$trait_value == max(dataset$trait_value)])
  
  # If the dataset has measurements before the thermal optimum, 
  # a starting value for E can be set as roughly the slope of the 
  # following regression for that subset of the data:
  #
  # ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
  #
  # Otherwise, just set E to 0.6 eV.
  
  dataset_before_peak <- dataset[dataset$temp < T_pk_start,]
  if ( 
    nrow(dataset_before_peak) < 2 || 
    length(unique(dataset_before_peak$temp)) < 2 || 
    length(unique(dataset_before_peak$trait_value)) < 2 
  )
  {
    E_start <- 0.6
  } else
  {
    y_vals <- log(dataset_before_peak$trait_value)
    x_vals <- 1/(k * (dataset_before_peak$temp + 273.15))
    
    # Take the absolute value.
    E_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
    
    if ( E_start >= 10 )
    {
      E_start <- 0.6
    }
  }
  
  # If the dataset has measurements after the thermal optimum, 
  # a starting value for E_D can be set as the slope of the following 
  # regression for that subset of the data:
  #
  # ln(trait_value) ~ intercept + slope * 1/(k * (T + 273.15))
  #
  # Otherwise, just set E_D to 3 eV.
  dataset_after_peak <- dataset[dataset$temp > T_pk_start,]
  if ( 
    nrow(dataset_after_peak) < 3 || 
    length(unique(dataset_after_peak$temp)) < 3 || 
    length(unique(dataset_after_peak$trait_value)) < 3 
  )
  {
    E_D_start <- 3
  } else
  {
    y_vals <- log(dataset_after_peak$trait_value)
    x_vals <- 1/(k * (dataset_after_peak$temp + 273.15))
    
    # Take the absolute value.
    E_D_start <- abs(lm(y_vals ~ x_vals)$coefficients[2])
    
    if ( E_D_start >= 50 )
    {
      E_D_start <- 3
    }
  }
  
  if ( E_start >= E_D_start )
  {
    E_start <- 0.9 * E_D_start
  }
  
  function_to_be_fitted <- function(B_0, E, T_pk, E_D, temp)
  {
    temp <- temp + 273.15
    
    # Set the Boltzmann constant.
    k <- 8.617 * 10^-5
    
    if ( E == E_D )
    {
      return(rep(1e10, length(temp)))
    } else
    {
      return(
        log(
          B_0 * temp * exp(-E * ((1/(k*temp)) - (1/(k*273.15)))) / ( 1 + ( E/( E_D - E ) ) * exp( (E_D/k) * (1/T_pk - 1/temp) ) )
        )
      )
    }
  }
  
  fit <- NULL
  
  set.seed(1)
  try(
    fit <- nls_multstart(
      log(trait_value) ~ function_to_be_fitted(
        B_0, E, T_pk, E_D, temp = temp
      ),
      data = dataset,
      iter = 1000,
      start_lower = c(
        B_0 = 0.5 * B_0_start,				E = 0.5 * E_start,
        T_pk = 0.5 * T_pk_start + 273.15,	E_D = 0.5 * E_D_start
      ),
      start_upper = c(
        B_0 = 1.5 * B_0_start,				E = 1.5 * E_start,
        T_pk = 1.5 * T_pk_start + 273.15,	E_D = 1.5 * E_D_start
      ),
      supp_errors = 'Y',
      convergence_count = FALSE,
      control = nls.lm.control(
        ftol = .Machine$double.eps, ptol = .Machine$double.eps, maxiter = 1024, 
        maxfev = 100000
      ),
      lower = c(0, 0, 273.15, 0),
      upper = c(Inf, 10, 273.15 + 150, 50)
    )
  )
  
  trait_values_string <- toString(dataset$trait_value)
  temps_string <- toString(dataset$temp)
  
  if ( is.null(fit) )
  {
    B_0_estimate <- NA
    E_estimate <- NA
    E_stderr_estimate <- NA
    T_pk_estimate <- NA
    E_D_estimate <- NA
    r_squared <- NA
    data_points_up_to_peak <- NA
  } else
  {
    
    B_0_estimate <- summary(fit)$parameters[1,1]
    E_estimate <- summary(fit)$parameters[2,1]
    E_stderr_estimate <- summary(fit)$parameters[2,2]
    T_pk_estimate <- summary(fit)$parameters[3,1]
    E_D_estimate <- summary(fit)$parameters[4,1]
    
    # Calculate the R-squared.
    rss <- sum((dataset$trait_value - exp(fitted(fit)))^2)
    tss <- sum((dataset$trait_value - mean(dataset$trait_value))^2)
    
    r_squared <- 1 - (rss/tss)
    
    data_points_up_to_peak <- nrow(
      dataset[dataset$temp < (T_pk_estimate - 273.15),]
    )
  }
  
  return(
    c(
      dataset$originalid[1],	dataset$interactor1[1],
      dataset$latitude[1],		dataset$longitude[1],
      trait_values_string,		dataset$originaltraitunit[1],
      temps_string,				    dataset$interactor1stage[1],
      dataset$interactor1sex[1],	dataset$citation[1],
      r_squared,					B_0_estimate,
      E_estimate,					E_stderr_estimate,
      T_pk_estimate - 273.15,		E_D_estimate,
      data_points_up_to_peak
    )
  )
}

# This function splits the data into subsets with unique IDs.
split_data <- function(dataset)
{
  split_dataset <- list()
  
  # Remove trait values <= 0 or NA.
  dataset$originaltraitvalue <- as.numeric(dataset$originaltraitvalue)
  dataset <- dataset[!is.na(dataset$originaltraitvalue) & dataset$originaltraitvalue > 0,]
  
  unique_IDs <- unique(dataset$originalid)
  
  for ( i in unique_IDs )
  {
    temp_dataset <- dataset[dataset$originalid == i,]
    
    temp_dataset$temp <- temp_dataset$interactor1temp
    temp_dataset$trait_value <- temp_dataset$originaltraitvalue
    
    # Exclude experimental TPCs with less than 4 unique temperatures.
    if ( length(unique(temp_dataset$temp)) < 3 )
    {
      next
    }
    
    # Create a (simplified) dataset with useful information only.
    split_dataset[[as.character(i)]] <- temp_dataset[,
                                                     c(
                                                       'originalid',			'interactor1',
                                                       'originaltraitname','originaltraitdef',
                                                       'latitude',				'longitude',
                                                       'trait_value',			'originaltraitunit',
                                                       'temp',					  'interactor1stage',
                                                       'interactor1sex',		'citation'
                                                     )
    ]
  }
  
  return(split_dataset)
}

#####################################
# M A I N  C O D E  -- All MR types #
#####################################

# Read the input dataset and filter out species with single data points

all_data <- as_tibble(read.csv('MRdata.csv', stringsAsFactors = FALSE)) %>%
  filter(interactor1 != 'Anax junius' &
         interactor1 != 'Brachymesia gravida' &
         interactor1 !='Conocephalus fasciatus' &
         interactor1 !='Erythemis simplicicollis' &
         interactor1 !='Erythrodiplax berenice' &
         interactor1 !='Formica exsecta' &
         interactor1 !='Lasius alienus' &
         interactor1 !='Lasius flavus' &
         interactor1 !='Ischnura elegans' &
         interactor1 !='Libellula auripenms' &
         interactor1 !='Libellula needhami' &
         interactor1 !='Pantala flavescens' &
         interactor1 !='Perithemis tenera' &
         interactor1 !='Tramea carolina' & 
         originaltraitname != 'body size')

# visualise data
MR_all <- ggplot(all_data) +
  geom_point(aes(interactor1temp, originaltraitvalue, col = interactor1)) +
  facet_wrap(~interactor1, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none')

ggsave("../Desktop/AllMRs.pdf", MR_all, width = 40, height = 80, units = "cm")


# Separate the data by ID.
split_dataset <- split_data(all_data)

results <- data.frame(
  ID = rep(NA, length(split_dataset)),
  Species = rep(NA, length(split_dataset)),
  latitude = rep(NA, length(split_dataset)),
  longitude = rep(NA, length(split_dataset)),
  trait_values = rep(NA, length(split_dataset)),
  trait_unit = rep(NA, length(split_dataset)),
  temps = rep(NA, length(split_dataset)),
  stage = rep(NA, length(split_dataset)),
  sex = rep(NA, length(split_dataset)),
  citation = rep(NA, length(split_dataset)),
  R_squared = rep(NA, length(split_dataset)),
  B_0 = rep(NA, length(split_dataset)),
  E = rep(NA, length(split_dataset)),
  E_stderr = rep(NA, length(split_dataset)),
  T_pk = rep(NA, length(split_dataset)),
  E_D = rep(NA, length(split_dataset)),
  data_points_up_to_peak = rep(NA, length(split_dataset))
)

# Fit the models to each experimental TPC.
for ( i in 1:length(split_dataset) )
{
  cat("Now at TPC ", i, "/", length(split_dataset), " ...\n", sep = '')
  
  results[i,] <- fit_Sharpe_Schoolfield(split_dataset[[i]])	
}

# remove r square values < 0.4
results <- results %>% filter(R_squared > 0.4)

write.csv(results, file = 'TPC_parameter_estimates.csv', row.names = FALSE)


