###############################
# Phylogenetic analysis script#
##############################
# 1. Load in/wrangle data  & open packages

# ~E~
# 2A. Trees for E
# 2B. std error for E
# 2C. Phylosig: E

# ~Latitude vs E~
# 3A. Model 1: Linear regression
# 3B. Model 2: second order polynomial
# 3C. Compare models

# ~Body Size vs E~
# 4A. Model 1: Linear regression
# 4B. Model 2: second order polynomial
# 4C. Compare models

# ~Body Size and latitude vs E~

# 6. Plotting

#########################
#1. Load in/wrangle data#
#########################
library(phytools)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(gridExtra)

# load data
data <- read.csv("TPC_parameter_estimates.csv", row.names=1)

## data wrangling##

# Filter out any larval groups
#data <- data[!grepl("larval", data$stage),]

# remove r square values < 0.5
data <- data %>% filter(R_squared > 0.5)

# remove E values < 2.5
data <- data %>% filter(E < 2.5)

# remove values where there are less than 3 observations up to the peak
data <- data %>% filter(data_points_up_to_peak > 3)

# load e estimates as vector
x<-as.numeric(as.matrix(data)[,12]) 
names(x) <- data[,1] # add names of species to  corresponding E value

y <- x

## note###
##########
# Below is removed as no longer necessary

# y <- x[names(x) != "Hippodamia convergens"] # removes observations where name =/= H. convergens
# 
# y[length(y) + 1] <- mean(x[names(x) == "Hippodamia convergens"]) # mean of both H. convergens E values
# 
# # add name
# names(y)[length(y)] <- 'Hippodamia convergens'
# # y is nw a vctor with e values and corresponding species names

########
# 2. E #
########

#################
# 2A. Trees for E
#################
# load tree
tree<-read.tree("final_calibrated_tree.phy")

# fix typos
tree$tip.label <- gsub('_', ' ', tree$tip.label) # This line changes underscores to spaces in species' names.

# Names that had to be replaced with synonyms. These are changed back
tree$tip.label[tree$tip.label == 'Hogna carolinensis'] <- 'Lycosa carolinensis' 
tree$tip.label[tree$tip.label == 'Novomessor cockerelli'] <- 'Aphaenogaster cockerelli'
tree$tip.label[tree$tip.label == 'Bootettix argentatus'] <- 'Bootettix punctatus'
tree$tip.label[tree$tip.label == 'Trimerotropis verruculata'] <- 'Trimerotropis suffusa'

# species only incl in the phylogeny 
z <- y[names(y) %in% tree$tip.label]

# prune tree
tree_pruned <- keep.tip(tree, names(z))

# make trees
dot <- dotTree(tree_pruned,z,length=10,ftype= "i", fsize=0.7,)

###############
# 2B. std error
###############
# load std error as a vector.
x_error <- as.numeric(as.matrix(data)[,13]) # load in std errors
names(x_error) <- data[,1] # assign names

y_error <- x_error

# get species only included in the phylogeny 
z_error <- y_error[names(y_error) %in% tree$tip.label]

#################
# 2C. Phylosig: E
#################

# Initialise empty vectors for model output
lambda <- as.vector(matrix(nrow=1000, ncol=1))
logL <- as.vector(matrix(nrow=1000, ncol=1))
MLE <- as.vector(matrix(nrow=1000, ncol=1))
P_val <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

# Run phylosig 1000 times and save output
for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, z, method="lambda", test=TRUE, nsim=1000, se=z_error, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL[i] <- tmp[[3]]
  
  # MLE 
  MLE[i] <- tmp[[2]]
  
  # P-values
  P_val[i] <- tmp[[7]]
}

# make a df for results
results <- matrix(NA, ncol = 4, nrow = 1000)
results <- data.frame(results)
# plot(density(results$lambda)) 

# rename cols
names(results) <- c("lambda", "logL", "MLE", "P_value")

# Add estimates to DF
as.numeric(results$lambda <- lambda)
results$logL <- logL
results$MLE <- MLE
results$P_value <- P_val

# select row with highest Log liklihood value
best_fit <- results[which.max(results$logL),]

# estimating 95% CI
#calulate 95th percintiles 
# manually limit density plot to only include 95% ci 
# df = . %>% filer(x<= (calulate where 2.5% is in data))
sd_ss <- sd(results$lambda)
mean_ss <- mean(results$lambda)
n_ss <- length(results$lambda)
se_ss <- sd_ss/sqrt(n_ss)
t_val_ss <- qt(0.975, df = 999)
x_ss <- se_ss*t_val_ss
upper_ss <- mean_ss+x_ss
lower_ss <- mean_ss-x_ss

###################
# 3. E vs latitude#
###################
# select unique values for latitude
lat_vec <- c()
for (i in 1:length(z)){
  lat_vec[i] <- unique(data$latitude[data$Species == names(z)[i]])
} # this line takes all the unique values for latitude of species

# make df with latitude, e estimates and error
lat_df <- data.frame(Species = names(z), latitude = lat_vec, E = z, error = z_error) 
lat_df <- lat_df %>% drop_na() # remove NAs

##################################
# 3. A: Model 1: Linear regression
##################################

# fit simple linear model with e as response var and abs latitude as predictor 
# e values, latitude and errors in dataframe. have row names as species names 

# linear model
lat_model1 <- lm(E ~ abs(latitude), data= lat_df,  weights = 1/(error^2))
# results suggest higher latitudes = smaller E

# Get residuals for model 1 
res_lat1 <- resid(lat_model1)

# get fitted values
fitted_lat1 <- fitted(lat_model1)

# Plot fitted v residuals: 
plot(fitted_lat1, res_lat1, pch=21, col="black", bg="rosybrown3",
     main = "Linear Regression (E ~ absolute latitude): Residuals vs fitted", cex = 1.2, cex.main= 0.9, xlab="Fitted value", ylab="Residual")

# get model predictions for latitudes 
#names(summary(lat_model1))
lat_intercept1 <- summary(lat_model1)$coefficients[1,1]
lat_slope1 <- summary(lat_model1)$coefficients[2,1]

# plot  this instead of the fitted values 
lat_predict1 <- data.frame(Latitude = seq(-50, 80, 1), E = lat_intercept1 + lat_slope1 * abs(seq(-50, 80, 1)))

# Lambda estimates
##################
# Initialise empty vectors for model output
lambda_lat <- as.vector(matrix(nrow=1000, ncol=1))
logL_lat <- as.vector(matrix(nrow=1000, ncol=1))
MLE_lat <- as.vector(matrix(nrow=1000, ncol=1))
P_val_lat <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, res_lat1, method="lambda", test=TRUE, nsim=1000, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda_lat[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL_lat[i] <- tmp[[2]]
  
  # P-values
  P_val_lat[i] <- tmp[[4]]
}

# make a df for results
results_lat <- data.frame(lambda = lambda_lat, logL = logL_lat, P_value = P_val_lat)

# select row with highest Log liklihood value
best_fit_lat <- results_lat[which.max(results_lat$logL),]
## results say little signal (results are true)

########################################
## checking the CI interval of residuals 
########################################
# modified confint
##################
# confint.double <- function(x, conf = 0.95, method = c("normal", "quantile")) {
#   method <- match.arg(method)
#   l <- (1 - conf)/2
#   l <- c(l, 1 - l)
#   nms <- paste0(round(100*l, 2), "%")
#   if(method == "normal") {
#     xbar <- mean(x)
#     se <- sd(x)
#     ci <- xbar + qnorm(l)*se
#     setNames(ci, nms)
#   } else {
#     quantile(x, l)
#   }
# }
# 
# confint.double(res_lat1, method = "q")
# 

######################################
# 3B. Model 2: second order polynomial
######################################

# fit model
lat_model2<- lm(E ~ poly(latitude, 2, raw = TRUE), weights = 1/(error^2), data=lat_df)

# get residuals
res_lat2 <- resid(lat_model2)

# Get fitted values
fitted_lat2 <- fitted(lat_model2)

# Plot fitted v residuals: 
plot(fitted_lat2, res_lat2, pch=21, col="black", bg="rosybrown3",
     main = "Polynomial Regression (E ~ latitude): Residuals vs fitted", cex.main=0.9, cex= 1.2, xlab="Fitted value", ylab="Residual")

# get intercept and slope for model prediction for latitudes between -40 to 40
#names(summary(lat_model2))
lat_intercept2 <- summary(lat_model2)$coefficients[1,1]
lat_slope2_1 <- summary(lat_model2)$coefficients[2,1]
lat_slope2_2 <- summary(lat_model2)$coefficients[3,1]
# higher = more likely to be a difference 
# p value correponds to getting a t value that low i there isnt a sig difference
# first degree: lower p value so more significant 


# make a df with latitude and model predictions
# test <- data.frame(Latitude = seq(0, 40, 1), E = intercept + slope * seq(0, 40, 1))
# plot  this instead of the fitted values 

lat_predict2 <- data.frame(Latitude = seq(-50, 80, 1), E = lat_intercept2 + lat_slope2_1 * seq(-50, 80, 1) + lat_slope2_2 * seq(-50, 80, 1)^2)

# Lamba estimates
#################
# Initialise empty vectors for model output
lambda_lat2 <- as.vector(matrix(nrow=1000, ncol=1))
logL_lat2 <- as.vector(matrix(nrow=1000, ncol=1))
MLE_lat2 <- as.vector(matrix(nrow=1000, ncol=1))
P_val_lat2 <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, res_lat2, method="lambda", test=TRUE, nsim=1000, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda_lat2[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL_lat2[i] <- tmp[[2]]
  
  # P-values
  P_val_lat2[i] <- tmp[[4]]
}

# make a df for results
results_lat2 <- data.frame(lambda = lambda_lat2, logL = logL_lat2, P_value = P_val_lat2)

# select row with highest Log liklihood value
best_fit_lat2 <- results_lat2[which.max(results_lat2$logL),]
## results say little signal (results are true)

####################
# 3C. Compare models
####################

# plotting models
#################
compare_lat <- 
  ggplot(aes(x=latitude, y =E), data=lat_df) +
  ggtitle("Model comparison: Latitude") +
  geom_point() + 
  geom_line(aes(x= abs(Latitude), y=E, colour="Linear"), data = lat_predict1) + 
  geom_line(aes(x = Latitude, y=E, colour="Polynomial"), data = lat_predict2) + 
              scale_colour_manual("",  
                                  breaks = c("Linear", "Polynomial"), 
                                  values = c("red","blue"))

compare_lat_LINEAR <- 
  ggplot(aes(x=abs(latitude), y =E), data=lat_df) +
  ggtitle("Linear Regression") +
  geom_point() + 
  xlim(0,80) +
  labs(x="Absolute latitude") +
  theme_classic()+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  geom_line(aes(x= abs(Latitude), y=E), colour="red",data = lat_predict1) 

compare_lat_POLY <- 
  ggplot(aes(x=latitude, y =E), data=lat_df) +
  ggtitle("Polynomial Regression") +
  geom_point() + 
  labs(x="Latitude") +
  xlim(-50,80) +
  theme_classic()+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  geom_line(aes(x = Latitude, y=E), colour="Blue", data = lat_predict2) 

compare_lat_grid <-plot_grid(compare_lat_LINEAR, compare_lat_POLY, labels=c("A", "B"), ncol = 2, nrow = 1)
  
# AIC
#####
AIC_lat_model1 <-AIC(lat_model1)
AIC_lat_model2 <-AIC(lat_model2) # model 2: better AIC

## BIC
######
BIC_lat_model1 <-BIC(lat_model1)
BIC_lat_model2 <-BIC(lat_model2) # model 2: better BIC

# note, the second degree of the polynomial model had a p value of approx 0.4 compare to the linear model which had a p value 0f 0.02 so linear is better for latitue

####################
# 4. E vs body size#
####################
# load in body sizes 
data_size <- read.csv("MRdata.csv")

# convery dry mass to wet mass by multiplying by 3
for (i in 1:nrow(data_size)){
  if (str_detect(data_size$interactor1sizetype[i], "dry") == TRUE){
    data_size$interactor1size[i] <- data_size$interactor1size[i] *3
    data_size$interactor1sizetype[i] <- "mean wet mass"
  }}

# convert mg to g
for (i in 1:nrow(data_size)){
  if (data_size$interactor1sizeunit[i] == "mg"){
    data_size$interactor1size[i] <- data_size$interactor1size[i]/1000
    data_size$interactor1sizeunit[i] <- "g"
  }
}

# Get mean body size per species
mean_mass <- aggregate(interactor1size ~ interactor1, data=data_size, FUN=mean) # gives mean of size based on species

mean_mass$interactor1size <-log(mean_mass$interactor1size) # log transform

colnames(mean_mass)[1] <- "Species" # Rename interactor1 col 
colnames(mean_mass)[2] <- "mean_body_mass" # Rename interactor1 col 

size_df_a <- right_join(mean_mass, data) # Join parameter estimate and mean mass DFs

# make df with size, e estimates and error
size_df <- data.frame(species = size_df_a$Species, E = size_df_a$E, size = size_df_a$mean_body_mass, error = size_df_a$E_stderr, row.names = 1) 
size_df <- size_df[!is.na(size_df$size), ] # remove NAs

################################
# 4A. Model 1: Linear regression
################################
# fit simple linear model with e as response var and bodysize as predictor 
# e values, latitude and errors in dataframe. have row names as species names 

# linear model
size_model1 <- lm(E ~ size, weights = 1/(error^2), data= size_df)
# results suggest higher latitudes = smaller E

# Get residuals for model 1
res_size1 <- resid(size_model1)

# get fitted values for model 1
fitted_size1 <- fitted(size_model1)

# Plot fitted v residuals: 
plot(fitted_size1, res_size1, pch=21, cex=1.2, cex.main=0.9, col="black", bg="slategray3",main = "Linear Regression (E ~ Mean body mass): Residuals vs fitted", xlab="Fitted value", ylab="Residual")

# get intercept and slope for model prediction
size_intercept1 <- summary(size_model1)$coefficients[1,1]
size_slope1 <- summary(size_model1)$coefficients[2,1]

# make a df with size and model predictions
size_predict1 <- data.frame(Size = seq(-9, 3, 0.5), E = size_intercept1 + size_slope1 * seq(-9, 3, 0.5))


size_df <- size_df[!is.na(size_df$error), ] # remove NAs

# Initialise empty vectors for model output
lambda_size1 <- as.vector(matrix(nrow=1000, ncol=1))
logL_size1 <- as.vector(matrix(nrow=1000, ncol=1))
MLE_size1 <- as.vector(matrix(nrow=1000, ncol=1))
P_val_size1 <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, res_size1, method="lambda", test=TRUE, nsim=1000, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda_size1[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL_size1[i] <- tmp[[2]]
  
  # P-values
  P_val_size1[i] <- tmp[[4]]
}

# make a df for results
results_size1 <- data.frame(lambda = lambda_size1, logL = logL_size1, P_value = P_val_size1)

# select row with highest Log liklihood value
best_fit_size1 <- results_size1[which.max(results_size1$logL),]

######################################
# 4B. Model 2: second order polynomial
######################################

# Fit model
size_model2<- lm(E ~ poly(size, 2, raw = TRUE), weights = 1/(error^2), data=size_df)

# Get residuals
res_size2 <- resid(size_model2)

# Get fitted values
fitted_size2 <- fitted(size_model2)

# Plot fitted v residuals: 
plot(fitted_size2, res_size2, cex=1.2, cex.main=0.9,pch=21, col="black", bg="slategray3",
     main = "Polynomial Regression (E ~ Mean body size): Residuals vs fitted", xlab="Fitted value", ylab="Residual")

size_intercept2 <- summary(size_model2)$coefficients[1,1]
size_slope2_1 <- summary(size_model2)$coefficients[2,1]
size_slope2_2 <- summary(size_model2)$coefficients[3,1]

# make a df with size and model predictions
size_predict2 <- data.frame(Size = seq(-9, 3, 0.5), E = size_intercept2 + size_slope2_1 * seq(-9, 3, 0.5) + size_slope2_2 * seq(-9, 3, 0.5)^2)

# Lamba estimates
#################
# Initialise empty vectors for model output
lambda_size2 <- as.vector(matrix(nrow=1000, ncol=1))
logL_size2 <- as.vector(matrix(nrow=1000, ncol=1))
MLE_size2 <- as.vector(matrix(nrow=1000, ncol=1))
P_val_size2 <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, res_size2, method="lambda", test=TRUE, nsim=1000, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda_size2[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL_size2[i] <- tmp[[2]]
  
  # P-values
  P_val_size2[i] <- tmp[[4]]
}

# make a df for results
results_size2 <- data.frame(lambda = lambda_size2, logL = logL_size2, P_value = P_val_size2)

# select row with highest Log liklihood value
best_fit_size2 <- results_size2[which.max(results_size2$logL),]
## results say little signal (results are true)

####################
# 4C. Compare models
####################

# plotting models
#################
compare_size <-
ggplot(data = size_df, aes(x=size, y =E)) +
  ggtitle("Model comparison: Size") +
  labs(x = "log(size)") +
  geom_point() +
  geom_line(aes(x= Size, y = E, colour = "Linear"), data = size_predict1) +
  geom_line(aes(x = Size, y=E, colour="Polynomial"), data = size_predict2) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold")) +
  scale_colour_manual("", 
                      breaks = c("Linear", "Polynomial"),
                      values = c("red","blue")) 

# AIC
#####
AIC_size_model1 <-AIC(size_model1) # model 1 better
AIC_size_model2 <-AIC(size_model2) 

## BIC
######
BIC_size_model1 <-BIC(size_model1) # model 1 vetter
BIC_size_model2 <-BIC(size_model2) 

#######################
# 5. Mass*latitude vs E
#######################
# latitude * body_size: how does e change in a response to a combination of these two traits
# log transform body size  
# take mean of body mass per species 
# make df with species, E, latitude, mean body mass and error
size_lat_df <- data.frame(Species = size_df_a$Species, E = size_df_a$E, latitude = size_df_a$latitude, size = size_df_a$mean_body_mass, error = size_df_a$E_stderr, row.names = 1) 

size_lat_df <- size_lat_df %>% drop_na() # remove NAs

############################
# 5A. Model 1: linear model#
############################
# size_lat_model1 <- lm(E ~ size*abs(latitude), weights = 1/(error^2), data= size_lat_df)
size_lat_model1 <- lm(E ~ size*abs(latitude), weights = 1/(error^2), data= size_lat_df)

# Get residuals for model 1
res_size_lat1 <- resid(size_lat_model1)

# get fitted values for model 1
fitted_size_lat1 <- fitted(size_lat_model1)

# Plot fitted v residuals: 
plot(fitted_size_lat1, res_size_lat1, pch=21, col="black", bg="palegreen", cex=1.2, cex.main=0.9,
     main = "Linear Regression (E ~ size*absolute latitude): Residuals vs fitted", xlab="Fitted value", ylab="Residual")

# get intercept and slope for model prediction
size_lat_intercept1 <- summary(size_lat_model1)$coefficients[1,1]
size_lat_slope1 <- summary(size_lat_model1)$coefficients[4,1]

# make a df with size*lat and model predictions
size_lat_predict1 <- data.frame(Size_x_latitude = seq(-400, 200, 10), E = size_lat_intercept1 + size_lat_slope1 * seq(-400, 200, 10))

# Lambda esimates
#################
# Initialise empty vectors for model output
lambda_size_lat1 <- as.vector(matrix(nrow=1000, ncol=1))
logL_size_lat1 <- as.vector(matrix(nrow=1000, ncol=1))
P_val_size_lat1 <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, res_size_lat1, method="lambda", test=TRUE, nsim=1000, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda_size_lat1[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL_size_lat1[i] <- tmp[[2]]
  
  # P-values
  P_val_size_lat1[i] <- tmp[[4]]
}

# make a df for results
results_size_lat1 <- data.frame(lambda = lambda_size_lat1, logL = logL_size_lat1, P_value = P_val_size_lat1)

# select row with highest Log liklihood value
best_fit_size_lat1 <- results_size_lat1[which.max(results_size_lat1$logL),]

######################################
# 5B. Model 2: second order polynomial
######################################
# Fit model
size_lat_model2<- lm(E ~ poly(size*latitude, 2, raw = TRUE), weights = 1/(error^2), data=size_lat_df)

# Get residuals
res_size_lat2 <- resid(size_lat_model2)

# Get fitted values
fitted_size_lat2 <- fitted(size_lat_model2)

# Plot fitted v residuals: 
plot(fitted_size_lat2, res_size_lat2, pch=21, col="black", bg="palegreen", cex=1.2, cex.main=0.9, xlim=c(-2,1),
     main = "Polynomail Regression (E ~ size*absolute latitude): Residuals vs fitted", xlab="Fitted value", ylab="Residual")

# get model predictions
size_lat_intercept2 <- summary(size_lat_model2)$coefficients[1,1]
size_lat_slope2_1 <- summary(size_lat_model2)$coefficients[2,1]
size_lat_slope2_2 <- summary(size_lat_model2)$coefficients[3,1]

size_lat_predict2 <- data.frame(Size_x_latitude = seq(-400, 200, 10), E = size_lat_intercept2 + size_lat_slope2_1 * seq(-400, 200, 10) + size_lat_slope2_2 * seq(-400, 200, 10)^2)

# Lamba estimates
#################
# Initialise empty vectors for model output
lambda_size_lat2 <- as.vector(matrix(nrow=1000, ncol=1))
logL_size_lat2 <- as.vector(matrix(nrow=1000, ncol=1))
P_val_size_lat2 <- as.vector(matrix(nrow=1000, ncol=1))

tmp <- list()

for (i in 1:1000){
  
  # run model 
  tmp <- phylosig(tree_pruned, res_size_lat2, method="lambda", test=TRUE, nsim=1000, start=NULL, control=list())
  
  # Phylogenetic signal estimates
  lambda_size_lat2[i] <- tmp[[1]]
  
  # Log Likelihood values
  logL_size_lat2[i] <- tmp[[2]]
  
  # P-values
  P_val_size_lat2[i] <- tmp[[4]]
}

# make a df for results
results_size_lat2 <- data.frame(lambda = lambda_size_lat2, logL = logL_size_lat2, P_value = P_val_size_lat2)

# select row with highest Log liklihood value
best_fit_size_lat2 <- results_size_lat2[which.max(results_size_lat2$logL),]
## results say little signal (results are true)

####################
# 5C. Compare models
####################

compare_size_lat <-
  ggplot(data = size_lat_df, aes(x=size*latitude, y =E)) +
  ggtitle("Model comparison: Size x Latitude") +
  ylim(0, 2) +
  xlim(-250, 150) +
  geom_point() +
  geom_line(aes(x= Size_x_latitude, y=E, colour="Linear"), data = size_lat_predict1) +
  geom_line(aes(x= Size_x_latitude, y = E, colour = "Polynomial"), data = size_lat_predict2) +
  scale_colour_manual("", 
                      breaks = c("Linear", "Polynomial"),
                      values = c("red","blue")) 

compare_LS_lin <-
  ggplot(data = size_lat_df, aes(x=size*latitude, y =E)) +
  ggtitle("Linear Regression") +
  labs(x="log(mean body mass) x absolute latitude") +
  ylim(0, 2.5) +
  xlim(-200, 200) +
  geom_point() +
  geom_line(aes(x= Size_x_latitude, y=E), colour="red", data = size_lat_predict1) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 9)) +
  theme(plot.title = element_text(size = 10, face = "bold"))

compare_LS_poly <-
  ggplot(data = size_lat_df, aes(x=size*latitude, y =E)) +
  ggtitle("Polynomial Regression") +
  labs(x="log(mean body mass) x latitude") +
  ylim(0, 2.5) +
  xlim(-200, 250) +
  geom_point() +
  geom_line(aes(x= Size_x_latitude, y = E), colour="blue",data = size_lat_predict2) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 9)) +
  theme(plot.title = element_text(size = 10, face = "bold"))

compare_ls_grid <-plot_grid(compare_LS_lin, compare_LS_poly, labels=c("A", "B"), ncol = 2, nrow = 1)
 

#####
AIC_size_lat_model1 <-AIC(size_lat_model1)# model 1: better AIC
AIC_size_lat_model2 <-AIC(size_lat_model2) 

## BIC
######
BIC_size_lat_model1 <-BIC(size_lat_model1) # model 1: better BIC
BIC_size_lat_model2 <-BIC(size_lat_model2) 

##############
## 6. plots ##
##############

# density: 
density_plot <- ggplot(results, aes(x = lambda, fill = lambda)) +  geom_density(alpha=.2, fill="#FF6666") 

# density plot with histogram:
density_hist_plot <- ggplot(results, aes(x = lambda, fill = lambda)) +  geom_histogram(colour="black", fill="white") +  geom_density(alpha=.2, fill="#FF6666") 

### E
#####
# calculate confidence intervals
sd <- sd(data$E)
mean <- mean(data$E)
n <- length(data$E)
se <- sd/sqrt(n)
t_val <- qt(0.975, df = 54)
x <- se*t_val 
upper <- mean+x
lower <- mean-x

# get points for density line for polygon
E.1 <- 
  ggplot(data, aes(x = E, fill = E)) +  
  geom_density(alpha = 0.6, fill= "#FF6666")  
density.est <- ggplot_build(E.1)$data[[1]][,c(1,2)]

# set points for lower bound
lower.bound <- density.est[density.est$x < lower,]
lower.bound <- rbind(lower.bound, c(0, 0.71258171)) # bind 0 and last value for x
lower.bound <- rbind(lower.bound, c(0, 0.08534317)) # bind 0 and firt value for x

#set points for upper bound
upper.bound <- density.est[density.est$x > upper,]
upper.bound <- rbind(upper.bound, c(0, 2.139951)) # last value for x
upper.bound <- rbind(upper.bound, c(0, 0.9618688))# first value for x

# E distribution
E <- 
  ggplot(data, aes(x = E, fill = E)) +  
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  ggtitle("Distribution of E estimates") +
  geom_vline(aes(xintercept = mean(E), col = "mean"), linetype="dashed", lwd = 0.8) +  
  geom_vline(aes(xintercept = median(E), col = "median"), linetype="dashed", lwd = 0.8) + 
  scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red")) +
  geom_density(alpha = 0.6, fill= "#FF6666") +  
  theme_classic() + 
  geom_polygon(data = lower.bound, 
               aes(x = x, y = z), fill = 'black', alpha =0.3) +
  geom_polygon(data = upper.bound, 
               aes(x = x, y = y), fill = 'black', alpha =0.3)

### lambda density
###################
# calculate confidence intervals
sd_L <- sd(results$lambda)
mean_L <- mean(results$lambda)
n_L <- length(results$lambda)
se_L <- sd_L/sqrt(n_L)
t_val_L <- qt(0.975, df = 999)
x_L <- se_L*t_val_L
upper_L <- mean_L+x_L
lower_L <- mean_L-x_L

# get points for density line for polygon
lam.1 <- 
  ggplot(results, aes(x = lambda, fill = lambda)) +  
  geom_density(alpha = 0.6, fill= "#FF6666")  
density.est_L <- ggplot_build(lam.1)$data[[1]][,c(1,2)]

# set points for lower bound
lower.bound_L <- density.est_L[density.est_L$x < lower_L,]
lower.bound_L <- rbind(lower.bound_L, c(0, 0.4835884309)) # bind 0 and last value for x in lower.bound_L
lower.bound_L<- rbind(lower.bound_L, c(0, 0.0007816912)) # bind 0 and first value for x in lower.bound_L

#set points for upper bound
upper.bound_L <- density.est_L[density.est_L$x > upper_L,]
upper.bound_L <- rbind(upper.bound_L, c(0, 0.9996248)) # last value for x
upper.bound_L <- rbind(upper.bound_L, c(0, 0.5207274))# first value for x

# plot
lam <- 
  ggplot(results, aes(x = lambda, fill = lambda)) +  
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  ggtitle("Distribution of lambda estimates") +
  geom_density(alpha=.5, fill="#FF6666") +  
  geom_vline(aes(xintercept = mean(lambda), col = "mean"), linetype="dashed", lwd = 0.8) +  
  geom_vline(aes(xintercept = median(lambda), col = "median"), linetype="dashed", lwd = 0.8) + 
  scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red")) +
  theme_classic() +
  geom_polygon(data = lower.bound_L, 
               aes(x = x, y = y), fill = 'black', alpha =0.3) + # lower CI
  geom_polygon(data = upper.bound_L, 
               aes(x = x, y = y), fill = 'black', alpha =0.3) # upper CI
  
## resting 
restsp <- subset(data, data$resting_or_active == "resting")

E_rest <- 
  ggplot(restsp, aes(x = E, fill = E)) +  
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  ggtitle("Distribution of E estimates: Resting") +
  geom_vline(aes(xintercept = mean(E), col = "mean"), linetype="dashed", lwd = 0.8) +  
  geom_vline(aes(xintercept = median(E), col = "median"), linetype="dashed", lwd = 0.8) + 
  scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red")) +
  geom_density(alpha=.2, fill="#FF6666") +  
  theme_classic() 

## actice
activesp <- subset(data, data$resting_or_active == "active")

E_active <- 
  ggplot(activesp, aes(x = E, fill = E)) +  
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  ggtitle("Distribution of E estimates: Active") +
  geom_vline(aes(xintercept = mean(E), col = "mean"), linetype="dashed", lwd = 0.8) +  
  geom_vline(aes(xintercept = median(E), col = "median"), linetype="dashed", lwd = 0.8) + 
  scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red")) +
  ylim(0, 2) +
  geom_density(alpha=.2, fill="#FF6666") +  
  theme_classic() 

#larva
larva <- tail(data, n=5)

E_adult<- 
  ggplot(adult, aes(x = E, fill = E)) +  
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  ggtitle("Distribution of E estimates: Adults") +
  geom_vline(aes(xintercept = mean(E), col = "mean"), linetype="dashed", lwd = 0.8) +  
  geom_vline(aes(xintercept = median(E), col = "median"), linetype="dashed", lwd = 0.8) + 
  scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red")) +
  ylim(0, 2) +
  geom_density(alpha=.2, fill="#FF6666") +  
  theme_classic() 

# density plot with histogram with median:
median <- ggplot(results, aes(x = lambda, fill = lambda)) +  geom_density(alpha=.2, fill="#FF6666") +  geom_vline(xintercept = median(x), col = "red", lwd = 0.5)

# box plot
box <- boxplot(results$lambda, 
               main="boxplot",
               ylab="lambda",
               col="orange",
               border="brown")


# bootstrap

bootobject <- boot(data$E, statistic = function_1, R = 1000)

function_1 <- function(data,indices) {
  x <- mean(data)
  return(x)
}

boot.ci(bootobject, conf = 0.95)
bootobject[2]


