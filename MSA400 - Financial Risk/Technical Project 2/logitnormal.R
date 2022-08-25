
library(logitnorm)
library(tidyverse)
library(xtable)


#### TASK 2 ####

#function to calculate correlation
correlation <- function(expectation,variance) {
  variance/(expectation*(1-expectation))
}

mu <- seq(-20,10,length.out = 100)
sigma <- seq(0,10,length.out = 100)

#matrix for storing mu, sigma , E(p(Z)) and Var(p(Z))
param <- matrix(0,nrow = length(mu)*length(sigma),ncol = 4)
colnames(param) <- c("mu","sigma","mean","var")
count = 0;

#grid search over mu and sigma to find the
# parameters of our models
for(i in 1:length(mu)){
  for(j in 1:length(sigma)){
    count = count + 1;
    param[count,3] <- momentsLogitnorm(mu[i],sigma[j])[1]
    param[count,4] <- momentsLogitnorm(mu[i],sigma[j])[2]
    param[count,1] <- mu[i]
    param[count,2] <- sigma[j]
  }
}

tol = 0.001

param <- as.data.frame(param)

#obtaining the combinations of mu and sigma that gives a mean of approximately 0.04
#and a variance depending on which correlation is desired
param_new <- param %>% filter(abs(mean-0.04) <= 0.001)

## correlation = 0.8
high_set <- param_new[2,]
correlation(high_set$mean,high_set$var)


## correlation = 0.1
low_set <- param_new[46,]
correlation(low_set$mean,low_set$var)

##### TASK 3 #####

DEFAULTS <- 0:35
obligors <- 35
loan <- 3
DEFAULT_LOSSES <- 0.6
loss <- obligors*DEFAULT_LOSSES*loan


x <- DEFAULTS/obligors

#inverse of the mixing distribution function
p_inverse <- function(mu,sigma,x){
  (1/sigma)*(log(x/(1-x))-mu)
}

# plots of loss distributions

# los distribution (LPA) N(p^-1(x)
loss_dist1 <- dnorm(p_inverse(high_set$mu,high_set$sigma,x),0,1)
loss_dist2 <- dnorm(p_inverse(low_set$mu,low_set$sigma,x),0,1)

plot(DEFAULTS,loss_dist1,"l",ylim = c(0,0.4))
lines(DEFAULTS,loss_dist2,"l")

plot(DEFAULTS,loss_dist1,"o",col = "green",bg = "green",lwd = 1,ylim = c(0,0.4),main = "Logit normal" ,xlab = "Number of defaults",ylab = "Probability")
lines(DEFAULTS,loss_dist2,"o",col = "red",bg = "red",lwd = 1)
lines(DEFAULTS,binomial,"o",col = "blue",bg="blue",lwd = 1)

color.names <- c("green","red","blue")
legend("top",c("Logit normal 1","Logit normal 2","Binomial"),fill = color.names)

## value at risk and expected shortfall 

quantiles <- c(0.95,0.99,0.999)

# function for calculating Value At Risk
ValueAr <- function(quantiles,mu,sigma,loss)
{
  valueAtRisk <- (exp(qnorm(quantiles)*sigma + mu)/(1+exp(qnorm(quantiles)*sigma + mu)))
}

# value at risk
VaR_model1 <- loss*ValueAr(quantiles,high_set$mu,high_set$sigma,loss)
VaR_model2 <- loss*ValueAr(quantiles,low_set$mu,low_set$sigma,loss)

#expected shortfall
ES_model1 <- numeric(3)
ES_model2 <- numeric(3)

  for(i in 1:length(quantiles)){
    ES_model1[i] <- (loss/(1-quantiles[i]))*integrate(ValueAr,lower = quantiles[i],upper = 1,high_set$mu,high_set$sigma,loss)$value
    ES_model2[i] <- (loss/(1-quantiles[i]))*integrate(ValueAr,lower = quantiles[i],upper = 1,low_set$mu,low_set$sigma,loss)$value
  }





