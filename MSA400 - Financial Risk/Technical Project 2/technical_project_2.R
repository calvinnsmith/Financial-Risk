#############################################
################## TASK 1 ###################
#############################################

OBLIGORS <- 35
LOAN <- 3
DEFAULT_LOSSES <- 0.6
DEFAULT_PROBABILITY <- 0.04



a <- 0.2
b <- 4.59

VaR_beta <- LOAN*DEFAULT_LOSSES*OBLIGORS*qbeta(alpha,a,b)

#############################################
################## TASK 3 ###################
#############################################
library("extraDistr")


DEFAULTS <- 0:35

# low correlation (corr = 0.1)
a1 <- 0.36
b1 <- 8.64

mod1 <- c(a1,b1,a1/(a1+b1),0.1)

# high correlation (corr = 0.8)
a2 <- 0.01
b2 <- 0.24

mod2 <- c(a2,b2,a2/(a2+b2),0.8)

summary <- rbind(mod1,mod2)
rownames(summary) <- c("Model 1","Model 2")
colnames(summary) <- c("a","b","p","rho")

# Fitting density functions
betamodel_1 <- dbbinom(DEFAULTS,OBLIGORS,alpha = a1,beta = b1)

betamodel_2 <- dbbinom(DEFAULTS,OBLIGORS,alpha = a2,beta = b2)

binomial <- dbinom(DEFAULTS,OBLIGORS,0.04)

plot(DEFAULTS,betamodel_1,"o",col = "green",bg = "green",lwd = 1,ylim = c(0,0.2),main = "Beta binomial" ,xlab = "Number of defaults",ylab = "Probability")
lines(DEFAULTS,betamodel_2,"o",col = "red",bg = "red",lwd = 1)
lines(DEFAULTS,binomial,"o",col = "blue",bg="blue",lwd = 1)

color.names <- c("green","red","blue")
legend("top",c("Mixed binomial 35, Beta(0.36,8.64)","Mixed binomial 35, Beta(0.01,0.24)","Binomial (35,0.04)"),fill = color.names)

# Value at Risk
quantiles <- c(0.95,0.99,0.999)

VaR_beta1 <- LOAN*DEFAULT_LOSSES*OBLIGORS*qbeta(quantiles,a1,b1)

VaR_beta2 <- LOAN*DEFAULT_LOSSES*OBLIGORS*qbeta(quantiles,a2,b2)

# Expected Shortfall
ES_model1 <- numeric(3)
ES_model2<- numeric(3)

for (i in 1:length(quantiles)){
  ES_model1[i] <-  ((LOAN*OBLIGORS*DEFAULT_LOSSES)/(1-quantiles[i]))*integrate(qbeta,lower = quantiles[i],upper = 1,shape1 = a1,shape2=b1)$value
  ES_model2[i] <-  ((LOAN*OBLIGORS*DEFAULT_LOSSES)/(1-quantiles[i]))*integrate(qbeta,lower = quantiles[i],upper = 1,shape1 = a2,shape2=b2)$value
  
}

# Making tables for LaTex
model1_var_es <- rbind(VaR_beta1,ES_model1)
rownames(model1_var_es) <- c("VaR","ES")
colnames(model1_var_es) <- c("95%","99%","99.9")

model2_var_es <- rbind(VaR_beta2,ES_model2)
rownames(model2_var_es) <- c("VaR","ES")
colnames(model2_var_es) <- c("95%","99%","99.9")

library(xtable)

xtable(model1_var_es)
xtable(model2_var_es)


##################################################################
################## TASK 3: MONTE CARLO SIMULATION ################
##################################################################

# function for Monte Carlo simulation
monte_carlo <- function(nsim,a,b,obligors,loan,default_losses){
  L <- numeric(nsim)
  for (j in 1:nsim){
    Z <- rbeta(1,shape1 = a,shape2 = b)
    U <- runif(obligors,0,1)
    X <- numeric(obligors)
    for (i in 1:obligors){
      
      if (U[i] <= Z){
        X[i] <- 1
      }
      else{
        X[i] <- 0
      }
    }
    L[j] <- sum(X)*loan*default_losses
  }
  return(L)
}

n <- 1000000

# Beta Binomial model 1, VaR
L_model1 <- monte_carlo(n,a1,b1,OBLIGORS,LOAN,DEFAULT_LOSSES)
MC_VaR_model1 <- quantile(L_model1,quantiles)

# Beta Binomial model 2, VaR
L_model2 <- monte_carlo(n,a2,b2,OBLIGORS,LOAN,DEFAULT_LOSSES)
MC_VaR_model2 <- quantile(L_model2,quantiles)


# function for calculating Expected Shortfall

ES <- function(vaR,losses){
  
  expectedShortfall <- numeric(length(vaR))
  
  for (i in 1:length(vaR)){
    expectedShortfall[i] <- sum(losses[which(losses >= vaR[i])])/(length(which(losses >= vaR[i])))
  }
  return(expectedShortfall)
}


# Beta Binomial Model 1, ES
MC_ES_model1 <- ES(MC_VaR_model1,L_model1)

# Beta Binomial Model 1, ES
MC_ES_model2 <- ES(MC_VaR_model2,L_model2)

# Creating tables for LaTex
summaryMC_1 <- rbind(MC_VaR_model1,MC_ES_model1)
rownames(summaryMC_1) <- c("VaR","ES")

summaryMC_2 <- rbind(MC_VaR_model2,MC_ES_model2)
rownames(summaryMC_1) <- c("VaR","ES")

# plots of MC distributions

hist(L_model1,freq = FALSE, breaks = 50,main = "MC loss distribution - Beta Model 1", xlab = "Loss")

hist(L_model2,freq = FALSE, breaks = 50, main = "MC loss distribution - Beta Model 2", xlab = "Loss")

hist(sort(L_model1)[900000:1000000],freq = FALSE,breaks = 50,main ="MC loss distribution(10% highest) - Beta Model 1", xlab = "Loss")

hist(sort(L_model2)[900000:1000000],freq = FALSE, breaks = 50, main = "MC loss distribution(10% highest - Beta  Model 2", xlab = "Loss")
