## Program for simulations in the rain skedasis project
# Programmed by Chen Zhou

library(kolmim)
## Parameter specification
# Simulation setup
# The number of stations
m <- 5
# The number of observations over time
n <- 2000
# The number of simulations
sim <- 500

# Model setup
# The dependence parameter, can be 0.1 or 0.05
beta <- 0.1
# Number of simulated BM
nBM <- 7
# The volatility of the c function over station
# zero indicates no difference across stations
vos <- 0.4
# The volatility of the c function over time
vot <- 0.5
# Shift in starting point across station
# zero indicats that all station starts and ends at C_j(1)
# Non zero means from station 2, there is a shift 
shift <- 0.5
# EVI of all station
gamma <- 1

## Simulate Brownian motion with parameter beta at location 1:m
simBM <- function(beta,m,nsim){
  eps <- rnorm(m*nsim, mean=0, sd=sqrt(beta))
  dim(eps) <- c(m,nsim)
  results <- apply(eps, 2, cumsum)
  return(results)
}

## Simulate a Brown-Resnick process with fixed number of BM
simBR<- function(beta,m,nBM){
  BMs <- simBM(beta,m,nBM)
  BMs <- BMs-(1:m)%*%t(rep(1,nBM))*beta/2
  expo <- rexp(nBM)
  cexp <- 1/cumsum(expo)
  results <- exp(BMs)%*%diag(cexp)
  results <- apply(results, 1, max)
  return(results)
}

## Pre-simulate Brown-Resnick processes
# To determine the number of Brownian motion needed
# For any beta and m, always run this functon first
# The optimal nBM should be either the maximum of the results
# or a very high quantile of it
presimBR <- function (beta, m, ntest = 1000, start = 5, end = 30) {
  results <- c()
  for(j in 1:ntest){
    BMs <- simBM(beta,m,end)
    BMs <- BMs-(1:m)%*%t(rep(1,end))*beta/2
    expo <- rexp(end)
    cexp <- 1/cumsum(expo)
    temp <- exp(BMs)%*%diag(cexp)
    tempbr <- apply(temp,1,max)
    temppos <- start:end
    for(l in temppos){
      partial <- apply(temp[,1:l],1,max)
      if(sum(partial!=tempbr)==0) break
    }
    results=c(results,l)
  }
  return(results)
}

## C function generating
cfun <- function(m,n,vos,vot,shift){
  Cfun <- 1/m*(1+ vos*sin((1:m)/m*(2*pi)))
  results <- rep(0,m*n)
  dim(results) <- c(n,m)
  for(i in 1:m){
    results[,i]=m*Cfun[i]*(1+vot*sin((1:n)/n*(2*pi)+shift*i))
  }
  return(results)
}

## Data generating process
# Data structure, n rows, m columns
gendata <- function(m,n,vos,vot,shift, gamma, beta, nBM){
  skedasis <- cfun(m,n,vos,vot,shift)
  brprocess <- replicate(n, simBR(beta,m,nBM))
  results <- (skedasis*t(brprocess))^gamma
  return (results)
}

## Statistical procedure
# Estimate the value of C_j(1)
estCfunone <- function(data,k){
  size <- dim(data)
  n <- size[1]
  m <- size[2]
  #alldata <- array(data)
  klen <- length(k)
  results <- array(dim=c(klen,m))
  rankalldata <- sort(data, decreasing = T)
  for(i in 1:klen){
    threshold <- rankalldata[k[i]+1]
    results[i,] <- apply(data>threshold,2,sum)/k[i]
  }
  return(cbind(k,results))
}

testCfunone <- function(data,k){
  size <- dim(data)
  n <- size[1]
  m <- size[2]
  #alldata <- array(data)
  klen <- length(k)
  results <- numeric(klen)
  rankalldata <- sort(data, decreasing = T)
  for(i in 1:klen){
    threshold <- rankalldata[k[i]+1]
    potdata <- data>threshold
    matsigma <- t(potdata)%*%potdata/k[i]
    #B <- diag(1,m,m)-matrix(1,m,m)/m
    B1 <- cbind(diag(1,m-1,m-1),rep(-1,m-1))
    #Sigma <- B1%*%B%*%matsigma%*%t(B)%*%t(B1)
    Sigma <- B1%*%matsigma%*%t(B1)
    Cones <- apply(potdata,2,sum)/k[i]
    compare <- B1%*%Cones
    results[i] <- k[i]*t(compare)%*%solve(Sigma)%*%(compare)
  }
  return(cbind(k,results))
}

testCovertime <- function(data,k){
  size <- dim(data)
  n <- size[1]
  m <- size[2]
  #alldata <- array(data)
  klen <- length(k)
  results <- numeric(klen*m)
  dim(results) <- c(klen,m)
  rankalldata <- sort(data, decreasing = T)
  for(i in 1:klen){
    threshold <- rankalldata[k[i]+1]
    tempc <- apply(data>threshold,2,cumsum)/k[i]
    cones <- tempc[n,]
    tempt <- ((1:n)/n)%*%t(rep(1,m))
    ksdist <- apply(abs(tempc-tempt),2,max)
    results[i,] <- ksdist*sqrt(k[i]*cones)
  }
  return(cbind(k,results))
}
#

## Simulations

# Simulation for point estimates
onesampleCfun <- function(m,n,vos,vot,shift, gamma, beta, nBM, k){
  data <- gendata(m,n,vos,vot,shift, gamma,beta,nBM)
  Cfun <-  estCfunone(data, k)
  return(Cfun[2:(m+1)])
}

multisampleCfun <- function(m,n,vos,vot,shift, gamma, beta, nBM, k, sim){
  results <- replicate(sim, onesampleCfun(m,n,vos,vot,shift, gamma, beta, nBM, k))
  c <- cfun(m,n,vos,vot,shift)
  Cone <- colMeans(c)/m
  boxplot(t(results))
  points(1:m, Cone, col=2, pch = 19, lwd = 3)
}

# Simulation for testing across stations
# Conduct simulation for testing C(1), single sample
onesampletest <- function(m,n,vos,vot,shift, gamma, beta, nBM, k){
  data <- gendata(m,n,vos,vot,shift, gamma,beta,nBM)
  Cfun <-  testCfunone(data, k)
  pvalue <- 1-pchisq(Cfun[,2],df=m-1)
  return(pvalue)
}

# Conduct simulation for testing C(1), multiple sample
multisampletest <- function(m,n,vos,vot,shift, gamma, beta, nBM, k, sim){
  results <- replicate(sim, onesampletest(m,n,vos,vot,shift, gamma, beta, nBM, k))
  return(results)
}

# Conduct simulation for testing constant c over time, single sample
onesampletest2 <- function(m,n,vos,vot,shift, gamma, beta, nBM, k){
  data <- gendata(m,n,vos,vot,shift, gamma,beta,nBM)
  Covertime <-  testCovertime(data, k)
  pvalue <- 1-pkolm(Covertime[,-1],n)
  return(pvalue)
}

# Conduct simulation for testing constant c over time, multiple sample
multisampletest2 <- function(m,n,vos,vot,shift, gamma, beta, nBM, k, sim){
  results <- replicate(sim, onesampletest2(m,n,vos,vot,shift, gamma, beta, nBM, k))
  return(results)
}

# Set estresult to TRUE to obtain the simulation results regarding estimation
# Figure 1 and Figure 3 in the Supplementary Material
estresult <- FALSE

if (estresult){
pdf(file="est1_k200.pdf", width=8, height=5);
multisampleCfun(m,n,vos,vot,shift, gamma, beta, nBM, 200, sim)
dev.off()

pdf(file="est1_k500.pdf", width=8, height=5);
multisampleCfun(m,n,vos,vot,shift, gamma, beta, nBM, 500, sim)
dev.off()

pdf(file="est2_k200.pdf", width=8, height=5);
multisampleCfun(m,n,0.2,vot,shift, gamma, beta, nBM, 200, sim)
dev.off()

pdf(file="est2_k500.pdf", width=8, height=5);
multisampleCfun(m,n,0.2,vot,shift, gamma, beta, nBM, 500, sim)
dev.off()
}

# Set simresult1 to TRUE to obtain the simulation results regarding testing C(1)
# Table 1 and Table 3 in the Supplementary Material
simresult1 <- FALSE
# Set plot to TRUE to obtain the simulation results regarding testing C(1)
# Figure 2 and Figure 4 in the Supplementary Material
plot <- TRUE

if(simresult1){
set.seed(198105)
k <- 200
sim <- 200
results <- numeric(sim*4)
dim(results) <- c(sim,4)

results[,1] <- multisampletest(m,n,0,0,0, gamma, beta, nBM, k, sim)
results[,2] <- multisampletest(m,n,0,vot,shift, gamma, beta, nBM, k, sim)
results[,3] <- multisampletest(m,n,0.2,vot,shift, gamma, beta, nBM, k, sim)
results[,4] <- multisampletest(m,n,0.4,vot,shift, gamma, beta, nBM, k, sim)

allres <- numeric(3*4)
dim(allres) <- c(3,4)
for (i in 1:4){
  allres[1,i] <- sum(results[,i]<0.01)/sim
  allres[2,i] <- sum(results[,i]<0.05)/sim
  allres[3,i] <- sum(results[,i]<0.1)/sim
}

if(plot){
pdf(file="QQ1_k200.pdf", width=5, height=5)
selectlow <- results[,1] < 0.1
lowvalues <- results[selectlow,1]
qqplot(qunif(ppoints(sum(selectlow)), max=0.1), lowvalues, main ="",xlab="Theoretical quantiles", ylab="Empirical quantiles", xlim=c(0,0.1), ylim=c(0,0.1))
abline(a=0,b=1,lty=2)
box()
dev.off()

pdf(file="QQ2_k200.pdf", width=5, height=5)
selectlow <- results[,2] < 0.1
lowvalues <- results[selectlow,2]
qqplot(qunif(ppoints(sum(selectlow)), max=0.1), lowvalues, main ="",xlab="Theoretical quantiles", ylab="Empirical quantiles", xlim=c(0,0.1), ylim=c(0,0.1))
abline(a=0,b=1,lty=2)
box()
dev.off()
}

k <- 500
sim <- 200
results <- numeric(sim*4)
dim(results) <- c(sim,4)

results[,1] <- multisampletest(m,n,0,0,0, gamma, beta, nBM, k, sim)
results[,2] <- multisampletest(m,n,0,vot,shift, gamma, beta, nBM, k, sim)
results[,3] <- multisampletest(m,n,0.2,vot,shift, gamma, beta, nBM, k, sim)
results[,4] <- multisampletest(m,n,0.4,vot,shift, gamma, beta, nBM, k, sim)

allres1 <- numeric(3*4)
dim(allres1) <- c(3,4)
for (i in 1:4){
  allres1[1,i] <- sum(results[,i]<0.01)/sim
  allres1[2,i] <- sum(results[,i]<0.05)/sim
  allres1[3,i] <- sum(results[,i]<0.1)/sim
}

if(plot){
  pdf(file="QQ1_k500.pdf", width=5, height=5)
  selectlow <- results[,1] < 0.1
  lowvalues <- results[selectlow,1]
  qqplot(qunif(ppoints(sum(selectlow)), max=0.1), lowvalues, main ="",xlab="Theoretical quantiles", ylab="Empirical quantiles", xlim=c(0,0.1), ylim=c(0,0.1))
  abline(a=0,b=1,lty=2)
  box()
  dev.off()
  
  pdf(file="QQ2_k500.pdf", width=5, height=5)
  selectlow <- results[,2] < 0.1
  lowvalues <- results[selectlow,2]
  qqplot(qunif(ppoints(sum(selectlow)), max=0.1), lowvalues, main ="",xlab="Theoretical quantiles", ylab="Empirical quantiles", xlim=c(0,0.1), ylim=c(0,0.1))
  abline(a=0,b=1,lty=2)
  box()
  dev.off()
}

save(allres,allres1,file="test.RData")
}