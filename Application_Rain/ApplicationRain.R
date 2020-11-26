## Program for application in the rain skedasis project
# (c) Chen Zhou, Claudia Neves, November 2020
#Choose season (Summer or Winter) and toggle option on/off

#===== Creates maps of the EVI estimates========
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(devtools)
library(stringr)

library(grid)
library(gridExtra)
require(eva)
require(stats)
require(reshape)

rm(list=ls())

#========== !! Input required in here !! =====================
title= "Winter"
datafile= "data/dataIW64.csv" # Toggle between "data/dataIS64.csv" for Summer and "data/dataIW64.csv" for Winter; There are 64 possible stations altogether
timefile="timeIW"
threshold= 80     #  Upper bound for altitude=elevation, application in the paper aimed at 80
k=1000    # Overall threshold
kmax= k*1.5
kmin= 300
j1=16 #j1 and j2 are stations intervening in pairwise dependence, if TRUE in the below
j2=21

compareC= TRUE
EVIunique= TRUE # MLE estimation of overall gamma; output in EVIresults.csv
EVIuniqueCI= !TRUE # to add CIs to the previous; output in EVIresults.csv
PairDependence= !TRUE
EVIresult= !TRUE
PValues= TRUE # FALSE returns observed values of test statistic
Maps= !TRUE
CIPlot= !TRUE
TestExtremal= TRUE
TEstConstEVIinStation= !TRUE
TEstConstEVIBetweenStation = !TRUE
#========== Initialise =====================

allstations<- read.csv("data/All64Stns.csv", header=TRUE)
data <- read.csv(paste(datafile), header = FALSE)

m= as.integer(length(data)) 
n= as.integer(length(data[[1]]))

testTimeTrend= TRUE
#==========Select stations with altitude below threshold===============

select <- allstations$elevation < threshold
stationno <- (1:m)[select]
m=sum(select)
print(paste("no. of stations entering analysis is m=", m))
data <- as.matrix(data[, select]) # filters out data from selected stations
A<- data.frame(allstations$no, allstations$lat, allstations$lon, allstations$elevation)
colnames(A)<- c("no","latitude","longitude","elevation")
Stns<- A[select==TRUE,]
rm(A)


#======== Randomize the data

set.seed(19810527)
error <- runif(n*m,min=-0.05,max=0.05)
dim(error) <- c(n,m)
data <- data+error 
rm(error)
print(paste("number of days per stations", n))
#======== Order all data
alldata <- array(data) #data is a matrix 3571xm if threshold= 80, yielding 49 stations
rankalldata <- sort(data, decreasing = T)

#======== Obtain Cj(1), plot later ========

Cj1 = matrix(0,m,0)
threshold = rankalldata[k+1]
print(paste("overall thershold is",threshold,"mm"))
Cj1= cbind(Cj1,apply(data>threshold,2,sum) )
maxC=max(Cj1)
Stns$C1=Cj1 
print(paste("max frequency is",maxC/k))

###########################################################
######### Begin of Statistical functions ##################

#======== Obtain alphaj1j2(j/k,l/k) ==== makes use of Cj1 from the above =======
alphaj1j2 <- function(j, l, k, evi){
  c1 = matrix(1,n,m)
  c2 = matrix(1,n,m)
  c3 = matrix(1,n,m)
  c1[data<rankalldata[j]]<- 0 # Select j exceedances above X_{N-j:N} ---
  c2[data<rankalldata[l]]<- 0 # Select l exceedances above X_{N-l:N} ---
  c3[data<rankalldata[k]]<- 0 # Select k exceedances above the baseline overall threshold X_{N-k:N} ---
  rjl=matrix(0,nrow = m, ncol = m)  #this is an mXm matrix ---
  for (j1 in 1:m){
    for (j2 in 1:m){
      rjl[j1,j2]= (j*l/k^2)^(-evi-1 )*sum(as.vector(c1[,j1]*c2[,j2]))/k + sum(as.vector(c3[,j1]*c3[,j2]))/k - (j/k)^(-evi-1) * sum(as.vector(c1[,j1]*c3[,j2]))/k - (l/k)^(-evi-1)* sum(as.vector(c3[,j1]*c2[,j2]))/k
    }
  }  
  for (j1 in 1:m){
    rjl[j1,j1]= ((j*l/k^2)^(-evi-1 ) * min(j,l)/k +1 - (j/k)^(-evi) - (l/k)^(-evi) )* Cj1[j1]/k
  }
  rm(c1,c2,c3)
  aux= sum(rjl) # dependence across space is assumed present
  #  aux= sum(rjl[!outer(1:m, 1:m, "-")]) # full independence setting
  rm(rjl)
  return(aux) 
}

#--- station-wise homogeneity test
testCfun <- function(k, thres){
  
  s <- seq(1 / n, 1, 1 / n)
  results <- array(dim=c(m,2))
  
  for(j in 1:m){
    tempdata <- data[,j]
    sdata <- sort(tempdata, decreasing = T)
    Cest <- cumsum((tempdata>thres))
    results[j,1] = sqrt(Stns$C1[j])* max(abs(Cest/Stns$C1[j]-s))
    results[j,2] <- Stns$C1[j]/(n) * (t(Cest/Stns$C1[j]-s)%*% (Cest/Stns$C1[j]-s))
  }
  return(results)
}

#--- overall station test for a given k, elapsed time ----
testCfunone <- function(k){
  size <- dim(data)
  n <- size[1]
  m <- size[2]
  klen <- length(k)
  results <- numeric(klen)
  
  for(i in 1:klen){
    threshold <- rankalldata[k[i]+1]
    potdata <- data>threshold
    matsigma <- t(potdata)%*%potdata/k[i]
    B1 <- cbind(diag(1,m-1,m-1),rep(-1,m-1))
    Sigma <- B1%*%matsigma%*%t(B1)
    Cones <- apply(potdata,2,sum)/k[i]
    compare <- B1%*%Cones
    results[i] <- k[i]*t(compare)%*%solve(Sigma)%*%(compare)
  }
  return(cbind(k,results))
}

#--- overall EVI estimation for a given k, elapse time ----
EVIpot <- function(k0){
  size <- dim(data)
  n = size[1]
  m = size[2]
  klen <- length(k0)
  threshold =rankalldata[k0+1]
  
  rankalldata[rankalldata<threshold]<- NA 
  rankalldata= rankalldata[is.finite(rankalldata)]
  gp<- gpdFit(rankalldata, nextremes=k0, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)
  
  return(as.double(gp$par.ests[2]))
}

#--- asymptotic confidence bands for the overall EVI estimator for a given k, elapse time ---

EVIpotStdErr<- function(k0, evi0){
  aux1= 0
  for (j in 1:(k0-1)){
    print(paste("current j is", j))
    for(l in 1:(k0-1)){
      inner= alphaj1j2(j,l,k0,evi0)
      aux1= aux1 + ((j/k0)^evi0 - (2*evi0 + 1)* (j/k0)^(2*evi0)) * inner * ((l/k0)^evi0 - (2*evi0 + 1)* (l/k0)^(2*evi0)) /k0^2
    }
  }
  evistdv=  (evi0+1)^2/evi0* sqrt( aux1/k0)
  return(evistdv)
}

#--- Function used for the sliding blocks estimator
### Input
## b: integer, block size
## Y: array of the observed process
## n: integer, number of observations
### Output
## theta, the estimate
slidingEstimator <- function(b,n,Y){
  Fhat = ecdf(Y)    # define empirical CDF
  
  # compute the sliding blocks estimator
  slSum = 0;
  for(i in 1:(n+1-b)){
    slSum = slSum + b*(1-Fhat(max(Y[(i):(i+b-1)])))
  }
  theta = (n+1-b)/(slSum)
  return(theta)
}

####### end of Statistical Functions ######
#############################################################

#-------- Obtain rj1j2(1,1) for selected pair of stations j1 and j2 ==== makes use of Cj1 from the above =======
if(PairDependence){
  
  ck = matrix(1,n,m)
  rj1j2<- numeric(kmax)
  kseq <- seq(from=1, to = kmax, by=1)
  kmin=300
  kmax=1500
  for(k0 in kmin:kmax){
    ck = matrix(1,n,m)
    ck[data<rankalldata[k0+1]]<- 0 
    rj1j2[k0]= sum(as.vector(ck[,j1]*ck[,j2]))/k0 * m 
    rm(ck)
  }
  dt0<- data.frame(kseq, rj1j2)
  colnames(dt0) <- c("K", "Rk")
  p0 <- ggplot(dt0, aes(K, Rk)) + geom_line(color="#FC4E07", size=2) + labs(x = "k", y = " ")+ scale_y_continuous(limits=c(0,1.0), breaks=seq(0, 1.0, 0.2)) +  scale_x_continuous(limits=c(kmin,kmax), breaks=seq(0,kmax, 200)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))  + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20))
  grid.arrange(p0, nrow = 1, top=textGrob(paste("Pairwise dependence: j1= Hude, j2= Schiffdorf --", title), gp=gpar(fontsize=20,font=3) ), left= textGrob("m x sigma j1,j2(k)", rot= 90, gp=gpar(fontsize=20) ) )
  
}

#-------- Test for serial dependence =======
if(TestExtremal){
  maxdata<-apply(data,1,FUN=max)
  seqb<-seq(from=10, to=100, by=1)
  lb<-length(seqb)
  thetaest<-numeric(lb)
  thetasd<-thetaest
  for (i in 1:lb){
    k<-n/seqb[i]
    thetaest[i]<-slidingEstimator(seqb[i],n,maxdata)
    thetasd[i]<-sqrt(0.2726)/sqrt(k)
  }
  
  if (timefile == 'timeIS'){
    plotfile <- "IS_extremal.pdf"
  } else{
    plotfile<- "IW_extremal.pdf"
  }
  
  pdf(plotfile, width = 8, height = 5)
  
  matplot(seqb, cbind(thetaest, thetaest-1.96*thetasd, thetaest+1.96*thetasd), type="l", lty=c(1,2,2), col=c(1,2,2), xlab="block size", ylab=expression(theta))
  abline(h=1,lty=2)
  
  dev.off()
}

#--- Testing C(1) across stations for different values of k & plot
if(compareC){
  print("Testing for a trend in space...")
  lseq <- seq(from=kmin, to = kmax, by=1)
  Cfun <- testCfunone(lseq) 
  pvalue <- 1-pchisq(Cfun[,2],df=m-1)
  dt1<- data.frame(lseq, pvalue)
  colnames(dt1) <- c("K", "pvalue")
  rm(lseq)
  p1 <- ggplot(dt1, aes(K, pvalue)) + geom_line(color="#FC4E07", size=2) + labs(x = "k", y = " ")+ scale_y_continuous(limits=c(0,0.07), breaks=seq(0, 1.0, 0.01))+ geom_hline(yintercept = 0.05, color="grey") +  scale_x_continuous(limits=c(kmin,kmax), breaks=seq(0,kmax, 200)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))  + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20)) 
  grid.arrange(p1, nrow = 1, top=textGrob(paste("Detecting a trend in space --", title), gp=gpar(fontsize=20,font=3) ), left= textGrob("Observed p-values for T_n", rot= 90, gp=gpar(fontsize=20) ) )
  
}

#--- test for constant gamma within station  & plot ----
if(TEstConstEVIinStation){
  b=7
  siglevel= 0.025/49
  inb= floor(n/b)
  kimin= 51
  kimax= 60
  
  c2=matrix(0,1,m)
  c0<- as.matrix(data[1:inb,])
  potdata<- apply(c0, 2, sort, decreasing=T)
  cb1=matrix(0,0,m)
  cb2=matrix(0,0,m)
  print("Testing for constant gamma within station...")
  for (kaux in kimin:kimax){
    for (j in 1:m){
      c2[j]=  as.double(gpdFit(potdata[,j], nextremes=kaux, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2])
    }
    cb1= rbind(cb1,c2)
    cb2= rbind(cb2,c2*c2)
  }
  for (l in 2:b)
  {
    print(paste("block=", l))
    c0<- as.matrix(data[((l-1)*inb+1):(l*inb),])
    potdata<- apply(c0, 2, sort, decreasing=T)
    c11=matrix(0,0,m)
    c12=matrix(0,0,m)
    for (kaux in kimin:kimax){
      for (j in 1:m){
        c2[j]=  as.double(gpdFit(potdata[,j], nextremes=kaux, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2])
      }
      c11= rbind(c11,c2)
      c12= rbind(c12,c2*c2)
    }
    cb1= cb1 + c11
    cb2= cb2 + c12
  }
  rm(c0,c2,c11,c12)
  sumb<- cb2-cb1*cb1/b
  kseq<- seq(from= kimin, to = kimax, by=1)
  sumb= sweep(sumb, 1, kseq, `*`)

  Station<- seq(from=1, to= m, by= 1)
  c0= matrix(0,0,m)
  c0= rbind(Station)
  c0= rbind(c0,sumb/((1+cb1/b)*(1+cb1/b)))
  dttt<- data.frame(t(c0))
  colnames(dttt)  <- c("Station", kseq)
  dttt.m= melt(dttt, id.vars ="Station", measure.vars = c("51", "52", "53", "54", "55", "56", "57", "58", "59", "60"))
  names(dttt.m)<- c("Station", "K_tilde", "Test_Statistic")
  p3<- ggplot(dttt.m, aes(x = Station, y = Test_Statistic, color = K_tilde)) + geom_point(size=3) + geom_hline(yintercept = qchisq(1-siglevel, df=b-1), color="grey") + theme(axis.text.y = element_text(size = 20))  + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20))
  grid.arrange(p3, nrow = 1, top=textGrob(paste(title), gp=gpar(fontsize=20,font=3) ) )
  rm(dttt)
}

#--- test for constant gamma BETWEEN stations & plot  ----
if(TEstConstEVIBetweenStation){
  siglevel= 0.025/49
  kimin= 30
  kimax= 252
  
  c2=matrix(0,1,m)
  clb=matrix(0,0,m)
  cub=matrix(0,0,m)
  c0<- as.matrix(data)
  potdata<- apply(c0, 2, sort, decreasing=T)
  rm(c0)
  print("Testing for constant gamma between stations...")
  for (kaux in kimin:kimax){
    for (j in 1:m){
      c2[j]=  as.double(gpdFit(potdata[,j], nextremes=kaux, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2])
    }
    clb= rbind(clb, c2 - qnorm(1-siglevel)/sqrt(kaux)* (1 + c2))
    cub= rbind(cub, c2 + qnorm(1-siglevel)/sqrt(kaux)* (1 + c2))
  }
  clb= as.vector(clb)
  cub= as.vector(cub)
  kseq= as.matrix(rep(seq(from= kimin, to = kimax, by=1), m), m*(kimax-kimin +1), 1)
  clb= cbind(kseq, clb)
  cub= cbind(kseq, cub)
  points= rbind(clb, cub)
  bound <- c(rep("LB", m*(kimax-kimin +1) ), rep("UB", m*(kimax-kimin +1) ) )
  dt1<- data.frame(bound, points)
  names(dt1)<- c("B","k","value")
  p4<- ggplot(dt1, aes(x = k, y = value, color = B)) + geom_point(size=2)+ theme(legend.position = "none") #+ scale_shape_manual(values = c(2, 6)) + scale_shape_manual(values = c(1, 2, 6))
  p4<- p4 + scale_y_continuous(limits=c(-0.6,0.8), breaks=seq(-1.0, 1.0, 0.2))+ scale_x_continuous(limits=c(80,250), breaks=seq(0, 2000, 20)) + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20))+ theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y=element_blank()) +theme(legend.title=element_blank(), legend.text=element_text(size = 12))
  grid.arrange(p4, nrow = 1, top=textGrob(paste(title), gp=gpar(fontsize=20,font=3) ) )
}

# --- Estimating common EVI for different values of k ranging from kmin to kmax & plot
if(EVIunique){
  Est<- numeric(kmax)
  lv<- numeric(kmax)
  UB <-numeric(kmax)
  LB <-numeric(kmax)
  UBInd <- numeric(kmax)
  LBInd <- numeric(kmax)
  for (i in kmin:kmax){
    Est[i] <- EVIpot(i)
    lv[i]<- i
    if(EVIuniqueCI){
      se= EVIpotStdErr(i, Est[i])
      UB[i] = Est[i] + qnorm(0.975)* se
      LB[i] = Est[i] - qnorm(0.975)* se
      UBInd[i] <- Est[i] + qnorm(0.975)* (1+Est[i])/sqrt(lv[i]) # Lower bound for the CI in the independence setting
      LBInd[i] <- Est[i] - qnorm(0.975)* (1+Est[i])/sqrt(lv[i]) # Upper bound for the CI in the independence setting
    }
  }
  dt2<- data.frame(lv,Est, LB, UB, LBInd, UBInd)
  print(paste("overall EVI estimate:", Est[k], LB[k], UB[k], LBInd[k], UBInd[k]))
  colnames(dt2) <- c("K", "EVI", "LB", "UB", "LBInd", "UBInd")
  write.csv(dt2,"EVIresults.csv")
  rm(lv)
  p2 <- ggplot(dt2, aes(K, EVI)) + geom_line(color="#C4961A", size=2) + labs(x = "k", y = " ")+ scale_y_continuous(limits=c(-0.04,0.1), breaks=seq(-1.0, 1.0, 0.01))+ geom_hline(yintercept = Est[k], color="grey")+ geom_vline(xintercept = 1000, color="grey")+  scale_x_continuous(limits=c(kmin,kmax), breaks=seq(0,kmax, 200))+ theme(axis.text.y = element_text(size = 20))  + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20))
  grid.arrange(p2, nrow = 1, top=textGrob(paste("Extreme-value index --", title), gp=gpar(fontsize=20,font=3) ), left= textGrob("EVI estimates", rot= 90, gp=gpar(fontsize=20) ) )
  if(EVIuniqueCI){
    p3 <- ggplot(dt2, aes(K, LB)) + geom_line(color="#C4961A", size=2) + labs(x = "k", y = " ")+ scale_y_continuous(limits=c(-0.04,0.1), breaks=seq(-1.0, 1.0, 0.01))+ geom_hline(yintercept = Est[k], color="grey")+ geom_vline(xintercept = 1000, color="grey")+  scale_x_continuous(limits=c(kmin,kmax), breaks=seq(0,kmax, 200))+ theme(axis.text.y = element_text(size = 20))  + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20))
    grid.arrange(p3, nrow = 1, top=textGrob(paste("Extreme-value index LB--", title), gp=gpar(fontsize=20,font=3) ), left= textGrob("EVI estimates", rot= 90, gp=gpar(fontsize=20) ) )
  }  
}

# --- Test constant scedasis at each station; simulates critical values from Brownian bridge
if(testTimeTrend){
  
  set.seed(19810527)
  grid <- 5000
  simbm <- 2000
  increments <- rnorm(grid*simbm)
  dim(increments) <- c(grid,simbm)
  simbms <- apply(increments, 2, cumsum)/sqrt(grid)
  
  gstep <- (1:grid)/grid
  
  simbms <- simbms - gstep %*% t(simbms[grid,])
  
  t1c <- apply(abs(simbms),2,max)
  t2c <- apply(simbms^2, 2, sum)/grid  # for the A-D goodness-of-fit test
  
  testing <- testCfun(k, threshold)
  pvalue <- testing
  
  if(PValues){  
    for(j in 1:m){
      pvalue[j,1] <- sum(t1c>testing[j,1])/simbm
      pvalue[j,2] <- sum(t2c>testing[j,2])/simbm
    }
  }
}
Stns$KStest= pvalue[,1] 
Stns$ADtest= pvalue[,2]

# --- Find estimated EVI for all stations ==========================
if(EVIresult){
  source("eviCZ.r")
  allevi <- numeric(m*2)
  dim(allevi) <- c(m,2)
  for (i in 1:m){
    temp <- data[,i]
    results <- mleest(temp, k)
    allevi[i,1] <- results[2]
    allevi[i,2] <- results[3]
  }
  
  Stns$EVI= allevi[,1] 
  Stns$StdEVI= allevi[,2]
}
# Write dataframe on to a file ====== Output =========

if (title =="Summer"){    
  evifile <- paste0("IS_k",k,".csv")
  print(evifile)
} else{
  evifile <- paste0("IW_k",k,".csv")
  print(evifile)
}
write.csv(Stns, file=evifile)
rm(Stns)


# ==== Plots ====== Output ===

#-------------------- EVI with confidence bands ----------------------
# Can be compiled on its own, all input is read form a GXXX.csv file

if(CIPlot){
# "#FC4E07" in geom_line for summer ; "#00AFBB" ibidem for winter
  stnEVI <- read.csv("G49W.csv", header = TRUE)
  if(title=="Summer"){
    stnEVI <- read.csv("G49S.csv", header = TRUE)
    print(paste("Title is", title))
  }
  LBInd<- stnEVI$EVI - qnorm(0.975)* (1+stnEVI$EVI)/ sqrt(stnEVI$K) # independence
  UBInd<- stnEVI$EVI + qnorm(0.975)* (1+stnEVI$EVI)/ sqrt(stnEVI$K) # independence
  LBound<- stnEVI$EVI - qnorm(0.975)* (1+stnEVI$EVI)^2/stnEVI$EVI* sqrt(stnEVI$AUX1/stnEVI$K) # dependence
  UBound<- stnEVI$EVI + qnorm(0.975)* (1+stnEVI$EVI)^2/stnEVI$EVI* sqrt(stnEVI$AUX1/stnEVI$K) # dependence
  p21 <- ggplot(stnEVI, aes(K, EVI)) + labs(x = "k", y = " ")+ scale_y_continuous(limits=c(-0.15,0.25), breaks=seq(-1.5, 1.5, 0.05))+  scale_x_continuous(limits=c(kmin+250,kmax-100), breaks=seq(0,kmax-1, 100))+ theme(axis.text.y = element_text(size = 20))  + theme(axis.text.x = element_text(size = 20)) + theme(axis.title.x = element_text(size = 20))
  p22 <- p21  + geom_ribbon(aes(ymin = LBound, ymax = UBound), fill = "#C4961A", alpha= 0.8) + geom_ribbon(aes(ymin = LBInd, ymax = UBInd), fill = "grey", alpha= 0.5) + geom_vline(xintercept = 1000, color="grey") +  geom_line(aes(y = EVI), color="#00AFBB", size=2) 
  grid.arrange(p22, nrow = 1, top=textGrob(paste(title), gp=gpar(fontsize=20,font=3) ), left= textGrob("EVI estimates", rot= 90, gp=gpar(fontsize=20) ) )
}

# ====  Maps  ===
if(Maps){
  register_google(key = "XXX")
  sq_map <- get_googlemap(center = c(lon=9.0, lat=53.0), zoom=7, maptype = "terrain", api_key = 'XXX')
  b<- c(0.05, 0.15, 0.25, 0.5, 0.75)
  mapt<- ggmap(sq_map) + geom_point(data = stn , mapping = aes(x = longitude, y = latitude, color = KStest), size=4)+ geom_text(data = stn, aes(x = longitude, y = latitude, label = elevation, angle = 0, hjust = 1.5))
  mapt1<- mapt + scale_color_gradient2(low = "blue", mid="red", high = "green", breaks= b)
  grid.arrange(mapt1, nrow = 1, top=paste("Testing --", title))
}
