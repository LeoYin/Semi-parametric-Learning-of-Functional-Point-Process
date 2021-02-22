.libPaths(c('/home/liang/R/x86_64-redhat-linux-gnu-library/3.5'))
command_line_arg = as.numeric(commandArgs(trailingOnly = TRUE)[1])
procid = command_line_arg+1;
setwd('/home/liang/Desktop/MFPCA/Simulation')
source("MFPCA_final.R")
lbd <- 0
ubd <- 1


Eigen.x <- list()
funx.1 <- function(t){
  1+t*0
}
funx.2 <- function(t){
  sqrt(3)*(1-2*t)
}

Eigen.x[[1]] <- funx.1
Eigen.x[[2]] <- funx.2


Eigen.y <- list()
funy.1 <- function(t){
  1+t*0
}
funy.2 <- function(t){
  sqrt(2)*sin(2*pi*t)
}

Eigen.y[[1]] <- funy.1
Eigen.y[[2]] <- funy.2

Eigen.z <- list()
funz.1 <- function(t){
  1+t*0
}
funz.2 <- function(t){
  sqrt(2)*sin(4*pi*t)
}

Eigen.z[[1]] <- funz.1
Eigen.z[[2]] <- funz.2

eta.x <- c(0.5,0.2)
eta.y <- c(0.5,0.2)
eta.z <- c(0.5,0.2)*0.1

lambda0 <- function(t){
  0.3*cos(2*pi*t)+1
}
h_seq <- seq(0.02,0.1,by=0.02)

filename = paste("./run/iteration=",procid,".csv",sep='')


for(n in c(50,100,150,200,250,300)){
  for(m in c(50,100,150,200,250,300)){
  
set.seed(procid)
data_twolevel <- Sim_twolevel(n,m,lbd,ubd,lambda0,Eigen.x,Eigen.y,Eigen.z,eta.x,eta.y,eta.z)

Process <- data_twolevel$Process
xi.x <- data_twolevel$xi.x
xi.y <- data_twolevel$xi.y
ngrid <- 100
obj.cv <- MFPCA_foldcv(Process,lbd,ubd,h_seq,ngrid,kern="epanechnikov",nlarge=1000,npc.x=2,npc.y=2,fold=10,cv.seed=procid,nbreaks=31)

time <- obj.cv$time
xi.xhat <- obj.cv$xi.x
xi.yhat <- obj.cv$xi.y

eta.xhat <- obj.cv$eta.x
eta.yhat <- obj.cv$eta.y
eta.zhat <- obj.cv$eta.z

Eigen.xhat <- obj.cv$Eigen.x
Eigen.yhat <- obj.cv$Eigen.y
Eigen.zhat <- obj.cv$Eigen.z

h.x <- obj.cv$h.x
h.y <- obj.cv$h.y
h.z <- obj.cv$h.z


covhat <- cor <- lad <- sign.x <- sign.y <- numeric()
grids <- seq(0,1,l=ngrid)

fhat <- Eigen.xhat[,1]
ftrue <- funx.1(grids)
pos <- ifelse(sum(abs(fhat-ftrue))<=sum(abs(fhat+ftrue)),1,-1)
sign.x[1] <- pos

lad[1] <- mean(abs(Eigen.xhat[,1]*pos-funx.1(grids)))

fhat <- Eigen.xhat[,2]
ftrue <- funx.2(grids)
pos <- ifelse(sum(abs(fhat-ftrue))<=sum(abs(fhat+ftrue)),1,-1)
sign.x[2] <- pos

lad[2] <- mean(abs(Eigen.xhat[,2]*pos-funx.2(grids)))

fhat <- Eigen.yhat[,1]
ftrue <- funy.1(grids)
pos <- ifelse(sum(abs(fhat-ftrue))<=sum(abs(fhat+ftrue)),1,-1)
sign.y[1] <- pos

lad[3] <- mean(abs(Eigen.yhat[,1]*pos-funy.1(grids)))


fhat <- Eigen.yhat[,2]
ftrue <- funy.2(grids)
pos <- ifelse(sum(abs(fhat-ftrue))<=sum(abs(fhat+ftrue)),1,-1)
sign.y[2] <- pos

lad[4] <- mean(abs(Eigen.yhat[,2]*pos-funy.2(grids)))


fhat <- Eigen.zhat[,1]
ftrue <- funz.1(grids)
pos <- ifelse(sum(abs(fhat-ftrue))<=sum(abs(fhat+ftrue)),1,-1)

lad[5] <- mean(abs(Eigen.zhat[,1]*pos-funz.1(grids)))


fhat <- Eigen.zhat[,2]
ftrue <- funz.2(grids)
pos <- ifelse(sum(abs(fhat-ftrue))<=sum(abs(fhat+ftrue)),1,-1)

lad[6] <- mean(abs(Eigen.zhat[,2]*pos-funz.2(grids)))

cor[1] <- cor(xi.x[,1],xi.xhat[,1]*sign.x[1])
cor[2] <- cor(xi.x[,2],xi.xhat[,2]*sign.x[2])
cor[3] <- cor(xi.y[,1],xi.yhat[,1]*sign.y[1])
cor[4] <- cor(xi.y[,2],xi.yhat[,2]*sign.y[2])

covhat[1] <- cov(xi.x[,1],xi.xhat[,1]*sign.x[1])
covhat[2] <- cov(xi.x[,2],xi.xhat[,2]*sign.x[2])
covhat[3] <- cov(xi.y[,1],xi.yhat[,1]*sign.y[1])
covhat[4] <- cov(xi.y[,2],xi.yhat[,2]*sign.y[2])


Result <- c(lad,cor,eta.xhat[1:2],eta.yhat[1:2],eta.zhat[1:2],time,h.x,h.y,h.z,n,m,procid,covhat)
write.table(t(Result),file=filename,sep=",",row.names=F,col.names=F,append = TRUE)
}
}

