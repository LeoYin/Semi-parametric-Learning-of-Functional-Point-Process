library(R.matlab)
library(spatstat)
library(compiler)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)

enableJIT(3)
kern.intensity <- function(process,bwd,r,kern="epanechnikov"){
  n.r <- length(r)
  edge <- pkernel(r/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((r-1)/bwd,kernel = kern,sd=sqrt(1/5))
  rho <- numeric()
  for(i in 1:n.r){
    ind <- abs(process-r[i]) < bwd  
    u <-   abs(process[ind]-r[i])/bwd
    rho[i] <- sum(dkernel(u,kernel = kern,sd=sqrt(1/5))/bwd/edge[i])
  }
  rho  
}


EigenScore <- function(process,Eigen,sig2,grids,Idx=NULL){
  process_all <- sort(unlist(process))  
  n <- length(process)
  p <- dim(Eigen)[2]
  freq_all <- table(cut(process_all,breaks=grids))
  ngrids <- length(grids)
  Eigen_mid <- (Eigen[-1,,drop=FALSE]+Eigen[-ngrids,,drop=FALSE])/2
  si2_mid <- (sig2[-1]+sig2[-ngrids])/2
  score <- NULL
  prob <- NULL
  loglik <- NULL
  if(is.null(Idx)) Idx <- 1:n
  
  for(i in Idx){
    process_sel <- process[[i]]
    freq_sel <- table(cut(process_sel,breaks=grids))
    
    fn <- function(par){
      tmp <-  exp(Eigen_mid%*%par-log(n-1)-si2_mid/2)
      prob_mid <-c(tmp/(tmp+1))
      val <- sum(log(prob_mid)*freq_sel)+sum(log(1-prob_mid)*(freq_all-freq_sel))
      -val
    }
    par0 <- rep(0,p)
    obj <- optim(par0,fn,method="BFGS")
    score <- rbind(score,obj$par)
    loglik <- c(loglik,obj$value)
    tmp <-  exp(Eigen%*%obj$par-log(n-1)-sig2/2)
    prob <- rbind(prob,c(tmp/(tmp+1)))
    if(i%%500==0) print(i)
  }
  list(score=score,loglik=loglik,prob=prob)
}



###Simulate 1-d inhomogenous poisson process
Sim_Poisson <- function(lbd,ubd,intensity,intensity.max=NULL){
  ###Get the rhomax if not specified
  if(is.null(intensity.max)){
    grids <- seq(lbd,ubd,l=1000)  
    intensity.max=max(intensity(grids))  
  }
  intensity.max <- max(intensity.max,15)
  ###Simulate total number of points
  count <- rpois(1,intensity.max*(ubd-lbd)) 
  if(count==0) process <- numeric(0) else{
    process_nothining <- sort(runif(count,lbd,ubd))  
    prob <- intensity(process_nothining)/intensity.max
    process <- process_nothining[prob>runif(count)]
  }
  
 process  
}

#########Simulate two-level data
Sim_twolevel <- function(n,m,lbd,ubd,lambda0,Eigen.x,Eigen.y,Eigen.z,eta.x,eta.y,eta.z,ngrid=100,rho=0.95){
logIntensity <-   Process <- matrix(rep(list(),m*n),n,m)
grids <- seq(lbd,ubd,l=ngrid)
  p.x <- length(eta.x)
  p.y <- length(eta.y)
  p.z <- length(eta.z)
  xi.x <- matrix(rnorm(p.x*n,0,1),n,p.x)
  xi.y <- matrix(rnorm(p.y*m,0,1),m,p.y)
  for(j in 1:(m-1)){
    xi.y[j+1,1] <- xi.y[j,1]*rho+rnorm(1)*sqrt(1-rho^2)  
  }
  
  xi.x <- xi.x%*%diag(sqrt(eta.x))
  xi.y <- xi.y%*%diag(sqrt(eta.y))
  
  
  for(i in 1:n){
    for(j in 1:m){
      xi.z <- rnorm(p.z,0,1)*sqrt(eta.z)
      
      intensity <- function(t){
        tmp.x <-tmp.y<-tmp.z<- 0
        for(k in 1:p.x)
          tmp.x <- tmp.x+Eigen.x[[k]](t)*xi.x[i,k]
        for(k in 1:p.y)
          tmp.y <- tmp.y+Eigen.y[[k]](t)*xi.y[j,k]
        for(k in 1:p.z)
          tmp.z <- tmp.z+Eigen.z[[k]](t)*xi.z[k]
        lambda0(t)*exp(tmp.x+tmp.y+tmp.z)
      }
      
      logIntensity[i,j] <- list(log(intensity(grids)))
      
      sample <- Sim_Poisson(lbd,ubd,intensity,intensity.max=NULL)
      if(length(sample)!=0)   Process[i,j] <- list(sample) 
    }
  } 
  list(Process=Process,xi.x=xi.x,xi.y=xi.y,logIntensity=logIntensity)
}


#########Simulate two-level data
Sim_twolevel_bivariate <- function(n,m,lbd,ubd,lambda0,Eigen.x,Eigen.y,Eigen.z,Eta.x,Eta.y,Eta.z,ngrid=100){
  logIntensity1 <-   Process1 <- matrix(rep(list(),m*n),n,m)
  logIntensity2 <-   Process2 <- matrix(rep(list(),m*n),n,m)
  grids <- seq(lbd,ubd,l=ngrid)
  p.x <- dim(Eta.x)[1]/2
  p.y <- dim(Eta.y)[1]/2
  p.z <- dim(Eta.z)[1]/2
  xi.x <- matrix(rnorm(p.x*n*2,0,1),n,p.x*2)
  xi.y <- matrix(rnorm(p.y*m*2,0,1),m,p.y*2)
  err.y <- matrix(rnorm(p.y*m*2,0,1),m,p.y*2)
 
#     for(j in 1:(m-1)){
#      xi.y[j+1,1] <- xi.y[j,1]*0.5+err.y[j,1]*sqrt(0.75)  
#    }
  
  R.x <- chol(Eta.x)
  R.y <- chol(Eta.y)
  R.z <- chol(Eta.z)
  xi.x <- xi.x%*%R.x
  xi.y <- xi.y%*%R.y


  
  xi.x1 <- xi.x[,1:p.x,drop=FALSE]
  xi.x2 <- xi.x[,-(1:p.x),drop=FALSE]
  
  xi.y1 <- xi.y[,1:p.y,drop=FALSE]
  xi.y2 <- xi.y[,-(1:p.y),drop=FALSE]
  
  for(i in 1:n){
    for(j in 1:m){
      xi.z <- rnorm(p.z*2,0,1)%*%R.z
      xi.z1 <- xi.z[1:p.z]
      xi.z2 <- xi.z[-(1:p.z)]
      
      intensity1 <- function(t){
        tmp.x <-tmp.y<-tmp.z<- 0
        for(k in 1:p.x)
          tmp.x <- tmp.x+Eigen.x[[k]](t)*xi.x1[i,k]
        for(k in 1:p.y)
          tmp.y <- tmp.y+Eigen.y[[k]](t)*xi.y1[j,k]
        for(k in 1:p.z)
          tmp.z <- tmp.z+Eigen.z[[k]](t)*xi.z1[k]
        lambda0(t)*exp(tmp.x+tmp.y+tmp.z)
      }
      
      intensity2 <- function(t){
        tmp.x <-tmp.y<-tmp.z<- 0
        for(k in 1:p.x)
          tmp.x <- tmp.x+Eigen.x[[k]](t)*xi.x2[i,k]
        for(k in 1:p.y)
          tmp.y <- tmp.y+Eigen.y[[k]](t)*xi.y2[j,k]
        for(k in 1:p.z)
          tmp.z <- tmp.z+Eigen.z[[k]](t)*xi.z2[k]
        lambda0(t)*exp(tmp.x+tmp.y+tmp.z)
      }
      
      logIntensity1[i,j] <- list(log(intensity1(grids)))
      logIntensity2[i,j] <- list(log(intensity2(grids)))
      
      sample1 <- Sim_Poisson(lbd,ubd,intensity1,intensity.max=NULL)
      sample2 <- Sim_Poisson(lbd,ubd,intensity2,intensity.max=NULL)
      if(length(sample1)!=0)   Process1[i,j] <- list(sample1) 
      if(length(sample2)!=0)   Process2[i,j] <- list(sample2) 
    }
  } 
  list(Process1=Process1,Process2=Process2,xi.x1=xi.x1,xi.y1=xi.y1,xi.x2=xi.x2,xi.y2=xi.y2,logIntensity1=logIntensity1,logIntensity2=logIntensity2)
}

###Estimate the covariance function of the MFPCA (use common bandwidth)
MFPCA_Cov <- function(Process,lbd,ubd,bwd,ngrid,kern="epanechnikov",nlarge=10^(10)){
  
  n <- dim(Process)[1]  
  m <- dim(Process)[2]  
  grids <- seq(lbd,ubd,l=ngrid)
  
  A2 <- A <- B <- C <- D <- matrix(0,ngrid,ngrid)
  edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
  edge <- outer(edge,edge,FUN="*")
  
  Kh <- function(t){
    dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
  }
  
  for(i in 1:n){
    for(j in 1:m){
      process <- Process[i,j][[1]] 
      if(!is.null(process)){
        cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
        tmp1 <- apply(tmp,2,sum)
        tmp2 <- outer(tmp1,tmp1,FUN="*")
        A <- A+(tmp2-Matrix::t(tmp)%*%tmp)
        A2 <- A2 + tmp2
      }
    }
    if(i%%100==0) print(i)
  }
  
  
  for(i in 1:n){
    process.tmp <- unlist(Process[i,])
    if(length(process.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process <- grid.mid
      freq <- c(table(cut(process.tmp,breaks=grid.large)))
    }else{
      process <- process.tmp
      freq <- rep(1,length(process.tmp))
    }
    cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")
    B <- B + tmp2
  }
  
  B <- B-A2
  
  for(j in 1:m){
    process.tmp <- unlist(Process[,j])
    if(length(process.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process <- grid.mid
      freq <- c(table(cut(process.tmp,breaks=grid.large)))
    }else{
      process <- process.tmp
      freq <- rep(1,length(process.tmp))
    }
    cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")
    C <- C + tmp2
  }
  
  C <- C-A2
  
  
  process.tmp <- unlist(Process)
  if(length(process.tmp)>nlarge){
    grid.large <- seq(lbd,ubd,l=nlarge+1)
    grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
    process <- grid.mid
    freq <- c(table(cut(process.tmp,breaks=grid.large)))
  }else{
    process <- process.tmp
    freq <- rep(1,length(process.tmp))
  }
  cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
  tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
  tmp1 <- apply(tmp,2,sum)
  tmp2 <- outer(tmp1,tmp1,FUN="*")
  
  D <- tmp2 - B - C - A2
  
  A <- A/(n*m*edge)
  B <- B/(n*m*(m-1)*edge)
  C <- C/(n*m*(n-1)*edge)
  D <- D/(n*m*(n-1)*(m-1)*edge)
  
  ####Find the covariance matrix
  Cov.x <- log(B)-log(D)
  Cov.y <- log(C)-log(D)
  Cov.z <- log(A)+log(D)-log(B)-log(C)
  
  ###Find Eigen functions and eigen values
  obj.x <- eigen(Cov.x)
  Eigen.x <- obj.x$vectors/sqrt(grids[2]-grids[1])
  obj.y <- eigen(Cov.y)
  Eigen.y <- obj.y$vectors/sqrt(grids[2]-grids[1])
  obj.z <- eigen(Cov.z)
  Eigen.z <- obj.z$vectors/sqrt(grids[2]-grids[1])
  
  list(A=A,B=B,C=C,D=D,A2=A2,Cov.x=Cov.x,Cov.y=Cov.y,Cov.z=Cov.z,Eigen.x=Eigen.x,Eigen.y=Eigen.y,Eigen.z=Eigen.z)
}


###Estimate the covariance function of the MFPCA (use common bandwidth)
MFPCA_Cov_diffh <- function(lbd,ubd,ngrid,A,B,C,D){
  grids <- seq(lbd,ubd,l=ngrid)
  
  ####Find the covariance matrix
  Cov.x <- log(B)-log(D)
  Cov.y <- log(C)-log(D)
  Cov.z <- log(A)+log(D)-log(B)-log(C)
  
  ###Find Eigen functions and eigen values
  obj.x <- eigen(Cov.x)
  Eigen.x <- obj.x$vectors/sqrt(grids[2]-grids[1])
  obj.y <- eigen(Cov.y)
  Eigen.y <- obj.y$vectors/sqrt(grids[2]-grids[1])
  obj.z <- eigen(Cov.z)
  Eigen.z <- obj.z$vectors/sqrt(grids[2]-grids[1])
  
  list(A=A,B=B,C=C,D=D,Cov.x=Cov.x,Cov.y=Cov.y,Cov.z=Cov.z,Eigen.x=Eigen.x,Eigen.y=Eigen.y,Eigen.z=Eigen.z)
}


###Choose the bandwidth using leave-a-block-out cross validation
MFPCA_clik_cv.matrix_parallel <- function(Process,lbd,ubd,bwd,ngrid,kern="epanechnikov",nlarge=10^(10),npc.x=2,npc.y=2,numCores=6){
  n <- dim(Process)[1]  
  m <- dim(Process)[2]  
  grids <- seq(lbd,ubd,l=ngrid)
  n.h <- length(h_seq)
  cv.rho<- cv.A <- cv.B <- cv.C<- cv.D <- numeric()
  count_all <- table(cut(unlist(Process),breaks=grids))
  grids_mid <- (grids[-1]+grids[-ngrid])/2
  
  A2 <- A <- B <- C <- D <- matrix(0,ngrid,ngrid)
  
  edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
  edge <- outer(edge,edge,FUN="*")
  
  Kh <- function(t){
    dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
  }
  
  # numCores <- detectCores()
  
  
  cl <- makeCluster(numCores)
  #  registerDoParallel(cl)    
  registerDoSNOW(cl)
  ###Create progress report
  iterations <- n
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  result <- foreach(i = 1:n,.options.snow = opts)%dopar%{
    library(spatstat)
    source("MFPCA_final.R")
    tmp.A <- tmp.A2 <- 0
    for(j in 1:m){
      process <- Process[i,j][[1]] 
      if(!is.null(process)){
        cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
        tmp1 <- apply(tmp,2,sum)
        tmp2 <- outer(tmp1,tmp1,FUN="*")
        tmp.A <- tmp.A+(tmp2-Matrix::t(tmp)%*%tmp)
        tmp.A2 <- tmp.A2 + tmp2
      }
    }
    list(tmp.A=tmp.A,tmp.A2=tmp.A2)
  }
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  ##tag1 
  print("Step 1 finished")
  
  A2.out <- list()
  A <- A2 <- 0
  
  for(i in 1:n){
    A2.out[[i]] <- result[[i]]$tmp.A2
    A <- A + result[[i]]$tmp.A
    A2 <- A2 + A2.out[[i]]
  }
  
  rm(result)
  
  cl <- makeCluster(numCores)
  #  registerDoParallel(cl)    
  registerDoSNOW(cl)
  ###Create progress report
  iterations <- m
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  result <-  foreach(j = 1:m,.options.snow = opts)%dopar%{
    library(spatstat)
    source("MFPCA_final.R")
    tmp.A2.y <- 0
    for(i in 1:n){
      process <- Process[i,j][[1]] 
      if(!is.null(process)){
        cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
        tmp1 <- apply(tmp,2,sum)
        tmp2 <- outer(tmp1,tmp1,FUN="*")
        tmp.A2.y <- tmp.A2.y + tmp2
      }
    }
    tmp.A2.y
  }
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  print("Step 2 finished")
  
  A2.out.y <- list()
  for(j in 1:m){
    A2.out.y[[j]] <- result[[j]]
  }
  rm(result)
  
  
  
  cl <- makeCluster(numCores)
  #  registerDoParallel(cl)    
  registerDoSNOW(cl)
  ###Create progress report
  iterations <- n
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  result <-  foreach(i = 1:n,.options.snow = opts)%dopar%{
    library(spatstat)
    source("MFPCA_final.R")
    process.tmp <- unlist(Process[i,])
    if(length(process.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process <- grid.mid
      freq <- c(table(cut(process.tmp,breaks=grid.large)))
    }else{
      process <- process.tmp
      freq <- rep(1,length(process.tmp))
    }
    cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")
    list(tmp.B=tmp2,tmp.rho=c(tmp1)/sqrt(diag(edge)),process=process.tmp)
  }
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  print("Step 3 finished")
  
  rho.out <- B.out <- list()
  B <- rho <- 0
  for(i in 1:n){
    B.out[[i]] <- (result[[i]]$tmp.B-A2.out[[i]])/(m*(m-1)*edge)
    B <- B + B.out[[i]]
    rho.out[[i]] <- result[[i]]$tmp.rho/m
    rho <- rho + rho.out[[i]]
  }
  
  rm(result)
  
  
  cl <- makeCluster(numCores)
  #  registerDoParallel(cl)    
  registerDoSNOW(cl)
  ###Create progress report
  iterations <- m
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  result <-  foreach(j = 1:m,.options.snow = opts)%dopar%{
    library(spatstat)
    source("MFPCA_final.R")
    process.tmp <- unlist(Process[,j])
    if(length(process.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process <- grid.mid
      freq <- c(table(cut(process.tmp,breaks=grid.large)))
    }else{
      process <- process.tmp
      freq <- rep(1,length(process.tmp))
    }
    cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")
    list(tmp.C=tmp2,tmp.rho.y=c(tmp1)/sqrt(diag(edge)),process=process.tmp)
  }
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  print("Step 4 finished")
  
  rho.out.y <- C.out <- list()
  C <- 0
  for(j in 1:m){
    C.out[[j]] <- (result[[j]]$tmp.C-A2.out.y[[j]])/(n*(n-1)*edge)
    C <- C + C.out[[j]]
    rho.out.y[[j]] <- result[[j]]$tmp.rho/n
  }
  
  rm(result)
  
  A <- as.matrix(A/((n*m-1)*edge))
  D <- outer(rho/(n),rho/(n),FUN="*")
  
  B <- as.matrix(B/n)
  C <- as.matrix(C/m)
  
  rho <- rho/n
  
    Cov.x <- as.matrix(log(B)-log(D))

    Cov.y <- as.matrix(log(C)-log(D))

   Cov.z <- log(A)+log(D)-log(B)-log(C)


###Find Eigen functions and eigen values
obj.x <- eigen(Cov.x)
Eigen.x <- obj.x$vectors/sqrt(grids[2]-grids[1])
eta.x <- obj.x$values*(grids[2]-grids[1])
obj.y <- eigen(Cov.y)
Eigen.y <- obj.y$vectors/sqrt(grids[2]-grids[1])
eta.y <- obj.y$values*(grids[2]-grids[1])

obj.z <- eigen(Cov.z)
Eigen.z <- obj.z$vectors/sqrt(grids[2]-grids[1])
eta.z <- obj.z$values*(grids[2]-grids[1])

sig2.x <- diag(Cov.x)
sig2.y <- diag(Cov.y)

  list(A=A,B=B,C=C,D=D,rho=rho,B.out=B.out,rho.out=rho.out,rho.out.y=rho.out.y,C.out=C.out,Cov.x=Cov.x,Cov.y=Cov.y,Cov.z=Cov.z,Eigen.x=Eigen.x,Eigen.y=Eigen.y,Eigen.z=Eigen.z,eta.x=eta.x,eta.y=eta.y,eta.z=eta.z)
}


###Choose the bandwidth using leave-a-fold-out cross validation
FoldCV <- function(process.x,B.out,rho.out,lbd,ubd,npc,fold=10,seed=10,nbreaks=21){
  
  n <- length(B.out)
  ngrid <- dim(B.out[[1]])[1]
  grids <- seq(lbd,ubd,l=ngrid)
  grids_mid <- (grids[-1]+grids[-ngrid])/2
  
  B_all <- rho_all <- 0
  for(i in 1:n){
    B_all <- B_all + B.out[[i]]
    rho_all <- rho_all + rho.out[[i]]
  }
  
  set.seed(seed)
  if(fold==1) folds <- rep(1,n) else
  folds <- cut(sample(n),breaks=fold,labels=FALSE)
  cv <- 0
  r <- quantile(unlist(process.x),probs = seq(0,1,l=nbreaks))
  
  for(k in 1:fold){
    idx <- which(folds==k,arr.ind=TRUE) 
    B_k <- B_all
    rho_k <- rho_all
    n.out <- length(idx)
    
    for(i in 1:n.out){
      B_k <- B_k - B.out[[idx[i]]] 
      rho_k <- rho_k - rho.out[[idx[i]]] 
    }
    
    B_k <- B_k/(n-n.out)
    rho_k <- c(rho_k/(n-n.out))
    D_k <- outer(rho_k,rho_k,FUN="*")
    
    Cov.x <- as.matrix(log(B_k)-log(D_k))
    Eigen.x <- eigen(Cov.x)$vectors/sqrt(grids[2]-grids[1])
    sig2 <- diag(Cov.x)
    Eigen <- Eigen.x[,1:npc]
    cv <- cv + GofTest(lbd,ubd,process.x[idx],Eigen,sig2,rho = rho_k,r)/fold
 #   print(k)
  }
  list(cv2=mean(cv),cv=sum(cv*diff(r,1)),cv_seq=cv,r=r) 
}


Aggre_process <- function(Process,dim=1,numCores=1){
  n <- dim(Process)[dim]
  cl <- makeCluster(numCores)
  #  registerDoParallel(cl)    
  registerDoSNOW(cl)
  ###Create progress report
  iterations <- n
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  if(dim==1){
    process <- foreach(i=1:n,.options.snow = opts) %dopar% {
      unlist(Process[i,])
    }
  } else if(dim==2){ 
    process <- foreach(i=1:n,.options.snow = opts) %dopar% {
      unlist(Process[,i])
    }
  }
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  process
}

GofTest <- function(lbd,ubd,Process,Eigen,sig2,rho,r,reference=FALSE){
  n <- length(Process)
  n.pc <- dim(Eigen)[2]
  ngrid <- dim(Eigen)[1]
  grids <- seq(lbd,ubd,l=ngrid)
  
  process_all <- unlist(Process)
  freq_all <- table(cut(process_all,breaks=grids))
  Eigen_mid <- (Eigen[-1,,drop=FALSE]+Eigen[-ngrid,,drop=FALSE])/2
  si2_mid <- (sig2[-1]+sig2[-ngrid])/2
  n.r <- length(r)
  lambda_est <- matrix(NA,n,ngrid)
  Count_emp <- matrix(NA,n,n.r-1)
  
  for(i in 1:n){
    process_sel <- Process[[i]]
    Count_emp[i,] <- table(cut(process_sel,breaks = r))
  }
  
  counts_all <- colSums(Count_emp)
  Prob.emp <- scale(Count_emp,center=FALSE,scale=counts_all)

  if(reference){
    nsim <- 200
    KL_ref <- numeric(n.r-1)
    for(j in 1:(n.r-1)){
      prob_tmp <-  rmultinom(nsim,counts_all[j],prob =Prob.emp[,j] )/counts_all[j]
      KL_ref[j] <- mean(colSums(prob_tmp*(log((prob_tmp+10^(-100))/(Prob.emp[,j]+10^(-100))))))
    }
    return(KL_ref)
  }else{
    for(i in 1:n){
      process_sel <- Process[[i]]
      freq_sel <- table(cut(process_sel,breaks=grids))
      
      fn <- function(par){
        tmp <-  exp(Eigen_mid%*%par-log(n-1)-si2_mid/2)
        prob_mid <-c(tmp/(tmp+1))
        val <- sum(log(prob_mid)*freq_sel)+sum(log(1-prob_mid)*(freq_all-freq_sel))
        -val
      }
      par0 <- rep(0,n.pc)
      obj <- optim(par0,fn,method="BFGS")
      lambda_est[i,] <- exp(Eigen%*%obj$par-sig2/2)*rho
    }
    lambda_all <- colSums(lambda_est)
    
    Prob_est <-  matrix(NA,n,n.r-1)
    
    for(i in 1:n){
      process <- Process[[i]]
      for(j in 1:(n.r-1)){
        ind <- (grids>=r[j])&(grids<=r[j+1])
        Prob_est[i,j] <- sum(lambda_est[i,ind])/sum(lambda_all[ind])
      }
      }
     Prob_est <- scale(Prob_est,center=FALSE,scale=colSums(Prob_est))
    val <- -colSums(Prob.emp*(log(Prob_est)-log(Prob.emp+10^(-100))))
    return(val)
  }
}


MFPCA_foldcv <- function(Process,lbd,ubd,h_seq,ngrid,kern="epanechnikov",nlarge=10^(10),npc.x=2,npc.y=2,fold=10,cv.seed,nbreaks=21){
  n <- dim(Process)[1]  
  m <- dim(Process)[2]  
  grids <- seq(lbd,ubd,l=ngrid)
  n.h <- length(h_seq)
  cv.x <- cv.y <- numeric()
  count_all <- table(cut(unlist(Process),breaks=grids))
  grids_mid <- (grids[-1]+grids[-ngrid])/2
  
  time <- numeric()
  for(l in 1:n.h){
    bwd <- h_seq[l]
    A2 <- A <- B <- C <- D <- matrix(0,ngrid,ngrid)
    
    edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
    edge <- outer(edge,edge,FUN="*")
    
    Kh <- function(t){
      dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
    }
    
    start_time <- Sys.time()
    
    A2.out <- list()
    A2.out.y <- vector(mode = "list", length = m)
    for(i in 1:n){
      tmp.A2 <- matrix(0,ngrid,ngrid)
      for(j in 1:m){
        process <- Process[i,j][[1]] 
        if(!is.null(process)){
          cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
          tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
          tmp1 <- apply(tmp,2,sum)
          tmp2 <- outer(tmp1,tmp1,FUN="*")
          A <- A+(tmp2-Matrix::t(tmp)%*%tmp)
          A2 <- A2 + tmp2
          tmp.A2 <- tmp.A2 + tmp2
          if(is.null(A2.out.y[[j]])) A2.out.y[[j]] <- tmp2 else
            A2.out.y[[j]] <- A2.out.y[[j]] +tmp2
        }
      }
      A2.out[[i]] <- tmp.A2
      if(i%%50==0) print(i)
    }
    
    process.x <- rho.out <- B.out <- list()
    rho <- 0
    for(i in 1:n){
      process.x[[i]] <- process.tmp <- unlist(Process[i,])
      if(length(process.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process <- grid.mid
        freq <- c(table(cut(process.tmp,breaks=grid.large)))
      }else{
        process <- process.tmp
        freq <- rep(1,length(process.tmp))
      }
      cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
      tmp1 <- apply(tmp,2,sum)
      tmp2 <- outer(tmp1,tmp1,FUN="*")
      B <- B + tmp2
      B.out[[i]] <- (tmp2-A2.out[[i]])/(m*(m-1)*edge)
      rho <- rho + c(tmp1)/sqrt(diag(edge))
      rho.out[[i]] <- c(tmp1)/sqrt(diag(edge))/m
    }
    
    B <- B-A2
    
    process.y <- rho.out.y <- C.out <- list()
    
    for(j in 1:m){
      process.y[[j]] <- process.tmp <- unlist(Process[,j])
      if(length(process.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process <- grid.mid
        freq <- c(table(cut(process.tmp,breaks=grid.large)))
      }else{
        process <- process.tmp
        freq <- rep(1,length(process.tmp))
      }
      cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
      tmp1 <- apply(tmp,2,sum)
      tmp2 <- outer(tmp1,tmp1,FUN="*")
      C <- C + tmp2
      C.out[[j]] <- (tmp2-A2.out.y[[j]])/(n*(n-1)*edge)
      rho.out.y[[j]] <- c(tmp1)/sqrt(diag(edge))/n
    }
    
    C <- C-A2
    
    rho <- rho/(m*n)
    A <- as.matrix(A/(n*m*edge))
    D <- outer(rho,rho,FUN="*")
    B <- as.matrix(B/((n*m*(m-1))*edge))
    C <- as.matrix(C/((n*m*(n-1))*edge))
    end_time <- Sys.time()
    
    time[l] <- difftime(end_time,start_time,units="secs")
    
    fold.x <- min(fold,floor(n/10))
    cv.x[l] <-  FoldCV(process.x,B.out,rho.out,lbd,ubd,npc.x,fold=fold.x,seed=cv.seed,nbreaks)$cv
    fold.y <- min(fold,floor(m/10))
    cv.y[l] <-  FoldCV(process.y,C.out,rho.out.y,lbd,ubd,npc.y,fold=fold.y,seed=cv.seed,nbreaks)$cv
    cv.z <- cv.x+cv.y
    
    if(which.min(cv.x)==l) 
    {
      Cov.x <- as.matrix(log(B)-log(D))
      h.x <- bwd
      #     Cov.z.x <- log(A)+log(D)-log(B)-log(C)
    }
    if(which.min(cv.y)==l)
    {
      Cov.y <- as.matrix(log(C)-log(D))
      h.y <- bwd
      #      Cov.z.y <- log(A)+log(D)-log(B)-log(C)
    }
    
    if(which.min(cv.z)==l) 
    {
      Cov.z <- log(A)+log(D)-log(B)-log(C)
      h.z <- bwd
      rho_opt <- rho
    }
    
    
  }
  
  #    if(which.min(cv.x)>which.min(cv.y)) Cov.z <- Cov.z.y else Cov.z <- Cov.z.x
  #    Cov.z <- (Cov.z.x+Cov.z.y)/2
  
  ###Find Eigen functions and eigen values
  obj.x <- eigen(Cov.x)
  Eigen.x <- obj.x$vectors/sqrt(grids[2]-grids[1])
  eta.x <- obj.x$values*(grids[2]-grids[1])
  obj.y <- eigen(Cov.y)
  Eigen.y <- obj.y$vectors/sqrt(grids[2]-grids[1])
  eta.y <- obj.y$values*(grids[2]-grids[1])
  
  obj.z <- eigen(Cov.z)
  Eigen.z <- obj.z$vectors/sqrt(grids[2]-grids[1])
  eta.z <- obj.z$values*(grids[2]-grids[1])
  
  sig2.x <- diag(Cov.x)
  sig2.y <- diag(Cov.y)
  
  xi.x <- EigenScore(process.x,Eigen.x[,1:npc.x,drop=FALSE],sig2.x,grids,Idx=NULL)$score
  xi.y <- EigenScore(process.y,Eigen.y[,1:npc.y,drop=FALSE],sig2.y,grids,Idx=NULL)$score
  
  list(cv.x=cv.x,cv.y=cv.y,cv.z=cv.z,h.x=h.x,h.y=h.y,h.z=h.z,rho=rho_opt,Cov.x=Cov.x,Cov.y=Cov.y,Cov.z=Cov.z,Eigen.x=Eigen.x,Eigen.y=Eigen.y,Eigen.z=Eigen.z,eta.x=eta.x,eta.y=eta.y,eta.z=eta.z,xi.x=xi.x,xi.y=xi.y,time=mean(time))
}

###Estimate the cross covariance function
  MFPCA_CrossCov <- function(Process1,Process2,lbd,ubd,bwd,ngrid,kern="epanechnikov",nlarge=10^(10)){
    n <- dim(Process1)[1]  
    m <- dim(Process1)[2]  
    grids <- seq(lbd,ubd,l=ngrid)
    grids_mid <- (grids[-1]+grids[-ngrid])/2

      A <- B <- C <- D <- matrix(0,ngrid,ngrid)
      
      edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
      edge <- outer(edge,edge,FUN="*")
      
      Kh <- function(t){
        dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
      }
     
      start_time <- Sys.time()
      
      for(i in 1:n){
        for(j in 1:m){
          process1 <- Process1[i,j][[1]] 
          process2 <- Process2[i,j][[1]] 
          
          if(!(is.null(process1)|is.null(process2))){
            cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
            cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
            
            tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d), dims = c(length(process1), ngrid))
            tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d), dims = c(length(process2), ngrid))
            
            tmp1 <- apply(tmp.1,2,sum)
            tmp2 <- apply(tmp.2,2,sum)
            tmp3 <- outer(tmp1,tmp2,FUN="*")
            
            A <- A+tmp3
          }
        }
        if(i%%50==0) print(i)
      }
      
      rho1 <- rho2 <- 0
      for(i in 1:n){
        process1.tmp <- unlist(Process1[i,])
        process2.tmp <- unlist(Process2[i,])
        
        if(length(process1.tmp)>nlarge){
          grid.large <- seq(lbd,ubd,l=nlarge+1)
          grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
         process1 <- grid.mid
          freq1 <- c(table(cut(process1.tmp,breaks=grid.large)))
        }else{
          process1 <- process1.tmp
          freq1 <- rep(1,length(process1.tmp))
        }
        
        if(length(process2.tmp)>nlarge){
          grid.large <- seq(lbd,ubd,l=nlarge+1)
          grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
          process2 <- grid.mid
          freq2 <- c(table(cut(process2.tmp,breaks=grid.large)))
        }else{
          process2 <- process2.tmp
          freq2 <- rep(1,length(process2.tmp))
        }
        
        cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d)*freq1[cp1$i], dims = c(length(process1), ngrid))
    
        cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d)*freq2[cp2$i], dims = c(length(process2), ngrid))
        
        tmp1 <- apply(tmp.1,2,sum)
        tmp2 <- apply(tmp.2,2,sum)
        
        tmp3 <- outer(tmp1,tmp2,FUN="*")
        B <- B + tmp3
        rho1 <- rho1 + c(tmp1)/sqrt(diag(edge))
        rho2 <- rho2 + c(tmp2)/sqrt(diag(edge))
      }
      
      B <- B-A
      

      for(j in 1:m){
        process1.tmp <- unlist(Process1[,j])
        process2.tmp <- unlist(Process2[,j])
        
        if(length(process1.tmp)>nlarge){
          grid.large <- seq(lbd,ubd,l=nlarge+1)
          grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
          process1 <- grid.mid
          freq1 <- c(table(cut(process1.tmp,breaks=grid.large)))
        }else{
          process1 <- process1.tmp
          freq1 <- rep(1,length(process1.tmp))
        }
        
        if(length(process2.tmp)>nlarge){
          grid.large <- seq(lbd,ubd,l=nlarge+1)
          grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
          process2 <- grid.mid
          freq2 <- c(table(cut(process2.tmp,breaks=grid.large)))
        }else{
          process2 <- process2.tmp
          freq2 <- rep(1,length(process2.tmp))
        }
        
        cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d)*freq1[cp1$i], dims = c(length(process1), ngrid))
        
        cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
        tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d)*freq2[cp2$i], dims = c(length(process2), ngrid))
        
        tmp1 <- apply(tmp.1,2,sum)
        tmp2 <- apply(tmp.2,2,sum)
        
        tmp3 <- outer(tmp1,tmp2,FUN="*")
        
        C <- C + tmp3
      }
      
      C <- C-A
      
      rho1 <- rho1/(m*n)
      rho2 <- rho2/(m*n)
      
      A <- as.matrix(A/(n*m*edge))
      D <- outer(rho1,rho2,FUN="*")
      B <- as.matrix(B/((n*m*(m-1))*edge))
      C <- as.matrix(C/((n*m*(n-1))*edge))
      end_time <- Sys.time()
      
    time <- difftime(end_time,start_time,units="secs")

    Q_x <- log(B)-log(D)
    Q_y <- log(C)-log(D)
    Q_z <- log(A)+log(D)-log(B)-log(C)
    
list(Q_x=Q_x,Q_y=Q_y,Q_z=Q_z)
}
  
  
  ###Estimate the cross covariance function
  MFPCA_CrossCov_parallel <- function(Process1,Process2,lbd,ubd,bwd,ngrid,kern="epanechnikov",nlarge=10^(10),numCores=6){
    n <- dim(Process1)[1]  
    m <- dim(Process1)[2]  
    grids <- seq(lbd,ubd,l=ngrid)
    grids_mid <- (grids[-1]+grids[-ngrid])/2
    
    A <- B <- C <- D <- matrix(0,ngrid,ngrid)
    
    edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
    edge <- outer(edge,edge,FUN="*")
    
    Kh <- function(t){
      dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
    }
    
    cl <- makeCluster(numCores)
    #  registerDoParallel(cl)    
    registerDoSNOW(cl)
    ###Create progress report
    iterations <- n
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    result <- foreach(i = 1:n,.options.snow = opts)%dopar%{
      library(spatstat)
      source("MFPCA_final.R")
      tmp.A <- 0
      for(j in 1:m){
        process1 <- Process1[i,j][[1]] 
        process2 <- Process2[i,j][[1]] 
        
        if(!(is.null(process1)|is.null(process2))){
          cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
          cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
          
          tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d), dims = c(length(process1), ngrid))
          tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d), dims = c(length(process2), ngrid))
          
          tmp1 <- apply(tmp.1,2,sum)
          tmp2 <- apply(tmp.2,2,sum)
          tmp3 <- outer(tmp1,tmp2,FUN="*")
          
          tmp.A <- tmp.A+tmp3
        }
      }
      tmp.A
    }
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  ##tag1 
  print("Step 1 finished")
  
  A <- 0
  for(i in 1:n){
    A <- A + result[[i]]
  }
  rm(result)
  
  cl <- makeCluster(numCores)
  #  registerDoParallel(cl)    
  registerDoSNOW(cl)
  ###Create progress report
  iterations <- n
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  result <-  foreach(i = 1:n,.options.snow = opts)%dopar%{
    library(spatstat)
    source("MFPCA_final.R")
      tmp.B <- 0
      tmp.rho1 <- tmp.rho2 <- 0
      process1.tmp <- unlist(Process1[i,])
      process2.tmp <- unlist(Process2[i,])
      
      if(length(process1.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process1 <- grid.mid
        freq1 <- c(table(cut(process1.tmp,breaks=grid.large)))
      }else{
        process1 <- process1.tmp
        freq1 <- rep(1,length(process1.tmp))
      }
      
      if(length(process2.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process2 <- grid.mid
        freq2 <- c(table(cut(process2.tmp,breaks=grid.large)))
      }else{
        process2 <- process2.tmp
        freq2 <- rep(1,length(process2.tmp))
      }
      
      cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d)*freq1[cp1$i], dims = c(length(process1), ngrid))
      
      cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d)*freq2[cp2$i], dims = c(length(process2), ngrid))
      
      tmp1 <- apply(tmp.1,2,sum)
      tmp2 <- apply(tmp.2,2,sum)
      
      tmp3 <- outer(tmp1,tmp2,FUN="*")
      tmp.B <- tmp.B + tmp3
      tmp.rho1 <- tmp.rho1 + c(tmp1)/sqrt(diag(edge))
      tmp.rho2 <- tmp.rho2 + c(tmp2)/sqrt(diag(edge))
      list(tmp.B=tmp.B,tmp.rho1=tmp.rho1,tmp.rho2=tmp.rho2)
    }
    
  close(pb)
  #Stop the cluster
  stopCluster(cl)
  print("Step 2 finished")
  
  B <- rho1 <- rho2 <- 0
  for(i in 1:n){
    B <- B + result[[i]]$tmp.B
    rho1 <- rho1 + result[[i]]$tmp.rho1
    rho2 <- rho2 + result[[i]]$tmp.rho2
  }
  rm(result)
  
    B <- B-A
    
    cl <- makeCluster(numCores)
    #  registerDoParallel(cl)    
    registerDoSNOW(cl)
    ###Create progress report
    iterations <- m
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    
    result <-  foreach(j = 1:m,.options.snow = opts)%dopar%{
      library(spatstat)
      source("MFPCA_final.R")
      tmp.C <- 0
      process1.tmp <- unlist(Process1[,j])
      process2.tmp <- unlist(Process2[,j])
      
      if(length(process1.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process1 <- grid.mid
        freq1 <- c(table(cut(process1.tmp,breaks=grid.large)))
      }else{
        process1 <- process1.tmp
        freq1 <- rep(1,length(process1.tmp))
      }
      
      if(length(process2.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process2 <- grid.mid
        freq2 <- c(table(cut(process2.tmp,breaks=grid.large)))
      }else{
        process2 <- process2.tmp
        freq2 <- rep(1,length(process2.tmp))
      }
      
      cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d)*freq1[cp1$i], dims = c(length(process1), ngrid))
      
      cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d)*freq2[cp2$i], dims = c(length(process2), ngrid))
      
      tmp1 <- apply(tmp.1,2,sum)
      tmp2 <- apply(tmp.2,2,sum)
      
      tmp3 <- outer(tmp1,tmp2,FUN="*")
      
      tmp.C <- tmp.C + tmp3
      tmp.C
    }
    
    close(pb)
    #Stop the cluster
    stopCluster(cl)
    print("Step 3 finished")
    
    C <- 0
    for(j in 1:m){
      C <- C + result[[j]]
    }
    
    rm(result)
    
    C <- C-A
    
    
    process1.tmp <- unlist(Process1)
    process2.tmp <- unlist(Process2)
    
    if(length(process1.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process1 <- grid.mid
      freq1 <- c(table(cut(process1.tmp,breaks=grid.large)))
    }else{
      process1 <- process1.tmp
      freq1 <- rep(1,length(process1.tmp))
    }
    
    if(length(process2.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process2 <- grid.mid
      freq2 <- c(table(cut(process2.tmp,breaks=grid.large)))
    }else{
      process2 <- process2.tmp
      freq2 <- rep(1,length(process2.tmp))
    }
    
    cp1 <- crosspairs(ppp(x=process1,y=process1*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp.1 <- Matrix::sparseMatrix(i = cp1$i, j = cp1$j, x = Kh(cp1$d)*freq1[cp1$i], dims = c(length(process1), ngrid))
    
    cp2 <- crosspairs(ppp(x=process2,y=process2*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp.2 <- Matrix::sparseMatrix(i = cp2$i, j = cp2$j, x = Kh(cp2$d)*freq2[cp2$i], dims = c(length(process2), ngrid))
    
    tmp1 <- apply(tmp.1,2,sum)
    tmp2 <- apply(tmp.2,2,sum)
    
    tmp3 <- outer(tmp1,tmp2,FUN="*")
    D <- (tmp3-B-C+A)/(m*n*(m-1)*(n-1)*edge)
    
    rho1 <- rho1/(m*n)
    rho2 <- rho2/(m*n)
    
    A <- as.matrix(A/(n*m*edge))
#    D <- outer(rho1,rho2,FUN="*")
    B <- as.matrix(B/((n*m*(m-1))*edge))
    C <- as.matrix(C/((n*m*(n-1))*edge))

    Q_x <- log(B)-log(D)
    Q_y <- log(C)-log(D)
    Q_z <- log(A)+log(D)-log(B)-log(C)
    
    list(Q_x=Q_x,Q_y=Q_y,Q_z=Q_z)
  }
  
  ####Measure the KL distance between observed process and a target intensity
  KLdist <- function(lbd,ubd,Process,r,Lambda){
    n <- length(Process)
    ngrid <- dim(Lambda)[2]
    grids <- seq(lbd,ubd,l=ngrid)
    
    process_all <- unlist(Process)
    n.r <- length(r)
    Count_emp <- matrix(NA,n,n.r-1)
    
    for(i in 1:n){
      process_sel <- Process[[i]]
      Count_emp[i,] <- table(cut(process_sel,breaks = r))
    }
    
    counts_all <- colSums(Count_emp)
    Prob.emp <- scale(Count_emp,center=FALSE,scale=counts_all)
    
      lambda_all <- colSums(Lambda)
      
      Prob_est <-  matrix(NA,n,n.r-1)
      
      for(i in 1:n){
        process <- Process[[i]]
        for(j in 1:(n.r-1)){
          ind <- (grids>=r[j])&(grids<=r[j+1])
          Prob_est[i,j] <- sum(Lambda[i,ind])/sum(lambda_all[ind])
        }
      }
      Prob_est <- scale(Prob_est,center=FALSE,scale=colSums(Prob_est))
      val <- -colSums(Prob.emp*(log(Prob_est)-log(Prob.emp+10^(-100))))
      return(val)
  }
  
  ####Make sign change of the Eigen functions
  SignMatch <- function(Eigen1,Eigen2){
    p <- dim(Eigen1)[2]
    Sign.change <- numeric()
    for(i in 1:p){
      Sign.change[i] <- ifelse(mean(abs(Eigen1[,i]-Eigen2[,i]))<mean(abs(Eigen1[,i]+Eigen2[,i])),1,-1)
    }
    Sign.change
  }
  
  
  ###Choose the bandwidth using leave-a-block-out cross validation
  MFPCA_parallel <- function(Process,lbd,ubd,bwd,ngrid,kern="epanechnikov",nlarge=10^(10),npc.x=2,npc.y=2,numCores=6){
    n <- dim(Process)[1]  
    m <- dim(Process)[2]  
    grids <- seq(lbd,ubd,l=ngrid)
    grids_mid <- (grids[-1]+grids[-ngrid])/2
    
    A2 <- A <- B <- C <- D <- matrix(0,ngrid,ngrid)
    
    edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
    edge <- outer(edge,edge,FUN="*")
    
    Kh <- function(t){
      dkernel(t/bwd,kernel = kern,sd=sqrt(1/5))/bwd
    }
    
    
    cl <- makeCluster(numCores)
    #  registerDoParallel(cl)    
    registerDoSNOW(cl)
    ###Create progress report
    iterations <- n
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    result <- foreach(i = 1:n,.options.snow = opts)%dopar%{
      library(spatstat)
      source("MFPCA_final.R")
      tmp.A <- tmp.A2 <- 0
      for(j in 1:m){
        process <- Process[i,j][[1]] 
        if(!is.null(process)){
          cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
          tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
          tmp1 <- apply(tmp,2,sum)
          tmp2 <- outer(tmp1,tmp1,FUN="*")
          tmp.A <- tmp.A+(tmp2-Matrix::t(tmp)%*%tmp)
          tmp.A2 <- tmp.A2 + tmp2
        }
      }
      list(tmp.A=tmp.A,tmp.A2=tmp.A2)
    }
    close(pb)
    #Stop the cluster
    stopCluster(cl)
    ##tag1 
    print("Step 1 finished")
    
    A <- A2 <- 0
    
    for(i in 1:n){
      A <- A + result[[i]]$tmp.A
      A2 <- A2 + result[[i]]$tmp.A2
    }
    
    rm(result)
    
    
    cl <- makeCluster(numCores)
    #  registerDoParallel(cl)    
    registerDoSNOW(cl)
    ###Create progress report
    iterations <- n
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    
    result <-  foreach(i = 1:n,.options.snow = opts)%dopar%{
      library(spatstat)
      source("MFPCA_final.R")
      process.tmp <- unlist(Process[i,])
      if(length(process.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process <- grid.mid
        freq <- c(table(cut(process.tmp,breaks=grid.large)))
      }else{
        process <- process.tmp
        freq <- rep(1,length(process.tmp))
      }
      cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
      tmp1 <- apply(tmp,2,sum)
      tmp2 <- outer(tmp1,tmp1,FUN="*")
      list(tmp.B=tmp2,tmp.rho=c(tmp1)/sqrt(diag(edge)))
    }
    close(pb)
    #Stop the cluster
    stopCluster(cl)
    print("Step 2 finished")
    
    B <- rho <- 0
    for(i in 1:n){
      B <- B + result[[i]]$tmp.B
      rho <- rho + result[[i]]$tmp.rho/m
    }
    
    B <- B-A2
    
    rm(result)
    
    
    cl <- makeCluster(numCores)
    #  registerDoParallel(cl)    
    registerDoSNOW(cl)
    ###Create progress report
    iterations <- m
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    
    result <-  foreach(j = 1:m,.options.snow = opts)%dopar%{
      library(spatstat)
      source("MFPCA_final.R")
      process.tmp <- unlist(Process[,j])
      if(length(process.tmp)>nlarge){
        grid.large <- seq(lbd,ubd,l=nlarge+1)
        grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
        process <- grid.mid
        freq <- c(table(cut(process.tmp,breaks=grid.large)))
      }else{
        process <- process.tmp
        freq <- rep(1,length(process.tmp))
      }
      cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
      tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
      tmp1 <- apply(tmp,2,sum)
      tmp2 <- outer(tmp1,tmp1,FUN="*")
      list(tmp.C=tmp2)
    }
    close(pb)
    #Stop the cluster
    stopCluster(cl)
    print("Step 3 finished")
    
    C <- 0
    for(j in 1:m){
      C <- C + result[[j]]$tmp.C
    }
    
    C <- C-A2
    rm(result)
    
    
    process.tmp <- unlist(Process)
    if(length(process.tmp)>nlarge){
      grid.large <- seq(lbd,ubd,l=nlarge+1)
      grid.mid <- (grid.large[-1]+grid.large[-(nlarge+1)])/2
      process <- grid.mid
      freq <- c(table(cut(process.tmp,breaks=grid.large)))
    }else{
      process <- process.tmp
      freq <- rep(1,length(process.tmp))
    }
    cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
    tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d)*freq[cp$i], dims = c(length(process), ngrid))
    tmp1 <- apply(tmp,2,sum)
    tmp2 <- outer(tmp1,tmp1,FUN="*")

    D <- (tmp2-B-C+A2)/(m*n*(m-1)*(n-1)*edge)
    
    A <- as.matrix(A/((n*m-1)*edge))
    #D <- outer(rho/(n),rho/(n),FUN="*")
    
    B <- as.matrix(B/((n*m*(m-1))*edge))
    C <- as.matrix(C/((n*m*(n-1))*edge))
    
    rho <- rho/n
    
    Cov.x <- as.matrix(log(B)-log(D))
    
    Cov.y <- as.matrix(log(C)-log(D))
    
    Cov.z <- log(A)+log(D)-log(B)-log(C)
    
    
    ###Find Eigen functions and eigen values
    obj.x <- eigen(Cov.x)
    Eigen.x <- obj.x$vectors/sqrt(grids[2]-grids[1])
    eta.x <- obj.x$values*(grids[2]-grids[1])
    obj.y <- eigen(Cov.y)
    Eigen.y <- obj.y$vectors/sqrt(grids[2]-grids[1])
    eta.y <- obj.y$values*(grids[2]-grids[1])
    
    obj.z <- eigen(Cov.z)
    Eigen.z <- obj.z$vectors/sqrt(grids[2]-grids[1])
    eta.z <- obj.z$values*(grids[2]-grids[1])
    
    sig2.x <- diag(Cov.x)
    sig2.y <- diag(Cov.y)
    
    list(A=A,B=B,C=C,D=D,rho=rho,Cov.x=Cov.x,Cov.y=Cov.y,Cov.z=Cov.z,Eigen.x=Eigen.x,Eigen.y=Eigen.y,Eigen.z=Eigen.z,eta.x=eta.x,eta.y=eta.y,eta.z=eta.z)
  }
  
  
  ####Estimate single level model using Wu et al. 2013
 FPCA_single <- function(process.x,lbd,ubd,ngrid,bwd,kern="epanechnikov",n.pc=2){
   
   n <- length(process.x)
   N <- sapply(process.x,length)
   grids <- seq(lbd,ubd,l=ngrid)
   ########Estimate the covariance kernel
   Kh <- function(t){
     dkernel(t/bwd,kernel =kern,sd=sqrt(1/5))/bwd
   }
   
   edge <- pkernel((grids-lbd)/bwd,kernel = kern,sd=sqrt(1/5))-pkernel((grids-ubd)/bwd,kernel = kern,sd=sqrt(1/5))
   edge <- outer(edge,edge,FUN="*")
   
   process_pool <- unlist(process.x)
   cp <- crosspairs(ppp(x=process_pool,y=process_pool*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
   tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process_pool), ngrid))
   f.mu <- apply(tmp,2,sum)/sum(N)/sqrt(diag(edge))
   
   A <- 0
   for(i in 1:n){
     process <- unlist(process.x[[i]])
     cp <- crosspairs(ppp(x=process,y=process*0),ppp(x=grids,y=grids*0),rmax = bwd,what = "ijd")
     tmp <- Matrix::sparseMatrix(i = cp$i, j = cp$j, x = Kh(cp$d), dims = c(length(process), ngrid))
     tmp1 <- apply(tmp,2,sum)
     tmp2 <- outer(tmp1,tmp1,FUN="*")
     A <- A+(tmp2-Matrix::t(tmp)%*%tmp)
   }
   Ghat <- A/sum(N*(N-1))/edge-outer(f.mu,f.mu,FUN="*")
   Eigen <- eigen(Ghat)$vectors/sqrt(grids[2]-grids[1])
   Eigen_mid <- (Eigen[-1,,drop=FALSE]+Eigen[-ngrid,,drop=FALSE])/2
   f.mu_mid <- (f.mu[-1]+f.mu[-ngrid])/2
   xi <- matrix(NA,n,n.pc)
   for(k in 1:n.pc){
     phi_mid <- Eigen_mid[,k]
     for(i in 1:n){
     freq <- table(cut(process.x[[i]],breaks=grids))
     xi[i,k] <- sum(phi_mid*freq)/sum(freq)-sum(f.mu_mid*phi_mid*diff(grids))
     }
   }
   
   intensity.hat <- matrix(NA,n,ngrid)
   
   for(i in 1:n){
     tmp <- f.mu
     for(k in 1:n.pc)
       tmp <- tmp + xi[i,k]*Eigen[,k]
     ind <- tmp<0
     tmp[ind] <- 0
     tmp <- tmp/sum(tmp*(grids[2]-grids[1]))
     intensity.hat[i,] <- tmp * N[i]
   }

   list(Eigen=Eigen,intensity.hat=intensity.hat,xi=xi,f.mu=f.mu,Ghat=Ghat)
 }
