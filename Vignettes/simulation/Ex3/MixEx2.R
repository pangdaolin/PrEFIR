library(MASS)
library(mgcv)
library(pROC)
library(parallel)
library(rmutil)
library(GFM)

source("~/pdl/project3_simu/eval_space.R")
source("~/pdl/project3_simu/gmf_p.R")
source("~/pdl/project3_simu/GIRL_h_parallel_1.r")

for (n in c(100,200,600)) {
  for (p in c(50,100,300)) {
    d=2; k=1; mixp=p/2
    
    mAUC_GFM <- numeric()
    
    alpha <- rep(0,p)#rnorm(p, mean = 0, sd = 1)rep(0,p)
    Phi <- rbind(matrix(c(-1,1,-1,1,-1,1,1,1,-1,1),1,p),matrix(c(0,0,0,0,1,0,0,1,1,1),1,p))
    Gamma <- matrix(c(1,1,1,1,0,-1,1,0,0,0),k,p)
    sigma0 <- rep(1,mixp)
    
    system.time(for (runs in 1:100) {
      set.seed(927017+runs)
      f <- runif(n*d,-2,2)
      ff <- matrix(f,n,d)
      Y <- rnorm(n, mean = 0, sd = 1)
      V = as.numeric(Y>0)
      stru1 <- (matrix(alpha, n, p, byrow = TRUE) + ff %*% (Phi) + V %*% (Gamma))
      stru1.norm <- stru1[,1:mixp]
      stru1.posi <- stru1[,-(1:mixp)]
      
      
      X <- matrix(0,n,p)
      for (i in 1:n) {
        while (sum(X[i,(1+mixp):(p)])==0) {
          for (j in 1:mixp) {
            X[i,j] <- rnorm(1, mean = stru1[i,j], sd = sigma0[j])#rdoublepois(1, m = stru1[i,], s = theta)
          }
          for (j in (1+mixp):(p)) {
            X[i,j] <- rpois(1, lambda = exp(stru1[i,j]))#rdoublepois(1, m = stru1[i,], s = theta)
          }
        }
      }
      
      X.train <- X[1:(n/2),];V.train <- V[1:(n/2)];ff.train <- ff[1:(n/2),]
      X.test <- X[(n/2+1):n,];V.test <- V[(n/2+1):n];ff.test <- ff[(n/2+1):n,]
      Y.train <- V.train
      Y.test<- V.test
      
      #############################
      X.norm <- X[,1:mixp]
      X.pois <- X[,(1+mixp):(p)]
      X.all_list <- list(X.norm,X.pois)
      gfm2 <- gfm(X.all_list, types=c('gaussian','poisson'), q=3, verbose = FALSE, maxIter=5)
      F.train <- gfm2[["hH"]][1:(n/2),]
      F.test <- gfm2[["hH"]][(n/2+1):n,]
      
      train.object <- data.frame(X=F.train,Y=Y.train)
      test.object <- data.frame(X=F.test)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2),data=train.object,family = binomial()),silent = T)
      pre1 <- predict(GAM1,test.object,type = "response")
      AUC1 <- roc(as.factor(Y.test),pre1)[["auc"]]
      mAUC_GFM <- c(mAUC_GFM,AUC1)   
      
      cat("Ber n =",n/2,"p =",p,"Iteration step",runs,"completed","\n")
      
    })
    save.image(paste0("~/pdl/review_p2/biometrics/simu/Mixd2/GFMMixEx2n",n/2,"p",p,".RData"))
  }
}