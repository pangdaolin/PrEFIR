library(MASS)
library(mgcv)
library(parallel)
library(randomForest)
library(pROC)
library(MAVE)
library(rmutil)
library(GFM)

source("~/pdl/project3_simu/eval_space.R")
source("~/pdl/review_p2/biometrics/gmf_mix_1.R")
source("~/pdl/project3_simu/GIRL_h_parallel_1.r")

for (n in c(100,200,600)) {
  for (p in c(50,100,300)) {
    # n=60;p=30
    d=2; k=1; mixp=p/2
    
    eval1 <- c(1,1,1)
    eval2 <- c(1,1,1)
    eval3 <- c(1,1,1)
    eval4 <- c(1,1,1)
    eval5 <- c(1,1,1)
    eval6 <- c(1,1,1)
    eval7 <- c(1,1,1)
    
    
    rmse1 <- numeric()
    rmse2 <- numeric()
    rmse3 <- numeric()
    rmse4 <- numeric()
    rmse5 <- numeric()
    rmse6 <- numeric()
    rmse7 <- numeric()
    rmse8 <- numeric()
    rmse9 <- numeric()
    rmse10 <- numeric()
    rmse1_0 <- numeric()
    rmse2_0 <- numeric()
    rmse3_0 <- numeric()
    rmse4_0 <- numeric()
    rmse6_0 <- numeric()
    rmse8_0 <- numeric()
    rmse9_0 <- numeric()
    errlist1 <- numeric()
    errlist2 <- numeric()
    errlist3 <- numeric()
    
    rmse7 <- numeric()
    rmse7_0 <- numeric()
    rmse7_1 <- numeric()   
    
    rmse_GFM <- numeric()
    
    errlist7 <- numeric()
    
    alpha <- rep(0,p)#rnorm(p, mean = 0, sd = 1)rep(0,p)
    Phi <- rbind(matrix(c(-1,1,-1,1,-1,1,1,1,-1,1),1,p),matrix(c(0,0,0,0,1,0,0,1,1,1),1,p))
    Gamma <- matrix(c(1,1,1,1,0,-1,1,0,0,0),k,p)
    sigma0 <- rep(1,mixp)
    
    
    
    
    system.time(for (runs in 1:100) {
      set.seed(927018+runs)
      
      
      f <- rnorm(n*d,mean = 0, sd = 1)
      ff <- matrix(f,n,d)
      # Y <- rnorm(n, mean = 0, sd = 1)
      # V = Y#as.numeric(Y>0)
      stru1 <- (matrix(alpha, n, p, byrow = TRUE) + ff %*% (Phi) )
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
      
      X.train <- X[1:(n/2),];ff.train <- ff[1:(n/2),]
      X.test <- X[(n/2+1):n,];ff.test <- ff[(n/2+1):n,]
      X.norm <- X[,1:mixp]
      X.pois <- X[,(1+mixp):(p)]
      
      b_1=c(1,-2,3)
      b_2=c(1,1,1,1,0,-1,1,0,0,0)
      b_1=c(b_1,rep(0,mixp-length(b_1)))
      b_2=c(b_2,rep(0,mixp-length(b_2)))
      
      Y = X.norm %*% b_1 + exp(10*(X.pois %*% b_2)/rowSums(X.pois)) + rnorm(n, mean = 0, sd = 1)
      
      Y.train <- Y[1:(n/2)]
      Y.test<- Y[(n/2+1):n]
      V.train <- Y.train
      V.test <- Y.test
      #############################
      fit1 <- try(gmf_mix(Y=as.matrix(X.train),X=as.matrix(Y.train),family=poisson(),
                          mixp=p/2,p=d,penaltyU=1,penaltyV=0,penaltyBeta=0,
                          intercept=T,
                          init="svd",
                          normalize_uv=F,
                          verbose=F,
                          method="airwls",
                          airwls.internal.steps = 1,
                          parallel=1), silent = TRUE)
      
      
      if("try-error" %in% class(fit1)){ 
        errlist1 <- c(errlist1,runs)
        print("error")
        print(fit1)
      }else{
        Gamma1 <- as.vector(fit1[["beta"]][-1,])
        Beta1 <- (fit1[["v"]])
        coef1 <- t(cbind(Gamma1,Beta1))
        rownames(coef1) <- NULL
      }
      
      dr.norm <- solve(Beta1[1:mixp,]%*%t(Beta1[1:mixp,])+diag(mixp))%*%Gamma1[1:mixp]
      # dr.norm <- t(t(dr.norm)/sqrt(colSums(dr.norm^2)))*sqrt(nrow(dr.norm))
      z1_train.norm <- X.train[,1:mixp]%*%dr.norm
      z1_test.norm <- X.test[,1:mixp]%*%dr.norm
      
      dr.posi <- coef1[,-(1:mixp)]
      # dr.posi <- t(t(dr.posi)/sqrt(colSums(dr.posi^2)))*sqrt(nrow(dr.posi))
      
      z1_train.posi <- X.train[,-(1:mixp)]%*%t(dr.posi)
      z1_test.posi <- X.test[,-(1:mixp)]%*%t(dr.posi)
      z1_train <- cbind(z1_train.norm,z1_train.posi)
      z1_test <- cbind(z1_test.norm,z1_test.posi)
      
      
      
      train.object <- data.frame(X=z1_train,Y=Y.train)
      test.object <- data.frame(X=z1_test)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2)+s(X.3)+s(X.4,k=3),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse1 <- c(rmse1,NA)
      }else{
        pre1 <- predict(GAM1,test.object)
        rmse1 <- c(rmse1,sqrt(sum((pre1-Y.test)^2)/length(Y.test)))
      }
      
      GAM1 <- try(gam(Y~s(X.1,X.2,X.3,X.4,k=16),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse1_0 <- c(rmse1_0,NA)
      }else{
        pre1_0 <- predict(GAM1,test.object)
        rmse1_0 <- c(rmse1_0,sqrt(sum((pre1_0-Y.test)^2)/length(Y.test)))
      }
      # sqrt(sum((pre1-Y.test)^2)/length(Y.test))
      # sqrt(sum((pre1_0-Y.test)^2)/length(Y.test))
      
      # ####################################
      fit2 <- try(gmf_mix(Y=as.matrix(X.train),X=as.matrix(Y.train),family=poisson(),
                          mixp=p/2,p=d,penaltyU=1,penaltyV=0,penaltyBeta=0,
                          intercept=T,
                          init="svd",
                          normalize_uv=F,
                          verbose=F,
                          method="quasi",
                          airwls.internal.steps = 1,
                          parallel=1), silent = TRUE)
      
      
      if("try-error" %in% class(fit2)){
        errlist2 <- c(errlist2,runs)
        print("error")
      }else{
        Gamma2 <- as.vector(fit2[["beta"]][-1,])
        Beta2 <- (fit2[["v"]])
        coef2 <- t(cbind(Gamma2,Beta2))
        rownames(coef2) <- NULL
      }
      
      
      dr.norm <- solve(Beta2[1:mixp,]%*%t(Beta2[1:mixp,])+diag(mixp))%*%Gamma2[1:mixp]
      # dr.norm <- t(t(dr.norm)/sqrt(colSums(dr.norm^2)))*sqrt(nrow(dr.norm))
      z1_train.norm <- X.train[,1:mixp]%*%dr.norm
      z1_test.norm <- X.test[,1:mixp]%*%dr.norm
      
      dr.posi <- coef2[,-(1:mixp)]
      # dr.posi <- t(t(dr.posi)/sqrt(colSums(dr.posi^2)))*sqrt(nrow(dr.posi))
      
      z1_train.posi <- X.train[,-(1:mixp)]%*%t(dr.posi)
      z1_test.posi <- X.test[,-(1:mixp)]%*%t(dr.posi)
      z1_train <- cbind(z1_train.norm,z1_train.posi)
      z1_test <- cbind(z1_test.norm,z1_test.posi)
      
      
      
      train.object <- data.frame(X=z1_train,Y=Y.train)
      test.object <- data.frame(X=z1_test)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2)+s(X.3)+s(X.4,k=3),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse2 <- c(rmse2,NA)
      }else{
        pre2 <- predict(GAM1,test.object)
        rmse2 <- c(rmse2,sqrt(sum((pre2-Y.test)^2)/length(Y.test)))
      }
      
      GAM1 <- try(gam(Y~s(X.1,X.2,X.3,X.4,k=12),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse2_0 <- c(rmse2_0,NA)
      }else{
        pre2_0 <- predict(GAM1,test.object)
        rmse2_0 <- c(rmse2_0,sqrt(sum((pre2_0-Y.test)^2)/length(Y.test)))
      }
      # sqrt(sum((pre2-Y.test)^2)/length(Y.test))
      # sqrt(sum((pre2_0-Y.test)^2)/length(Y.test))
      
      ##############
      coef3 <- numeric()
      for (otus in 1:mixp) {
        fit3 <- glm(X.train[,otus]~V.train,family = gaussian())
        coef3[otus] <- coef(fit3)[-1]
      }
      for (otus in (1+mixp):p) {
        fit3 <- glm(X.train[,otus]~V.train,family = poisson())
        coef3[otus] <- coef(fit3)[-1]
      }
      
      # z2_train <- X.train%*%(coef3)#[1:mixp]
      # z2_test <- X.test%*%(coef3)
      # res1 = data.frame(Y=Y.train,X1=(z2_train))#,m=scale(m)
      # # res21 = data.frame(Y=Y.test,X1=(z11_test),X2=(z12_test))
      # train.object <- data.frame(X1=z2_train)
      # test.object <- data.frame(X1=z2_test)
      # Z2 <- as.matrix(train.object)
      # GAM1 <- gam(Y~s(X1),data=data.frame(res1))
      # # GAM1 <- lm(Y~X1,data=data.frame(res1))
      # pre3 <- predict(GAM1,test.object)
      
      z2_train.norm <- X.train[,1:mixp]%*%(coef3[1:mixp])
      z2_test.norm <- X.test[,1:mixp]%*%(coef3[1:mixp])
      z2_train.posi <- X.train[,-(1:mixp)]%*%(coef3[-(1:mixp)])
      z2_test.posi <- X.test[,-(1:mixp)]%*%(coef3[-(1:mixp)])
      z2_train=cbind(z2_train.norm,z2_train.posi)
      z2_test=cbind(z2_test.norm,z2_test.posi)
      res1 = data.frame(Y=Y.train,X=(z2_train))#,m=scale(m)
      # res21 = data.frame(Y=Y.test,X1=(z11_test),X2=(z12_test))
      train.object <- data.frame(X=z2_train)
      test.object <- data.frame(X=z2_test)
      Z2 <- as.matrix(train.object)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2),data=data.frame(res1)),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse3 <- c(rmse3,NA)
      }else{
        pre3 <- predict(GAM1,test.object)
        rmse3 <- c(rmse3,sqrt(sum((pre3-Y.test)^2)/length(Y.test)))
      }
      GAM1 <- try(gam(Y~s(X.1,X.2),data=data.frame(res1)),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse3_0 <- c(rmse3_0,NA)
      }else{
        pre3_0 <- predict(GAM1,test.object)
        rmse3_0 <- c(rmse3_0,sqrt(sum((pre3_0-Y.test)^2)/length(Y.test)))
      }
      # sqrt(sum((pre3-Y.test)^2)/length(Y.test))
      # sqrt(sum((pre3_0-Y.test)^2)/length(Y.test))
      
      ##############
      coef4 <- matrix(0,d+1,p)
      for (otus in 1:mixp) {
        fit4 <- glm(X.train[,otus]~V.train+ff.train,family = gaussian())
        coef4[,otus] <- coef(fit4)[-1]
      }
      for (otus in (1+mixp):p) {
        fit4 <- glm(X.train[,otus]~V.train+ff.train,family = poisson())
        coef4[,otus] <- coef(fit4)[-1]
      }
      
      dr.norm <- solve(coef4[2,1:mixp]%*%t(coef4[2,1:mixp])+diag(mixp))%*%coef4[1,1:mixp]
      # dr.norm <- t(t(dr.norm)/sqrt(colSums(dr.norm^2)))*sqrt(nrow(dr.norm))
      # dr.norm <- c(1,1,1,1,0)
      z1_train.norm <- X.train[,1:mixp]%*%dr.norm
      z1_test.norm <- X.test[,1:mixp]%*%dr.norm
      # dr.posi <- t(spantrue[-(1:mixp),])
      
      dr.posi <- coef4[,-(1:mixp)]
      # dr.posi <- t(t(dr.posi)/sqrt(colSums(dr.posi^2)))*sqrt(nrow(dr.posi))
      
      z1_train.posi <- X.train[,-(1:mixp)]%*%t(dr.posi)
      z1_test.posi <- X.test[,-(1:mixp)]%*%t(dr.posi)
      z1_train <- cbind(z1_train.norm,z1_train.posi)
      z1_test <- cbind(z1_test.norm,z1_test.posi)
      
      
      
      train.object <- data.frame(X=z1_train,Y=Y.train)
      test.object <- data.frame(X=z1_test)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2)+s(X.3)+s(X.4,k=3),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse4 <- c(rmse4,NA)
      }else{
        pre4 <- predict(GAM1,test.object)
        rmse4 <- c(rmse4,sqrt(sum((pre4-Y.test)^2)/length(Y.test)))
      }
      
      GAM1 <- try(gam(Y~s(X.1,X.2,X.3,X.4,k=12),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse4_0 <- c(rmse4_0,NA)
      }else{
        pre4_0 <- predict(GAM1,test.object)
        rmse4_0 <- c(rmse4_0,sqrt(sum((pre4_0-Y.test)^2)/length(Y.test)))
      }
      # sqrt(sum((pre4-Y.test)^2)/length(Y.test))
      # sqrt(sum((pre4_0-Y.test)^2)/length(Y.test))
      
      #####################################
      fit5_1 <- try(Girl_h(y=V.train, X=X.train[,1:mixp], G = 0,family="gaussian", method = "BIC", sparseGamma = F, parallel = F),silent=T)
      fit5_2 <- try(Girl_h(y=V.train, X=X.train[,(1+mixp):p], G = 0,family="poisson", method = "BIC", sparseGamma = F, parallel = F),silent=T)
      
      Gamma5 <- c(c(fit5_1[["Gamma.PGM"]]),c(fit5_2[["Gamma.PGM"]]))
      coef5 <- Gamma5
      z5_train <- X.train%*%(Gamma5)
      z5_test <- X.test%*%(Gamma5)
      train.object <- data.frame(X=z5_train,Y=Y.train)
      test.object <- data.frame(X=z5_test)
      GAM1 <- try(gam(Y~s(X),data=train.object),silent = T)
      if("try-error" %in% class(GAM1)){ 
        rmse5 <- c(rmse5,NA)
      }else{
        pre5 <- predict(GAM1,test.object)
        rmse5 <- c(rmse5,sqrt(sum((pre5-Y.test)^2)/length(Y.test)))
      }
      
      # # sqrt(sum((pre10-Y.test)^2)/length(Y.test))
      
      
      
      #######################
      fit7=mave(Y.train~X.train,max.dim=d+1)
      
      if(is.na(fit7[["dir"]][[d+1]][1,1])){
        errlist7 <- c(errlist7,runs)
        rmse7 <- c(rmse7,NA)
        rmse7_0 <- c(rmse7_0,NA)
        rmse7_1 <- c(rmse7_1,NA)
        print("error7")
      }else{
        coef7 <- coef(fit7, dim=d+1)
        rownames(coef7) <- NULL
        colnames(coef7) <- NULL
        coef7 <- t(coef7)
        z7_train <- X.train%*%t(coef7)
        z7_test <- X.test%*%t(coef7)
        
        train.object <- data.frame(X=z7_train,Y=Y.train)
        test.object <- data.frame(X=z7_test)
        
        GAM1 <- gam(Y~s(X.1)+s(X.2)+s(X.3),data=train.object)
        pre7 <- predict(GAM1,test.object)
        rmse7 <- c(rmse7,sqrt(sum((pre7-Y.test)^2)/length(Y.test)))
        GAM1 <- gam(Y~s(X.1,X.2,X.3,k=10),data=train.object)
        pre7_0 <- predict(GAM1,test.object)
        rmse7_0 <- c(rmse7_0,sqrt(sum((pre7_0-Y.test)^2)/length(Y.test)))
        pre7_1 <- predict(fit7, newx=X.test, dim=d+1)
        rmse7_1 <- c(rmse7_1,sqrt(sum((pre7_1-Y.test)^2)/length(Y.test)))
        
      }
      
      # z1_train <- X.train%*%t(coef4)
      # z1_test <- X.test%*%t(coef4)
      # train.object <- data.frame(X=z1_train,Y=Y.train)
      # test.object <- data.frame(X=z1_test)
      # GAM1 <- gam(Y~s(X.1)+s(X.2),data=train.object)
      # pre4 <- predict(GAM1,test.object)
      # GAM1 <- gam(Y~s(X.1,X.2,k=12),data=train.object)
      # pre4_0 <- predict(GAM1,test.object)
      
      ########################################
      # fit10 <- Girl_h(y=V.train, X=X.train, G = 0, method = "BIC", sparseGamma = F, parallel = F)
      # Gamma10 <- fit10[["Gamma.PGM"]]
      # z7_train <- X.train%*%(Gamma10)
      # z7_test <- X.test%*%(Gamma10)
      # train.object <- data.frame(X=z7_train,Y=Y.train)
      # test.object <- data.frame(X=z7_test)
      # # Z2 <- as.matrix(train.object)
      # GAM1 <- gam(Y~s(X),data=train.object)
      # # GAM1 <- lm(Y~X1,data=data.frame(res1))
      # pre10 <- predict(GAM1,test.object)
      # # sqrt(sum((pre10-Y.test)^2)/length(Y.test))
      
      #############################
      X.norm <- X[,1:mixp]
      X.pois <- X[,(1+mixp):(p)]
      X.all_list <- list(X.norm,X.pois)
      gfm2 <- gfm(X.all_list, types=c('gaussian','poisson'), q=3, verbose = FALSE, maxIter=5)
      F.train <- gfm2[["hH"]][1:(n/2),]
      F.test <- gfm2[["hH"]][(n/2+1):n,]
      
      train.object <- data.frame(X=F.train,Y=Y.train)
      test.object <- data.frame(X=F.test)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2),data=train.object),silent = T)
      pre1 <- stats::predict(GAM1,test.object)
      rmse_GFM <- c(rmse_GFM,sqrt(sum((pre1-Y.test)^2)/length(Y.test)))
      
      GAM1 <- try(gam(Y~s(X.1,X.2,k=12),data=train.object),silent = T)
      pre1 <- stats::predict(GAM1,test.object)
      rmse_GFM_0 <- c(rmse_GFM_0,sqrt(sum((pre1-Y.test)^2)/length(Y.test)))
      
      
      cat("Iteration step",runs,"completed","\n")
      
    })
    save.image(paste0("~/pdl/review_p2/biometrics/simu/Ex4/MixEx1n",n/2,"p",p,".RData"))
  }
}



mean(rmse1)
mean(rmse1_0)
mean(rmse2)
mean(rmse2_0)
mean(rmse3)
mean(rmse3_0)
mean(rmse4)
mean(rmse4_0)
mean(rmse5)
mean(rmse7)
mean(rmse7_0)
mean(rmse7_1)
mean(rmse_GFM)
mean(rmse_GFM_0)

mean(rmse10)

sd(rmse1)
sd(rmse1_0)
sd(rmse2)
sd(rmse2_0)
sd(rmse3)
sd(rmse4)
sd(rmse4_0)
sd(rmse10)


rmsedata <- data.frame(rmse=c(rmse1,rmse2,rmse3,rmse4,rmse10),
                       methods=factor(c(rep("FADR-IRLS",runs),rep("FADR-quasi",runs),rep("EFDR",runs),rep("Oracle",runs),rep("EFDR-GM",runs)),
                                      level = c("EFDR","EFDR-GM","FADR-IRLS","FADR-quasi","Oracle")))

p1 = ggplot(rmsedata, aes(x=methods, y=rmse, color = methods)) +
  geom_boxplot()+labs(title="", x="", y="RMSE")+ylim(0.10,1.9)
p1
#####################################################
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggsci)

for (ni in c(100,200,600)) {
  for (pj in c(50,100,300)) {
    for (Exk in 1:1) {
      load(paste0("~/pdl/review_p2/biometrics/simu/Ex4/Mix","Ex",Exk,"n",ni/2,"p",pj,".RData"))
      
      if(Exk==1&pj==50&ni==100){
        dataEx <- data.frame(n=rep(n/2,5*100),p=rep(paste0("p=",p),5*100),
                             Method=c(rep("EFIR",100),rep("GraphIR",100),rep("PrEFIR",100),rep("MAVE",100),rep("GFM",100)),
                             Example=rep("Ex1",5*100),
                             Error=c((rmse3),(rmse5),(rmse2),(rmse7),rmse_GFM))
      } else{
        
          dataEx <- rbind(dataEx,data.frame(n=rep(n/2,5*100),p=rep(paste0("p=",p),5*100),
                                            Method=c(rep("EFIR",100),rep("GraphIR",100),rep("PrEFIR",100),rep("MAVE",100),rep("GFM",100)),
                                            Example=rep("Ex1",5*100),
                                            Error=c((rmse3),(rmse5),(rmse2),(rmse7),rmse_GFM)))
      }
    }
  }
}

dataEx%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("GFM","EFIR","GraphIR","MAVE","PrEFIR"))), outlier.alpha = 0.1)+#
  labs(x = "", y = "Root mean squared error")+  
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  scale_fill_discrete(limits=c("GFM","EFIR","GraphIR","MAVE","PrEFIR"))+
  scale_y_continuous(breaks = c(0.5,1.0,1.5,2.0),limits = c(0.1,20.1))+
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_text(size=14,face="bold.italic"),
        strip.background = element_rect(color="white",size=0),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom"
  ) + scale_color_d3() ->p1

p1
p1a<- p1 + theme(legend.position = "none")


dataEx%>%filter(Example=="Ex2")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("EFIR","GraphIR","MAVE","PrEFIR"))), outlier.alpha = 0.1)+#
  labs(x = expression(bolditalic(n)), y = "Area under the ROC curve")+  
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  scale_fill_discrete(limits=c("EFIR","GraphIR","MAVE","PrEFIR"))+
  scale_y_continuous(breaks = c(0.4,0.6,0.8,1.0),limits = c(0.4,1.0))+
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_blank(),
        strip.background = element_rect(color="white",size=0),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom"
  ) + scale_color_d3() ->p2

p2
p2a<- p2 + theme(legend.position = "none")



pp <- cowplot::plot_grid(p1a,p2a,nrow = 2,align = "vh")
title <- ggdraw() + 
  draw_label(
    "",
    fontface = 'bold',size=15)
legend1 <- get_legend(p1)
plot_grid(
  title, pp,legend1,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 2)
)







# dataExp=dataEx$p
# dataExp[dataEx$p=="p=50"]=expression(italic(p)==50)
# dataExp[dataEx$p=="p=100"]=expression(italic(p)==100)
# dataExp[dataEx$p=="p=300"]=expression(italic(p)==300)
# dataEx$p=dataExp

# dataEx$Method=factor(dataEx$Method,levels = c("EFIR","GraphIR","MAVE","PrEFIR"))
# dataEx$p=factor(dataEx$Method,levels = c("p=50","p=100","p=300"))
# 
# dataEx%>%filter(Example=="Ex1")%>%
#   ggplot(mapping = aes(x=factor(n)))+
#   geom_boxplot(mapping=aes(y=Error,color=Method))+#
#   labs(x = "", y = "Root mean squared error")+  
#   facet_wrap(~(p),nrow = 1,scales="free")+
#   scale_fill_discrete(limits=c("EFIR","GraphIR","MAVE","PrEFIR"))+
#   scale_y_continuous(breaks = c(0.5,1.0,1.5,2.0),limits = c(0.1,2.1))+
#   theme_classic()+
#   # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
#   theme(legend.title=element_blank(),
#         legend.text=element_text(size=10),
#         plot.caption = element_text(hjust=0.5, size=15),
#         strip.text = element_blank(),
#         strip.background = element_rect(color="white",size=0),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"),
#         legend.position = "bottom"
#   )+ scale_color_jco()->p1





dataEx%>%filter(Example=="Ex1")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("EFIR","GraphIR","MAVE","PrEFIR"))), outlier.alpha = 0.1)+#
  labs(x = "", y = "Root mean squared error")+  
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  scale_fill_discrete(limits=c("EFIR","GraphIR","MAVE","PrEFIR"))+
  scale_y_continuous(breaks = c(0.5,1.0,1.5,2.0),limits = c(0.1,2.1))+
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=12),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_text(size=14,face="bold.italic"),
        strip.background = element_rect(color="white",size=0),
        # axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom"
  ) + scale_color_d3() ->p1
# + scale_color_jco()
# +  scale_color_manual(values = c("INDEP"="#55B7E6","NBSEL"="#193E8F","MAVE"="#F09739","FADR"="#E53528"))
#scale_y_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits = c(0.2,1.0))+
#+  scale_color_manual(values = c("IND"="red","NODE"="red","FADR-quasi"="blue","MAVE"="red","Oracle"="black"))+
#  scale_linetype_manual(values = c("IND"="solid","NODE"="dotted","FADR-quasi"="solid","MAVE"="dashed","Oracle"="solid"))->p1

p1
p1a<- p1 + theme(legend.position = "none")


dataEx%>%filter(Example=="Ex2")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("EFIR","GraphIR","MAVE","PrEFIR"))), outlier.alpha = 0.1)+#
  labs(x = expression(bolditalic(n)), y = "Area under the ROC curve")+  
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  scale_fill_discrete(limits=c("EFIR","GraphIR","MAVE","PrEFIR"))+
  scale_y_continuous(breaks = c(0.4,0.6,0.8,1.0),limits = c(0.4,1.0))+
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=12),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_blank(),
        strip.background = element_rect(color="white",size=0),
        # axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom"
  ) + scale_color_d3() ->p2

p2
p2a<- p2 + theme(legend.position = "none")



pp <- cowplot::plot_grid(p1a,p2a,nrow = 2,align = "vh")
title <- ggdraw() + 
  draw_label(
    "",
    fontface = 'bold',size=15)
legend1 <- get_legend(p1)
plot_grid(
  title, pp,legend1,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 2)
)

