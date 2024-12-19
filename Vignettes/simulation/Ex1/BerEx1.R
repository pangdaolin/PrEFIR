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
    d=2; k=1
    
    
    
    rmse_GFM <- numeric()
    
    time1 <- numeric()
    time2 <- numeric()
    time3 <- numeric()
    
    alpha <- rep(0,p)#rnorm(p, mean = 0, sd = 1)rep(0,p)
    Phi <- rbind(matrix(c(-1,1,-1,1,-1,1,1,1,-1,1),1,p),matrix(c(0,0,0,0,1,0,0,1,1,1),1,p))
    Gamma <- matrix(c(1,1,1,1,0,-1,1,0,0,0),k,p)
    
    system.time(for (runs in 1:100) {
      set.seed(927017+runs)
      f <- runif(n*d,-2,2)
      ff <- matrix(f,n,d)
      Y <- rnorm(n, mean = 0, sd = 1)
      V = Y#as.numeric(Y>0)
      stru1 <- matrix(alpha, n, p, byrow = TRUE) + ff %*% (Phi) + V %*% (Gamma)
      q <- exp(stru1)/(exp(stru1)+1)
      
      X <- matrix(0,n,p)
      for (i in 1:n) {
        while (sum(X[i,])==0) {
          for (j in 1:p) {
            X[i,j] <- rbinom(1,1,prob = q[i, j])
          }
        }
      }
      
      X.train <- X[1:(n/2),];V.train <- V[1:(n/2)];ff.train <- ff[1:(n/2),]
      X.test <- X[(n/2+1):n,];V.test <- V[(n/2+1):n];ff.test <- ff[(n/2+1):n,]
      Y.train <- V.train
      Y.test<- V.test
      
      #############################
      X.all_list <- list(rbind(X.train,X.test))
      gfm2 <- gfm(X.all_list, 'binomial', q=3, verbose = FALSE, maxIter=5)
      F.train <- gfm2[["hH"]][1:(n/2),]
      F.test <- gfm2[["hH"]][(n/2+1):n,]
      
      train.object <- data.frame(X=F.train,Y=Y.train)
      test.object <- data.frame(X=F.test)
      GAM1 <- try(gam(Y~s(X.1)+s(X.2),data=train.object),silent = T)
      pre1 <- stats::predict(GAM1,test.object)
      rmse_GFM <- c(rmse_GFM,sqrt(sum((pre1-Y.test)^2)/length(Y.test)))
      
      cat("Ber n =",n/2,"p =",p,"Iteration step",runs,"completed","\n")
      
    })
    save.image(paste0("~/pdl/review_p2/biometrics/simu/Berd2/GFMBerEx1n",n/2,"p",p,".RData"))
  }
}


##########################
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggsci)

for (ni in c(100,200,600)) {
  for (pj in c(50,100,300)) {
    for (Exk in 1:2) {
      load(paste0("~/pdl/review_p2/biometrics/simu/Berd2/Ber","Ex",Exk,"n",ni/2,"p",pj,".RData"))
      load(paste0("~/pdl/review_p2/biometrics/simu/Berd2/GFMBer","Ex",Exk,"n",ni/2,"p",pj,".RData"))
      
      if(Exk==1&pj==50&ni==100){
        dataEx <- data.frame(n=rep(n/2,5*100),p=rep(paste0("p=",p),5*100),
                             Method=c(rep("EFIR",100),rep("GraphIR",100),rep("PrEFIR",100),rep("MAVE",100),rep("GFM",100)),
                             Example=rep(if(length(levels(factor(V)))==2){"Ex2"}else{"Ex1"},5*100),
                             Error=c((rmse3),(rmse5),(rmse2),(rmse7),rmse_GFM))
      } else{
        Ex <- if(length(levels(factor(V)))==2){"Ex2"}else{"Ex1"}
        if(Ex=="Ex2"){
          dataEx <- rbind(dataEx,data.frame(n=rep(n/2,5*100),p=rep(paste0("p=",p),5*100),
                                            Method=c(rep("EFIR",100),rep("GraphIR",100),rep("PrEFIR",100),rep("MAVE",100),rep("GFM",100)),
                                            Example=rep(Ex,5*100),
                                            Error=c((mAUC3),(mAUC5),(mAUC2),(mAUC7),mAUC_GFM)))
        }else{
          dataEx <- rbind(dataEx,data.frame(n=rep(n/2,5*100),p=rep(paste0("p=",p),5*100),
                                            Method=c(rep("EFIR",100),rep("GraphIR",100),rep("PrEFIR",100),rep("MAVE",100),rep("GFM",100)),
                                            Example=rep(Ex,5*100),
                                            Error=c((rmse3),(rmse5),(rmse2),(rmse7),rmse_GFM)))
        } 
      }
    }
  }
}


dataEx%>%filter(Example=="Ex1")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("GFM","EFIR","GraphIR","MAVE","PrEFIR"))), outlier.alpha = 0.1)+#
  labs(x = "", y = "Root mean squared error")+  
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  scale_fill_discrete(limits=c("EFIR","GraphIR","MAVE","PrEFIR"))+
  scale_y_continuous(breaks = c(0.5,1.0,1.5,2.0),limits = c(0.1,2.1))+
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
# + scale_color_jco()
# +  scale_color_manual(values = c("INDEP"="#55B7E6","NBSEL"="#193E8F","MAVE"="#F09739","FADR"="#E53528"))
#scale_y_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits = c(0.2,1.0))+
#+  scale_color_manual(values = c("IND"="red","NODE"="red","FADR-quasi"="blue","MAVE"="red","Oracle"="black"))+
#  scale_linetype_manual(values = c("IND"="solid","NODE"="dotted","FADR-quasi"="solid","MAVE"="dashed","Oracle"="solid"))->p1

p1
p1a<- p1 + theme(legend.position = "none")


dataEx%>%filter(Example=="Ex2")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("GFM","EFIR","GraphIR","MAVE","PrEFIR"))), outlier.alpha = 0.1)+#
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
legend1 <- ggpubr::get_legend(p1)
plot_grid(
  title, pp,legend1,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 2)
)







dataEx%>%filter(Example=="Ex1")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("INDEP","NBSEL","MAVE","FADR"))))+#
  labs(x = "", y = "Root mean squared error")+  
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  scale_fill_discrete(limits=c("INDEP","NBSEL","MAVE","FADR"))+
  scale_y_continuous(breaks = c(0.5,1.0,1.5,2.0),limits = c(0.2,2.0))+
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=10),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_blank(),
        strip.background = element_rect(color="white",size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom"
  )+ scale_color_jco()->p1
#scale_y_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits = c(0.2,1.0))+
#+  scale_color_manual(values = c("IND"="red","NODE"="red","FADR-quasi"="blue","MAVE"="red","Oracle"="black"))+
#  scale_linetype_manual(values = c("IND"="solid","NODE"="dotted","FADR-quasi"="solid","MAVE"="dashed","Oracle"="solid"))->p1

p1
p1a<- p1 + theme(legend.position = "none")


dataEx%>%filter(Example=="Ex2")%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_boxplot(mapping=aes(y=Error,color=factor(Method,levels = c("INDEP","NBSEL","MAVE","FADR"))))+#
  facet_wrap(~factor(p,levels=c("p=50","p=100","p=300")),nrow = 1,scales="free")+
  labs(x = "n", y = "Area under the ROC curve")+  
  scale_fill_discrete(limits=c("INDEP","NBSEL","MAVE","FADR"))+
  scale_y_continuous(breaks = c(0.4,0.6,0.8,1.0),limits = c(0.3,1.0))+
  theme_classic()+
  # theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=10),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_blank(),
        strip.background = element_rect(color="white",size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "bottom"
  )+ scale_color_jco()->p2
#scale_y_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits = c(0.2,1.0))+
#+  scale_color_manual(values = c("IND"="red","NODE"="red","FADR-quasi"="blue","MAVE"="red","Oracle"="black"))+
#  scale_linetype_manual(values = c("IND"="solid","NODE"="dotted","FADR-quasi"="solid","MAVE"="dashed","Oracle"="solid"))->p1

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





mean(rmse1)
mean(rmse1_0)
mean(rmse1_1)
mean(rmse2)
mean(rmse2_0)
mean(rmse2_1)
mean(rmse3)
mean(rmse3_1)
mean(rmse4)
mean(rmse4_0)
mean(rmse4_1)
mean(rmse5)
mean(rmse5_1)
mean(rmse6)


median(rmse1)
median(rmse1_0)
median(rmse2)
median(rmse2_0)
median(rmse3)
median(rmse4)
median(rmse4_0)
median(rmse5)
median(rmse7)

sd(rmse1)
sd(rmse1_0)
sd(rmse2)
sd(rmse2_0)
sd(rmse3)
sd(rmse4)
sd(rmse4_0)
sd(rmse10)

library(ggplot2)
library(cowplot)
library(dplyr)




rmsedata <- data.frame(rmse=c(rmse1,rmse2,rmse3,rmse4,rmse5),
                       methods=factor(c(rep("FADR-IRLS",runs),rep("FADR-quasi",runs),rep("EFDR",runs),rep("Oracle",runs),rep("EFDR-GM",runs)),
                                      level = c("EFDR","EFDR-GM","FADR-IRLS","FADR-quasi","Oracle")))

p1 = ggplot(rmsedata, aes(x=methods, y=rmse, color = methods)) +
  geom_boxplot()+labs(title="", x="", y="RMSE")+ylim(0.20,1.5)
p1


dataEx3 <- data.frame(n=rep(n/2,5),p=rep(paste0("p=",p),5),
                      Method=c("EFDR","EFDR-GM","FADR-IRLS","FADR-quasi","Oracle"),
                      ltpye=c(1,2,1,2,1),
                      Error=c(median(rmse3),median(rmse10),median(rmse1),median(rmse2),median(rmse4)))

dataEx3 <- rbind(dataEx1,data.frame(n=rep(n/2,5),p=rep(paste0("p=",p),5),
                                    Method=c("EFDR","EFDR-GM","FADR-IRLS","FADR-quasi","Oracle"),
                                    ltpye=c(1,2,1,2,1),
                                    Error=c(median(rmse3),median(rmse10),median(rmse1),median(rmse2),median(rmse4))))
library(ggplot2)
library(cowplot)
library(dplyr)

dataEx3%>%
  ggplot(mapping = aes(x=factor(n)))+
  geom_line(mapping=aes(y=Error,group=factor(Method),color=factor(Method),linetype=factor(Method)),size=0.8)+
  facet_wrap(~p,scales="free")+
  xlab("Sample Size")+
  ylab("RMSE in Example 3")+
  scale_y_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7),limits = c(0.2,0.8))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=10),
        plot.caption = element_text(hjust=0.5, size=15),
        strip.text = element_text(size=14,face="bold.italic",color="white"),
        strip.background = element_rect(color="white",size=0),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "bottom"
  )+
  scale_color_manual(values = c("EFDR"="red","EFDR-GM"="red","FADR-IRLS"="blue","FADR-quasi"="blue","Oracle"="black"))+
  scale_linetype_manual(values = c("EFDR"="solid","EFDR-GM"="dotted","FADR-IRLS"="solid","FADR-quasi"="dotted","Oracle"="solid"))->p3

p3

pp <- cowplot::plot_grid(p1,p2,p3,nrow = 3)
title <- ggdraw() + 
  draw_label(
    "Poisson model",
    fontface = 'bold',size=15)
plot_grid(
  title, pp,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 2)
)
