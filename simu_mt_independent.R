# simulation of bivs-multi-task under violated conditions
# comparison of FDR, power and auc and error; signal/noise 0.3, 0.5 and 0.8;correlation 0.1, 0.3 and 0.5;
#
library(reshape2)
library(ggplot2)
library(pROC)
# library(grpreg)
library(mvtnorm)
library(bivas)
# library(IGESS)
library(glmnet)
library(varbvs)
set.seed(10)

trial <- 30
useCore <- 1

# #setting1: similar sample size
# nn  <- c(6,5,4) * 100

#setting2: large difference between sample size
nn  <- c(6,2,1) * 100

K <- 2000
l <- 3
p <- K * l

# fix the unrelated settings
# pi_true <- .05
# alpha_true <- 0.8
sb2_true <- 0.1

snr <- c(0.5,1,2,5,10)

# #total sparsity 0.04
pi_true <- c(0.05,0.1,0.2,0.4,0.8)
alpha_true <- c(0.8,0.4,0.2,0.1,0.05)

#total sparsity 0.004
# pi_true <- c(0.005,0.01,0.02,0.04,0.1,0.2,0.4,0.8)
# alpha_true <- c(0.8,0.4,0.2,0.1,0.04,0.02,0.01,0.005)

sparsity <- data.frame(pi=pi_true,alpha=alpha_true)

# out.bivas <- data.frame(total_Error=numeric(0),task1_Error=numeric(0),task2_Error=numeric(0),task3_Error=numeric(0),snr=character(0),sparsity=character(0))
# out.VB <- out.ridge <- out.lasso <- out.bivas

out.bivas <- matrix(,0,l+1+3)
out.VB <- out.ridge <- out.lasso <- out.bivas

for(r in 1:length(snr)){
  for(s in 1: length(pi_true)){
    for(q in 1: trial) {
      cat(r,"/",length(snr),"SNR; ",s,"/",length(pi_true),"Sparsity; ",q,"/",trial,"trial","\n")


      eta       <- rbinom(K,1,pi_true[s])
      gamma     <- matrix(rbinom(p,1,alpha_true[s]),K,l)
      beta_true <- matrix(rep(0,p),K,l)
      beta_true[(eta*gamma)==1] <- rnorm(sum(eta*gamma==1),0,sqrt(sb2_true))

      #corelated desien matrix
      X <- list()
      y <- list()

      for(t in 1:l){
        X[[t]] <- replicate(K,rnorm(nn[t]))
        mu       <- 5
        y0       <- mu+X[[t]]%*%beta_true[,t]
        se2_true <- var(y0) / snr[r]
        y[[t]]   <- y0 + rnorm(nn[t],0,sqrt(se2_true))
      }

      # Run ridge and lasso models.
      # set.seed(1337)
      ridge.fit <- mapply(cv.glmnet,X,y,alpha=0,standardize=F,SIMPLIFY = F)
      # set.seed(1337)
      lasso.fit <- mapply(cv.glmnet,X,y,alpha=1,standardize=F,SIMPLIFY = F)

      # Run individual level VB
      VB.fit <- mapply(varbvs,X=X,y=y,Z=vector("list",l),verbose = F,SIMPLIFY = F)

      # Run bivas analysis
      bivas_mt.fit <- bivas_mt(X=X,y=y,coreNum = useCore,verbose = F)


      #ridge error
      ridge.coef <- sapply(ridge.fit, function(x,s) drop(coef(x,s=s)), s="lambda.min")[-1,]
      e_t <- mean((ridge.coef - beta_true)^2)
      e_e <- colMeans((ridge.coef - beta_true)^2)
      out.ridge <- rbind(out.ridge,c(e_t,e_e,pi_true[s],alpha_true[s],snr[r]))

      #lasso error
      lasso.coef <- sapply(lasso.fit, function(x,s) drop(coef(x,s=s)), s="lambda.min")[-1,]
      e_t <- mean((lasso.coef - beta_true)^2)
      e_e <- colMeans((lasso.coef - beta_true)^2)
      out.lasso <- rbind(out.lasso,c(e_t,e_e,pi_true[s],alpha_true[s],snr[r]))

      #varbvs error
      VB.coef <- sapply(VB.fit,function(x) with(x,(alpha*mu)%*%normalizelogweights(logw)))
      e_t <- mean((VB.coef - beta_true)^2)
      e_e <- colMeans((VB.coef - beta_true)^2)
      out.VB <- rbind(out.VB,c(e_t,e_e,pi_true[s],alpha_true[s],snr[r]))

      #bivas error
      e_t <- mean((coef(bivas_mt.fit)$beta - beta_true)^2)
      e_e <- colMeans((coef(bivas_mt.fit)$beta - beta_true)^2)
      out.bivas <- rbind(out.bivas,c(e_t,e_e,pi_true[s],alpha_true[s],snr[r]))

    }
  }
}


# out <- data.frame(total_err=numeric(0),t1_err=numeric(0),t2_err=numeric(0),t3_err=numeric(0),pi=numeric(0),alpha=numeric(0),method=character(0))
out <- rbind(data.frame(out.bivas,method="BIVAS"),data.frame(out.VB,method="varbvs"),data.frame(out.ridge,method="Ridge"),data.frame(out.lasso,method="Lasso"))
names(out) <- c("Overall","Task1","Task2","Task3","pi","alpha","SNR","method")
out$Sparsity <- paste(out$pi,out$alpha,sep = ":")
out <- out[,!(names(out)%in%c("pi","alpha"))]
out1 <- melt(out,id=c("SNR","method","Sparsity"))

# save(out1,file="simu_mt_independent(sparsity0.04,sample1032).RData")


out1 <- out1[out1$SNR%in%c(0.5,1,2,10) & out1$Sparsity!="0.8:0.05" & out1$method!="Ridge" ,]
P <- ggplot(out1,aes(x=Task,y=error,color=method)) + geom_boxplot()  + facet_grid(SNR~Sparsity,labeller = label_both,scales = "free_y") + ylab("Mean Squared Error") + theme(legend.position="bottom")
