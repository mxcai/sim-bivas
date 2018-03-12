# method to be compared: bivas, BSGS-SS, BGL-SS
# comparison of FDR, power and auc and error; sparsity 0.05:0.8, 0.1:0.4, 0.2:0.2, 0.4:0.1, 0.8:0.05
library(pROC)
# library(grpreg)
library(mvtnorm)
library(bivas)
library(IGESS)
library(MBSGS)
# set.seed(10)

trial <- 30
useCore <- 10

n <- 200
p <- 1000
K <- 100
l <- rep(p/K,K)
g <- rep(c(1:K),l)

#make auto-correlation
makcov2 <- function(p, rho){
  sigma <- matrix(0,p,p)
  sigma <- rho^abs(row(sigma)-col(sigma))
  return(sigma)
}


# fix the unrelated settings
snr <- 1
corr <- c(-0.5,-0.3,0,0.3,0.5)

pi_true <- c(0.05,0.1,0.2,0.4,0.8)
alpha_true <- c(0.8,0.4,0.2,0.1,0.05)
sparsity <- data.frame(pi=pi_true,alpha=alpha_true)


out <- data.frame(time=numeric(0),error=numeric(0),method=character(0),sparsity=character(0),corr=character(0))

for(i in 1:nrow(sparsity)) {
  for(t in 1:length(corr)){
    for(q in 1: trial) {
      cat(i,"/",nrow(sparsity),"sparsity; ",t,"/",length(corr),"sparsity; ",q,"/",trial,"trial","\n")

      sb2_true <- 0.1
      eta0 <- rbinom(K,1,sparsity$pi[i])
      eta <- rep(eta0,l)
      gamma <- rbinom(p,1,sparsity$alpha[i])
      beta_true <- rep(0,p)
      beta_true[(eta*gamma)==1] <- rnorm(sum(eta*gamma==1),0,sqrt(sb2_true))

      #corelated desien matrix
      X <- matrix(nrow=n,ncol=0)
      for (k in 1:K) {
        sigma <- makcov2(l[k],corr[t])
        X_k <- rmvnorm(n,mean=rep(0,l[k]),sigma=sigma)
        X <- cbind(X,X_k)
      }
      # X <- replicate(p,rnorm(n))
      # X <- scale(X,center=TRUE,scale=FALSE)
      mu <- 0#rnorm(1)

      y0 <- mu+X%*%beta_true
      se2_true <- var(y0) / snr
      y <- y0 + rnorm(n,0,sqrt(se2_true))
      # y <- y-mean(y)

      # out_varbvs <- varbvs(X=X,y=y,Z=NULL,verbose = F)
      # out_igess  <- iGess(X,y,opts = list(max_iter=600, display_gap = 100))
      time_bivas <- system.time(out_bivas  <- bivas(y,X,group=g,coreNum = useCore,verbose = F))[3]
      time_BSGS <- system.time(out_BSGS   <- BSGSSS(y,X,group_size=l,niter=500,burnin=100,num_update = 20,niter.update = 20))[3]
      # time_BGL <- system.time(out_BGL    <- BGLSS(y,X,group_size=l,niter=500,burnin=100,num_update = 20,niter.update = 20))[3]

      # FDR_varbvs <- 1 - sum((fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))<0.1)+(gamma*eta==1)==2) / sum(fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))<0.1)
      # FDR_bivas  <- 1 - sum(fdr(out_bivas,control="global")$FDR+(gamma*eta==1)==2) / sum(fdr(out_bivas,control="global")$FDR)
      # FDR_igess  <- 1 - sum((fdr2FDR(1-out_igess$gammas)<0.1)+(gamma*eta==1)==2) / sum(fdr2FDR(1-out_igess$gammas)<0.1)
      # FDR_BSGS  <- 1 - sum((out_BSGS$pos_median!=0) + (gamma*eta==1)==2) / sum((out_BSGS$pos_median!=0))
      # FDR_BGL  <- 1 - sum((out_BGL$pos_median!=0) + (gamma*eta==1)==2) / sum((out_BGL$pos_median!=0))

      # auc_varbvs <- auc(eta*gamma,as.vector(with(out_varbvs,alpha%*%normalizelogweights(logw))))
      # auc_bivas  <- auc(eta*gamma,as.vector(getPos(out_bivas)$var_pos))
      # auc_igess  <- auc(eta*gamma,as.vector(out_igess$gammas))
      # auc_BSGS  <-

      # power_varbvs <- sum((fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))<0.1)+(gamma*eta==1)==2) / sum(gamma*eta==1)
      # power_bivas  <- sum(fdr(out_bivas,control="global")$FDR+(gamma*eta==1)==2) / sum(gamma*eta==1)
      # power_igess  <- sum((fdr2FDR(1-out_igess$gammas)<0.1)+(gamma*eta==1)==2) / sum(gamma*eta==1)
      # power_BSGS  <- sum((out_BSGS$pos_median!=0)+(gamma*eta==1)==2) / sum(gamma*eta==1)
      # power_BGL  <- sum((out_BGL$pos_median!=0)+(gamma*eta==1)==2) / sum(gamma*eta==1)

      # error_varbvs <- mean((beta_true-with(out_varbvs,(alpha*mu)%*%normalizelogweights(logw)))^2)
      error_bivas  <- mean((beta_true-coef(out_bivas)$beta)^2)
      # error_igess  <- mean((beta_true-out_igess$gammas*out_igess$mu)^2)
      error_BSGS   <- mean((beta_true-out_BSGS$pos_median)^2)
      # error_BGL    <- mean((beta_true-out_BGL$pos_median)^2)

      out <- rbind(out,data.frame(time=time_bivas,error=error_bivas,method="BIVAS",sparsity=paste(sparsity$pi[i],":",sparsity$alpha[i],sep = ""),corr=paste(corr[t])))
      out <- rbind(out,data.frame(time=time_BSGS,error=error_BSGS,method="BSGS-SS",sparsity=paste(sparsity$pi[i],":",sparsity$alpha[i],sep = ""),corr=paste(corr[t])))
      # out <- rbind(out,data.frame(FDR=FDR_BGL,time=time_BGL,power=power_BGL,error=error_BGL,method="BGL-SS",sparsity=paste(sparsity$pi[i],":",sparsity$alpha[i],sep = "")))
    }
  }
}

save(out,file="B2.RData")
