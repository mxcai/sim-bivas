# method to be compared: bivas, BSGS-SS, BGL-SS
# comparison of FDR, power and auc and error; sparsity 0.05:0.8, 0.1:0.4, 0.2:0.2, 0.4:0.1, 0.8:0.05
library(pROC)
library(mvtnorm)
library(bivas)
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

      time_bivas <- system.time(out_bivas  <- bivas(y,X,group=g,coreNum = useCore,verbose = F))[3]
      time_BSGS <- system.time(out_BSGS   <- BSGSSS(y,X,group_size=l,niter=500,burnin=100,num_update = 20,niter.update = 20))[3]


      error_bivas  <- mean((beta_true-coef(out_bivas)$beta)^2)
      error_BSGS   <- mean((beta_true-out_BSGS$pos_median)^2)

      out <- rbind(out,data.frame(time=time_bivas,error=error_bivas,method="BIVAS",sparsity=paste(sparsity$pi[i],":",sparsity$alpha[i],sep = ""),corr=paste(corr[t])))
      out <- rbind(out,data.frame(time=time_BSGS,error=error_BSGS,method="BSGS-SS",sparsity=paste(sparsity$pi[i],":",sparsity$alpha[i],sep = ""),corr=paste(corr[t])))
    }
  }
}

save(out,file="B2.RData")
