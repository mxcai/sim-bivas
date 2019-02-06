# Analyzing the influence of fixed effects on BIVAS
library(pROC)
library(grpreg)
library(mvtnorm)
library(bivas)
library(varbvs)
# set.seed(10)

trial <- 30
useCore <- 10

n <- 1000
p <- 5000
K <- 250
l <- rep(p/K,K)
g <- rep(c(1:K),l)

rr <- 5

#make auto-correlation
makcov2 <- function(p, rho){
  sigma <- matrix(0,p,p)
  sigma <- rho^abs(row(sigma)-col(sigma))
  return(sigma)
}


snr <- 1#c(0.5,1,2)
corr <- 0#c(-0.5,-0.3,0,0.3,0.5)
pi_true <- c(0.05,0.1,0.2,0.4,0.8)
alpha_true <- c(0.8,0.4,0.2,0.1,0.05)
sparsity <- data.frame(pi=pi_true,alpha=alpha_true)

out <- data.frame(time=numeric(0),FDR=numeric(0),auc=numeric(0),groupAUC=numeric(0),power=numeric(0), error=numeric(0),
                  intercept=numeric(0),C1=numeric(0),C2=numeric(0),C3=numeric(0),C4=numeric(0),C5=numeric(0),
                  corr=character(0),SNR=character(0),sparsity=character(0))

get_gAUC <- function(eta0,coef,group){
  K <- max(group)
  gPos <- rep(0,K)
  for(i in 1:K){
    gPos[i] <- mean((coef[group==i])^2)
  }
  gAUC <- ifelse(max(gPos)==0,NA,auc(eta0,gPos/max(gPos)))
  return(gAUC)
}

for(i in 1:length(snr)) {
  for(j in 1:length(corr)) {
    for(r in 1:nrow(sparsity)){

      for(q in 1: trial) {
        cat(i,"/",length(snr),"snr; ",j,"/",length(corr),"corr; ",r,"/",nrow(sparsity),"sparsity; ",q,"/",trial,"trial","\n")

        sb2_true <- 0.1
        eta0 <- rbinom(K,1,sparsity$pi[r])
        eta <- rep(eta0,l)
        gamma <- rbinom(p,1,sparsity$alpha[r])
        beta_true <- rep(0,p)
        beta_true[(eta*gamma)==1] <- rnorm(sum(eta*gamma==1),0,sqrt(sb2_true))

        #corelated desien matrix
        X <- matrix(nrow=n,ncol=0)
        for (k in 1:K) {
          sigma <- makcov2(l[k],corr[j])
          X_k <- rmvnorm(n,mean=rep(0,l[k]),sigma=sigma)
          X <- cbind(X,X_k)
        }
        # X <- replicate(p,rnorm(n))
        # X <- scale(X,center=TRUE,scale=FALSE)

        Z <- matrix(rnorm(n*rr),n,rr)
        omega <- 1:rr
        mu <- Z %*% omega

        y0 <- mu+X%*%beta_true
        se2_true <- var(X%*%beta_true) / snr[i]
        y <- y0 + rnorm(n,0,sqrt(se2_true))
        # y <- y-mean(y)

        time_bivas  <- system.time( out_bivas  <- bivas(y,X,Z,group=g,coreNum = useCore,verbose = F) )[3]

        table_bivas  <- table(fdr(out_bivas,control="global")$FDR,gamma*eta)
        if(nrow(table_bivas)==1) table_bivas <- rbind(table_bivas,0)

        FDR_bivas  <- table_bivas[2,1]/sum(table_bivas[2,])

        AUC_bivas  <- auc(eta*gamma,as.vector(getPos(out_bivas)$var_pos))

        gAUC_bivas  <- pROC::auc(eta0,as.vector(getPos(out_bivas)$group_pos))

        power_bivas  <- sum(fdr(out_bivas,control="global")$FDR+(gamma*eta==1)==2) / sum(gamma*eta==1)

        power_bivas  <- table_bivas[2,2]/sum(table_bivas[,2])

        coef_bivas <- coef(out_bivas)
        error_bivas  <- mean((beta_true-coef_bivas$beta)^2)

        out <- rbind(out,data.frame(time=time_bivas,FDR=FDR_bivas,AUC=AUC_bivas,groupAUC=gAUC_bivas,power=power_bivas,error=error_bivas,intercept=coef_bivas$cov[1],C1=coef_bivas$cov[2],C2=coef_bivas$cov[3],C3=coef_bivas$cov[4],C4=coef_bivas$cov[5],C5=coef_bivas$cov[6],corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=paste(sparsity$pi[r],":",sparsity$alpha[r],sep = "")))
      }
    }
  }
}

save(out,file="/home/share/mingxuan/bivas/simulation_results/Bivas_FE_SNR1_rho0_pi.RData")
