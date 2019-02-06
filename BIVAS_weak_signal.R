# Analyzing the performances of BIVAS and varbvs under weak signal -- demonstrate the ability of identifying group when all variables are unidentified
library(pROC)
library(grpreg)
library(mvtnorm)
library(bivas)
library(varbvs)
# set.seed(10)

trial <- 1
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
pi_true <- 0.05#c(0.05,0.1,0.2,0.4,0.8)
alpha_true <- 0.8#c(0.8,0.4,0.2,0.1,0.05)
sparsity <- data.frame(pi=pi_true,alpha=alpha_true)

out <- data.frame(time=numeric(0),FDR=numeric(0),auc=numeric(0),groupAUC=numeric(0),power=numeric(0),
                  error=numeric(0),method=character(0),corr=character(0),SNR=character(0),sparsity=character(0))

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
        X <- replicate(p,rnorm(n))

        mu <- 0

        y0 <- mu+X%*%beta_true
        se2_true <- var(X%*%beta_true) / snr[i]
        y <- y0 + rnorm(n,0,sqrt(se2_true))

        time_varbvs <- system.time( out_varbvs <- varbvs(X=X,y=y,Z=NULL,verbose = F) )[3]
        time_bivas  <- system.time( out_bivas  <- bivas(y,X,group=g,coreNum = useCore,verbose = F) )[3]

        table_varbvs <- table((fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))<0.1),gamma*eta)
        if(nrow(table_varbvs)==1) table_varbvs <- rbind(table_varbvs,0)
        table_bivas  <- table(fdr(out_bivas,control="global")$FDR,gamma*eta)
        if(nrow(table_bivas)==1) table_bivas <- rbind(table_bivas,0)

        FDR_varbvs <- table_varbvs[2,1]/sum(table_varbvs[2,])
        FDR_bivas  <- table_bivas[2,1]/sum(table_bivas[2,])

        AUC_varbvs <- auc(eta*gamma,as.vector(with(out_varbvs,alpha%*%normalizelogweights(logw))))
        AUC_bivas  <- auc(eta*gamma,as.vector(getPos(out_bivas)$var_pos))

        gAUC_varbvs <- get_gAUC(eta0,with(out_varbvs,(alpha*mu)%*%normalizelogweights(logw)),g)
        gAUC_bivas  <- pROC::auc(eta0,as.vector(getPos(out_bivas)$group_pos))

        power_varbvs <- sum((fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))<0.1)+(gamma*eta==1)==2) / sum(gamma*eta==1)
        power_bivas  <- sum(fdr(out_bivas,control="global")$FDR+(gamma*eta==1)==2) / sum(gamma*eta==1)

        power_varbvs <- table_varbvs[2,2]/sum(table_varbvs[,2])
        power_bivas  <- table_bivas[2,2]/sum(table_bivas[,2])

        error_varbvs <- mean((beta_true-with(out_varbvs,(alpha*mu)%*%normalizelogweights(logw)))^2)
        coef_bivas <- coef(out_bivas)
        error_bivas  <- mean((beta_true-coef_bivas$beta)^2)

        out <- rbind(out,data.frame(time=time_bivas,FDR=FDR_bivas,AUC=AUC_bivas,groupAUC=gAUC_bivas,power=power_bivas,error=error_bivas,method="BIVAS",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=paste(sparsity$pi[r],":",sparsity$alpha[r],sep = "")))
        out <- rbind(out,data.frame(time=time_varbvs,FDR=FDR_varbvs,AUC=AUC_varbvs,groupAUC=gAUC_varbvs,power=power_varbvs,error=error_varbvs,method="varbvs",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=paste(sparsity$pi[r],":",sparsity$alpha[r],sep = "")))
      }
    }
  }
}

save(out_bivas,out_varbvs,eta,eta0,out,file="/home/share/mingxuan/bivas/simulation_results/Bivas_weak.RData")

