# A joint simulation for alyzing SNR-CORR-Sparsity effects on EMVS and varbvs
library(EMVS)
library(varbvs)
library(bivas)
library(mvtnorm)
# set.seed(10)

trial <- 30

n <- 500
p <- 1000


#make auto-correlation
makcov2 <- function(p, rho){
  sigma <- matrix(0,p,p)
  sigma <- rho^abs(row(sigma)-col(sigma))
  return(sigma)
}


snr <- c(0.5,1,2,10)
corr <- 0#c(0,0.5,0.95)
alpha_true <- 0.05#c(0.05,0.1,0.2,0.4,0.8)

out <- data.frame(time=numeric(0),FDR=numeric(0),auc=numeric(0),power=numeric(0),error=numeric(0),method=character(0),corr=character(0),SNR=character(0),sparsity=character(0))
out
for(i in 1:length(snr)) {
  for(j in 1:length(corr)) {
    for(r in 1:length(alpha_true)){

      for(q in 1: trial) {
        cat(i,"/",length(snr),"snr; ",j,"/",length(corr),"corr; ",r,"/",length(alpha_true),"sparsity; ",q,"/",trial,"trial","\n")

        sb2_true <- 1
        gamma <- rbinom(p,1,alpha_true[r])
        beta_true <- rep(0,p)
        beta_true[gamma==1] <- rnorm(sum(gamma==1),0,sqrt(sb2_true))

        #corelated desien matrix
        # sigma <- makcov2(p,corr[j])
        # X <- rmvnorm(n,mean=rep(0,p),sigma=sigma)

        X <- replicate(p,rnorm(n))
        # X <- scale(X,center=TRUE,scale=FALSE)

        mu <- 0

        y0 <- mu+X%*%beta_true
        se2_true <- var(X%*%beta_true) / snr[i]
        y <- y0 + rnorm(n,0,sqrt(se2_true))
        # y <- y-mean(y)

        v0 = exp(seq(-10, -1, length.out = 20))
        v1 = 1
        beta_init = rep(1,p)
        sigma_init = var(y)
        a = b = 1
        epsilon = 10^{-5}

        time_varbvs <- system.time( out_varbvs <- varbvs(X=X,y=y,Z=NULL,verbose = F) )[3]
        # time_EMVS  <- system.time( out_EMVS  <- EMVS(Y = y, X = X, v0 = v0, v1 = v1, type = "betabinomial",
        #                                               independent = TRUE, beta_init = beta_init,
        #                                               sigma_init = sigma_init, direction = "backward",
        #                                               a = a, b = b, log_v0 = TRUE) )[3]
        time_EMVS1  <- system.time( out_EMVS1  <- EMVS(Y = y, X = X, v0 = v0, v1 = v1, type = "betabinomial",
                                                     independent = FALSE, beta_init = beta_init,
                                                     sigma_init = sigma_init, direction = "backward",
                                                     a = a, b = b, log_v0 = TRUE) )[3]
        time_EMVS2  <- system.time( out_EMVS2  <- EMVS(Y = y, X = X, v0 = v0, v1 = v1*10, type = "betabinomial",
                                                       independent = FALSE, beta_init = beta_init,
                                                       sigma_init = sigma_init, direction = "backward",
                                                       a = a, b = b, log_v0 = TRUE) )[3]
        time_EMVS3  <- system.time( out_EMVS3  <- EMVS(Y = y, X = X, v0 = v0, v1 = v1*100, type = "betabinomial",
                                                       independent = FALSE, beta_init = beta_init,
                                                       sigma_init = sigma_init, direction = "backward",
                                                       a = a, b = b, log_v0 = TRUE) )[3]
        time_EMVS4  <- system.time( out_EMVS4  <- EMVS(Y = y, X = X, v0 = v0, v1 = v1*1000, type = "betabinomial",
                                                       independent = FALSE, beta_init = beta_init,
                                                       sigma_init = sigma_init, direction = "backward",
                                                       a = a, b = b, log_v0 = TRUE) )[3]

        table_varbvs <- table((fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))<0.1),gamma)
        if(nrow(table_varbvs)==1) table_varbvs <- rbind(table_varbvs,0)

        # table_EMVS  <- table((fdr2FDR(1-out_EMVS$prob_inclusion[1,])<0.1),gamma)
        # if(nrow(table_EMVS)==1) table_EMVS <- rbind(table_EMVS,0)

        idx1 <- which.max(EMVSsummary(out_EMVS1)$log_g_function)
        table_EMVS1  <- table((fdr2FDR(1-out_EMVS1$prob_inclusion[idx1,])<0.1),gamma)
        if(nrow(table_EMVS1)==1) table_EMVS1 <- rbind(table_EMVS1,0)

        idx2 <- which.max(EMVSsummary(out_EMVS2)$log_g_function)
        table_EMVS2  <- table((fdr2FDR(1-out_EMVS2$prob_inclusion[idx2,])<0.1),gamma)
        if(nrow(table_EMVS2)==1) table_EMVS2 <- rbind(table_EMVS2,0)

        idx3 <- which.max(EMVSsummary(out_EMVS3)$log_g_function)
        table_EMVS3  <- table((fdr2FDR(1-out_EMVS3$prob_inclusion[idx3,])<0.1),gamma)
        if(nrow(table_EMVS3)==1) table_EMVS3 <- rbind(table_EMVS3,0)

        idx4 <- which.max(EMVSsummary(out_EMVS4)$log_g_function)
        table_EMVS4  <- table((fdr2FDR(1-out_EMVS4$prob_inclusion[idx4,])<0.1),gamma)
        if(nrow(table_EMVS4)==1) table_EMVS4 <- rbind(table_EMVS4,0)

        FDR_varbvs <- table_varbvs[2,1]/sum(table_varbvs[2,])
        # FDR_EMVS  <- table_EMVS[2,1]/sum(table_EMVS[2,])
        FDR_EMVS1  <- table_EMVS1[2,1]/sum(table_EMVS1[2,])
        FDR_EMVS2  <- table_EMVS2[2,1]/sum(table_EMVS2[2,])
        FDR_EMVS3  <- table_EMVS3[2,1]/sum(table_EMVS3[2,])
        FDR_EMVS4  <- table_EMVS4[2,1]/sum(table_EMVS4[2,])

        AUC_varbvs <- auc(gamma,as.vector(with(out_varbvs,alpha%*%normalizelogweights(logw))))
        # AUC_EMVS  <- auc(gamma,as.vector(out_EMVS$prob_inclusion[1,]))
        AUC_EMVS1  <- auc(gamma,as.vector(out_EMVS1$prob_inclusion[idx1,]))
        AUC_EMVS2  <- auc(gamma,as.vector(out_EMVS2$prob_inclusion[idx2,]))
        AUC_EMVS3  <- auc(gamma,as.vector(out_EMVS3$prob_inclusion[idx3,]))
        AUC_EMVS4  <- auc(gamma,as.vector(out_EMVS4$prob_inclusion[idx4,]))

        power_varbvs <- table_varbvs[2,2]/sum(table_varbvs[,2])
        # power_EMVS  <- table_EMVS[2,2]/sum(table_EMVS[,2])
        power_EMVS1  <- table_EMVS1[2,2]/sum(table_EMVS1[,2])
        power_EMVS2  <- table_EMVS2[2,2]/sum(table_EMVS2[,2])
        power_EMVS3  <- table_EMVS3[2,2]/sum(table_EMVS3[,2])
        power_EMVS4  <- table_EMVS4[2,2]/sum(table_EMVS4[,2])

        error_varbvs <- mean((beta_true-with(out_varbvs,(alpha*mu)%*%normalizelogweights(logw)))^2)
        # error_EMVS  <- mean((beta_true-out_EMVS$betas[1,])^2)
        error_EMVS1  <- mean((beta_true-out_EMVS1$betas[idx1,])^2)
        error_EMVS2  <- mean((beta_true-out_EMVS2$betas[idx2,])^2)
        error_EMVS3  <- mean((beta_true-out_EMVS3$betas[idx3,])^2)
        error_EMVS4  <- mean((beta_true-out_EMVS4$betas[idx4,])^2)

        # out <- rbind(out,data.frame(time=time_EMVS,FDR=FDR_EMVS,AUC=AUC_EMVS,power=power_EMVS,error=error_EMVS,method="EMVS",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=alpha_true[r]))
        out <- rbind(out,data.frame(time=time_EMVS1,FDR=FDR_EMVS1,AUC=AUC_EMVS1,power=power_EMVS1,error=error_EMVS1,method="EMVS1",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=alpha_true[r]))
        out <- rbind(out,data.frame(time=time_EMVS2,FDR=FDR_EMVS2,AUC=AUC_EMVS2,power=power_EMVS2,error=error_EMVS2,method="EMVS2",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=alpha_true[r]))
        out <- rbind(out,data.frame(time=time_EMVS3,FDR=FDR_EMVS3,AUC=AUC_EMVS3,power=power_EMVS3,error=error_EMVS3,method="EMVS3",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=alpha_true[r]))
        out <- rbind(out,data.frame(time=time_EMVS4,FDR=FDR_EMVS4,AUC=AUC_EMVS4,power=power_EMVS4,error=error_EMVS4,method="EMVS4",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=alpha_true[r]))
        out <- rbind(out,data.frame(time=time_varbvs,FDR=FDR_varbvs,AUC=AUC_varbvs,power=power_varbvs,error=error_varbvs,method="varbvs",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=alpha_true[r]))
      }
    }
  }
}

save(out,file="/home/share/mingxuan/bivas/simulation_results/EMVS_varbvs_weak.RData")

