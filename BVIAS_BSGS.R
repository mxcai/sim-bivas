require(BSGS)
library(bivas)
set.seed(2)
Num.Of.Iteration <- 5000
Num.of.Iter.Inside.CompWise <- 100
num.of.obs <- 50
num.of.covariates <- 100
num.of.group.var <- 10
Group.Index <- rep(1:10, each = 10)

pair.corr <- 0.
Sigma <- matrix(pair.corr, num.of.covariates, num.of.covariates)
diag(Sigma) <- rep(1,num.of.covariates)
X <- mvrnorm(n = num.of.obs, rep(0, num.of.covariates), Sigma)
beta.true <- rep(0, num.of.covariates)
beta.true[c(7, 8, 9, 11, 12, 43, 77)] <- c(3.2, 3.2, 3.2, 1.5, 1.5, -1.5, -2)
beta.true <- cbind(beta.true)
r.true <- (beta.true != 0) * 1
sigma2.true <-1
Y <- rnorm(num.of.obs, X %*% beta.true, sigma2.true)

## hyperparameters
nu <- 0
lambda <- 1

tau2.value <- rep(5, num.of.covariates)
rho.value <- rep(0.5, num.of.covariates)
theta.value <- rep(0.5, num.of.group.var)

## Initial values and stopping point
r.value <- rbinom(num.of.covariates, 1, 0.5)
eta.value <- rbinom(num.of.group.var, 1, 0.5)
beta.value <- cbind(c(t(solve(t(X) %*% X +
                                diag(1/5, num.of.covariates)) %*% t(X) %*% Y) )) # beta.true
sigma2.value <- 1
MCSE.Sigma2.Given <- 0.5


## Apply two sampling approaches to generate samples
outputSimple <- BSGS.Simple(Y, X, Group.Index, r.value, eta.value, beta.value,
                            tau2.value, rho.value, theta.value, sigma2.value, nu, lambda,
                            Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given)
outputSample <- BSGS.Sample(Y, X, Group.Index, r.value, eta.value, beta.value,
                            tau2.value, rho.value, theta.value, sigma2.value, nu, lambda,
                            Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given)
outputBIVAS <- bivas(Y,X,group=Group.Index,coreNum = 2,verbose = F, alpha=0.5, sb2 = 5,se2 = 1)

PE_simple <- BSGS.PE(outputSimple)
PE_sample <- BSGS.PE(outputSample)
coef_bivas <- coef(outputBIVAS)

save(PE_simple,PE_sample,coef_bivas,outputSimple,outputSample,outputBIVAS,file="/home/share/mingxuan/bivas/simulation_results/BIVAS_BSGS_beta.RData")






# # method to be compared: bivas, BSGS
# # comparison of FDR, power and auc and error; sparsity 0.1:0.4, 0.2:0.2, 0.4:0.1
# library(pROC)
# # library(grpreg)
# library(mvtnorm)
# library(bivas)
# library(IGESS)
# library(BSGS)
# # set.seed(1)
#
# trial <- 1
# useCore <- 10
#
# n <- 50
# p <- 100
# K <- 10
# l <- rep(p/K,K)
# g <- rep(c(1:K),l)
#
# #make auto-correlation
# makcov2 <- function(p, rho){
#   sigma <- matrix(0,p,p)
#   sigma <- rho^abs(row(sigma)-col(sigma))
#   return(sigma)
# }
#
#
# # fix the unrelated settings
# snr <- 10#c(0.5,1,2)
# corr <- 0#c(-0.5,-0.3,0,0.3,0.5)
#
# pi_true <- 0.4#c(0.2,0.4,0.5)
# alpha_true <- 0.25#c(0.5,0.25,0.2)
# sparsity <- data.frame(pi=pi_true,alpha=alpha_true)
#
# Num.Of.Iteration <- 2000
# Num.of.Iter.Inside.CompWise <- 100
#
# out <- data.frame(time=numeric(0),FDR=numeric(0),gFDR=numeric(0),auc=numeric(0),groupAUC=numeric(0),power=numeric(0),gpower=numeric(0),
#                   error=numeric(0),method=character(0),corr=character(0),SNR=character(0),sparsity=character(0))
#
# for(i in 1:length(snr)) {
#   for(j in 1:length(corr)) {
#     for(r in 1:nrow(sparsity)){
#
#       for(q in 1: trial) {
#         cat(i,"/",length(snr),"snr; ",j,"/",length(corr),"corr; ",r,"/",nrow(sparsity),"sparsity; ",q,"/",trial,"trial","\n")
#
#         sb2_true <- 1
#         # eta0 <- rbinom(K,1,sparsity$pi[i])
#         eta0 <- rep(0,K)
#         idx_eta0 <- sample(1:K,sparsity$pi[i]*K)
#         eta0[idx_eta0] <- 1
#         eta <- rep(eta0,l)
#         # gamma <- rbinom(p,1,sparsity$alpha[i])
#         gamma <- rep(0,p)
#         idx_gamma <- sample(1:p,sparsity$alpha[i]*p)
#         gamma[idx_gamma] <- 1
#         beta_true <- rep(0,p)
#         beta_true[(eta*gamma)==1] <- rnorm(sum(eta*gamma==1),0,sqrt(sb2_true))
#
#         #corelated desien matrix
#         # X <- matrix(nrow=n,ncol=0)
#         # for (k in 1:K) {
#         #   sigma <- makcov2(l[k],corr[t])
#         #   X_k <- rmvnorm(n,mean=rep(0,l[k]),sigma=sigma)
#         #   X <- cbind(X,X_k)
#         # }
#         X <- replicate(p,rnorm(n))
#         # X <- scale(X,center=TRUE,scale=FALSE)
#         mu <- 0#rnorm(1)
#
#         y0 <- mu+X%*%beta_true
#         se2_true <- var(y0) / snr
#         y <- y0 + rnorm(n,0,sqrt(se2_true))
#
#
#         ## hyperparameters
#         nu <- 0
#         lambda <- 1
#
#         tau2.value <- rep(1, p)
#         rho.value <- rep(0.5, p)
#         theta.value <- rep(0.5, K)
#
#         ## Initial values and stopping point
#         r.value <- rbinom(p, 1, 0.5)
#         eta.value <- rbinom(K, 1, 0.5)
#         beta.value <- cbind(c(t(solve(t(X) %*% X + diag(1/5, p)) %*% t(X) %*% y) )) # beta.true
#         sigma2.value <- 1
#         MCSE.Sigma2.Given <- 0.9
#
#
#         time_BSGS1 <- system.time(out_BSGS1   <- BSGS.Simple(y, X, g, r.value, eta.value, beta.value,
#                                                              tau2.value, rho.value, theta.value, sigma2.value, nu, lambda,
#                                                              Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given))[3]
#         cat("BSGS-simple finished \n")
#         time_BSGS2 <- system.time(out_BSGS2   <- BSGS.Sample(y, X, g, r.value, eta.value, beta.value,
#                                                              tau2.value, rho.value, theta.value, sigma2.value, nu, lambda,
#                                                              Num.of.Iter.Inside.CompWise, Num.Of.Iteration, MCSE.Sigma2.Given))[3]
#         cat("BSGS-sample finished \n")
#         time_bivas <- system.time(out_bivas  <- bivas(y,X,group=g,coreNum = useCore,verbose = F, alpha = 1/log(p+1),se2=var(y)/2,sb2=var(y)/2))[3]
#         cat("BIVAS finished \n")
#
#         PE_BSGS1 <- BSGS.PE(out_BSGS1)
#         PE_BSGS2 <- BSGS.PE(out_BSGS2)
#
#         table_BSGS1 <- table((fdr2FDR(1-PE_BSGS1$r.est)<0.1),gamma*eta)
#         if(nrow(table_BSGS1)==1) table_BSGS1 <- rbind(table_BSGS1,0)
#
#         table_BSGS2 <- table((fdr2FDR(1-PE_BSGS2$r.est)<0.1),gamma*eta)
#         if(nrow(table_BSGS2)==1) table_BSGS2 <- rbind(table_BSGS2,0)
#
#         table_bivas  <- table(fdr(out_bivas,control="global")$FDR,gamma*eta)
#         if(nrow(table_bivas)==1) table_bivas <- rbind(table_bivas,0)
#
#         gtable_BSGS1 <- table((fdr2FDR(1-PE_BSGS1$eta.est)<0.1),eta0)
#         if(nrow(gtable_BSGS1)==1) gtable_BSGS1 <- rbind(gtable_BSGS1,0)
#
#         gtable_BSGS2 <- table((fdr2FDR(1-PE_BSGS2$eta.est)<0.1),eta0)
#         if(nrow(gtable_BSGS2)==1) gtable_BSGS2 <- rbind(gtable_BSGS2,0)
#
#         gtable_bivas  <- table(fdr(out_bivas,control="global")$groupFDR,eta0)
#         if(nrow(gtable_bivas)==1) gtable_bivas <- rbind(gtable_bivas,0)
#
#         FDR_BSGS1 <- table_BSGS1[2,1]/sum(table_BSGS1[2,])
#         FDR_BSGS2 <- table_BSGS2[2,1]/sum(table_BSGS2[2,])
#         FDR_bivas  <- table_bivas[2,1]/sum(table_bivas[2,])
#
#         gFDR_BSGS1 <- gtable_BSGS1[2,1]/sum(gtable_BSGS1[2,])
#         gFDR_BSGS2 <- gtable_BSGS2[2,1]/sum(gtable_BSGS2[2,])
#         gFDR_bivas  <- gtable_bivas[2,1]/sum(gtable_bivas[2,])
#
#         AUC_BSGS1 <- auc(eta*gamma,as.vector(PE_BSGS1$r.est))
#         AUC_BSGS2 <- auc(eta*gamma,as.vector(PE_BSGS2$r.est))
#         AUC_bivas  <- auc(eta*gamma,as.vector(getPos(out_bivas)$var_pos))
#
#         gAUC_BSGS1 <- pROC::auc(eta0,as.vector(PE_BSGS1$eta.est))
#         gAUC_BSGS2 <- pROC::auc(eta0,as.vector(PE_BSGS2$eta.est))
#         gAUC_bivas  <- pROC::auc(eta0,as.vector(getPos(out_bivas)$group_pos))
#
#         power_BSGS1 <- table_BSGS1[2,2]/sum(table_BSGS1[,2])
#         power_BSGS2 <- table_BSGS2[2,2]/sum(table_BSGS2[,2])
#         power_bivas  <- table_bivas[2,2]/sum(table_bivas[,2])
#
#         gpower_BSGS1 <- gtable_BSGS1[2,2]/sum(gtable_BSGS1[,2])
#         gpower_BSGS2 <- gtable_BSGS2[2,2]/sum(gtable_BSGS2[,2])
#         gpower_bivas  <- gtable_bivas[2,2]/sum(gtable_bivas[,2])
#
#         error_BSGS1 <- mean((beta_true-PE_BSGS1$beta.est)^2)
#         error_BSGS2 <- mean((beta_true-PE_BSGS2$beta.est)^2)
#         error_bivas  <- mean((beta_true-coef(out_bivas)$beta)^2)
#
#         out <- rbind(out,data.frame(time=time_BSGS1,FDR=FDR_BSGS1,gFDR=gFDR_BSGS1,AUC=AUC_BSGS1,groupAUC=gAUC_BSGS1,power=power_BSGS1,gpower=gpower_BSGS1,error=error_BSGS1,method="BSGS1",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=paste(sparsity$pi[r],":",sparsity$alpha[r],sep = "")))
#         out <- rbind(out,data.frame(time=time_BSGS2,FDR=FDR_BSGS2,gFDR=gFDR_BSGS2,AUC=AUC_BSGS2,groupAUC=gAUC_BSGS2,power=power_BSGS2,gpower=gpower_BSGS2,error=error_BSGS2,method="BSGS2",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=paste(sparsity$pi[r],":",sparsity$alpha[r],sep = "")))
#         out <- rbind(out,data.frame(time=time_bivas,FDR=FDR_bivas,gFDR=gFDR_bivas,AUC=AUC_bivas,groupAUC=gAUC_bivas,power=power_bivas,gpower=gpower_bivas,error=error_bivas,method="BIVAS",corr=paste(corr[j]),SNR=paste(snr[i]),sparsity=paste(sparsity$pi[r],":",sparsity$alpha[r],sep = "")))
#       }
#     }
#   }
# }
# save(out,file="/home/share/mingxuan/bivas/simulation_results/BIVAS_BSGS2.RData")
