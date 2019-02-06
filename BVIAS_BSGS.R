# Compare BVIAS posterior mean estimates with that of BSGS
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
