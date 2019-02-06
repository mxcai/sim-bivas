library(ggplot2)

#group BIVAS vs varbvs, cMCP, GEL
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/BF_Can.RData")
names(out) <- c("Time","FDR","AUC","groupAUC","Power","Error","Method","CORR","SNR","Sparsity")
levels(out$Method) <- c("BIVAS","varbvs","cMCP","GEL")
out <- out[out$CORR%in%c("-0.5","0","0.5"),]
x_grid <- c("SNR=0.5","SNR=1","SNR=2")
names(x_grid) <- levels(out$SNR)
y_grid <- c("CORR=-0.5","CORR=0","CORR=0.5")
names(y_grid) <- levels(droplevels(out$CORR))


P_BF_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "Supp_BF_AUC.pdf",width = 12,height = 9)
P_BF_AUC
dev.off()

P_BF_AUC   <- ggplot(subset(out,CORR==0),aes(x=Method,y=AUC,fill=Method)) + geom_boxplot() +
  facet_grid(SNR~Sparsity,labeller = labeller(SNR=x_grid))+
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  theme_grey(base_size = 15) + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "BF_AUC.pdf",width = 12,height = 9)
P_BF_AUC
dev.off()


P_BF_gAUC <- ggplot(out[out$Method!="varbvs",],aes(x=Sparsity,y=groupAUC,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylab("group-AUC") +
  scale_fill_manual(values=c("#CC6666", "#E69F00", "#56B4E9"),drop=TRUE) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "Supp_BF_gAUC.pdf",width = 12,height = 9)
P_BF_gAUC
dev.off()

P_BF_gAUC <- ggplot(subset(out,CORR==0 & Method!="varbvs"),aes(x=Method,y=groupAUC,fill=Method)) + geom_boxplot() +
  facet_grid(SNR~Sparsity,labeller = labeller(SNR=x_grid)) + ylab("group-AUC") +
  scale_fill_manual(values=c("#CC6666", "#E69F00", "#56B4E9"),drop=TRUE) +
  theme_grey(base_size = 15) + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "BF_gAUC.pdf",width = 12,height = 9)
P_BF_gAUC
dev.off()


P_BF_error <- ggplot(out,aes(x=Sparsity,y=Error,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,0.02)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "Supp_BF_error.pdf",width = 12,height = 9)
P_BF_error
dev.off()

P_BF_error <- ggplot(subset(out,CORR==0),aes(x=Method,y=Error,fill=Method)) + geom_boxplot() +
  facet_grid(SNR~Sparsity,labeller = labeller(SNR=x_grid)) + ylim(c(0,0.015)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "BF_error.pdf",width = 12,height = 9)
P_BF_error
dev.off()


P_BF_time <- ggplot(out,aes(x=Sparsity,y=Time,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,75)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  ylab("Time(second)") + theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "Supp_BF_time.pdf",width = 12,height = 9)
P_BF_time
dev.off()
P_BF_time <- ggplot(subset(out,CORR==0),aes(x=Method,y=Time,fill=Method)) + geom_boxplot() +
  facet_grid(SNR~Sparsity,labeller = labeller(SNR=x_grid)) + ylim(c(0,75)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  ylab("Time(second)") + theme_grey(base_size = 15) + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "BF_time.pdf",width = 12,height = 9)
P_BF_time
dev.off()

pdf(file = "B1_FDR.pdf",width = 12,height = 9)
P_B1_FDR <- ggplot(out[out$Method%in%c("BIVAS","varbvs"),],aes(x=Sparsity,y=FDR,fill=Method)) +
  geom_boxplot() + facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
P_B1_FDR
dev.off()

pdf(file = "B1_power.pdf",width = 12,height = 9)
P_B1_power <- ggplot(out[out$Method%in%c("BIVAS","varbvs"),],aes(x=Sparsity,y=Power,fill=Method)) +
  geom_boxplot() + facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
P_B1_power
dev.off()


#group BIVAS vs BSGS-SS
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/B2.RData")
names(out) <- c("Time","Error","Method","Sparsity","CORR")
out <- out[out$CORR==0.5,]

pdf(file="B2_time.pdf",width = 8,height = 8)
P_B2_time   <- ggplot(out,aes(x=Method,y=Time,fill=Method)) + geom_boxplot() +
  ylab(expression(Timing("log"[10]*"(Second)"))) + theme(legend.position="bottom") +  scale_y_log10() +
  facet_grid(~Sparsity,labeller = label_both,scale="free_y")
P_B2_time
dev.off()

pdf(file="B2_error.pdf",width = 8,height = 8)
P_B2_error <- ggplot(out,aes(x=Method,y=Error,fill=Method)) +
  geom_boxplot() + ylab("Mean Squared Error")  + theme(legend.position="bottom") +
  facet_grid(~Sparsity,labeller = label_both,scale="free_y")
P_B2_error
dev.off()



#multi-task simulation
load(file="/Users/cmx/Desktop/Bivas/R-package/bivas/simulation/simu_mt_independent(sparsity0.04,sample654).RData")
names(out1) <- c("SNR","Method","Sparsity","Task","Error")
out1 <- out1[out1$SNR%in%c(0.5,1,2) & out1$Method!="Ridge",]
pdf(file = "mt_error.pdf",width = 12,height = 6)
P_mt_error <- ggplot(out1,aes(x=Task,y=Error,fill=Method)) + geom_boxplot() + ylab("Mean Squared Error")  + theme(legend.position="bottom") + facet_grid(SNR~Sparsity,labeller = label_both,scale="free_y")
P_mt_error
dev.off()



#Manhattan plot of NFBC
load(file = "/Users/cmx/Desktop/Bivas/R-package/bivas/NFBC.RData")
info <- NFBC$info
rm(NFBC)
levels(info$V3) <- sub("([_].*)","",levels(info$V3))
levels(info$V3) <- sub("(chr)","",levels(info$V3))
levels(info$V3)[23] <- "23"
info$V3 <- as.numeric(levels(info$V3))[info$V3]
pheno <- c("CRP","Glucose","Insulin","TC","HDL","LDL","TG","BMI","SysBP","DiaBP")
pos_igess <- list()
pos_bivas <- list()


png(file="NFBC_HDL.png",res=400,height=1000,width=2500)
par(mfcol=c(1,3))

load(file = paste("/Users/cmx/Desktop/Bivas/R-package/NFBC/igess_NFBC_HDL.RData",sep=""))
load(file = paste("/Users/cmx/Desktop/Bivas/R-package/NFBC/bivas_NFBC_HDL_sp.RData",sep=""))

fdr <- fdr(fit_NFBC_HDL)
info[which(fdr$FDR==1),]    #show selected SNPs and genes

pos_igess <- igess_NFBC_HDL$gammas

pos_bivas <- getPos(fit_NFBC_HDL)

thr <- max(pos_igess,pos_bivas$var_pos,pos_bivas$group_pos)
thr <- ifelse((1-thr)<1e-8,1-1e-8,thr)

pos_igess <- ifelse(pos_igess>thr,thr,pos_igess)
pos_bivas$var_pos <- ifelse(pos_bivas$var_pos>thr,thr,pos_bivas$var_pos)
pos_bivas$group_pos <- ifelse(pos_bivas$group_pos>thr,thr,pos_bivas$group_pos)


manhattan(data.frame(pos=1-pos_igess,info),chr="V3",bp="V4",p="pos",snp="V8",
          genomewideline = -log10(0.1),suggestiveline = -log10(0.05),ylim=c(0,8.5),
          ylab=expression("-log"[10]*"(fdr)"),main="HDL: varbvs-SNP")

manhattan(data.frame(pos=1-pos_bivas$var_pos,info),chr="V3",bp="V4",p="pos",snp="V8",
          genomewideline = -log10(0.1),suggestiveline = -log10(0.05),ylim=c(0,8.5),
          ylab=expression("-log"[10]*"(fdr)"),main="HDL: BIVAS-SNP")

idx <- (with(info,match(levels(as.factor(V2)),V2)))
info_g <- info[idx,]

manhattan(data.frame(pos=1-pos_bivas$group_pos,info_g),chr="V3",bp="V4",p="pos",snp="V8",
          genomewideline = -log10(0.1),suggestiveline = -log10(0.05),ylim=c(0,8.5),
          ylab=expression("-log"[10]*"(fdr)"),main="HDL: BIVAS-gene")
dev.off()


#Manhattan plot of WTCCC
pheno <- c("RA","T1D")
pos_igess <- list()
pos_bivas <- list()

png(file="WTCCC_RA_T1D.png",res=400,height=2000,width=2500)
par(mfrow=c(2,3))
for(i in 1:length(pheno)){
  load(file = paste("/Users/cmx/Desktop/Bivas/R-package/WTCCC_result/igess_WTCCC_",pheno[i],".RData",sep=""))
  load(file = paste("/Users/cmx/Desktop/Bivas/R-package/WTCCC_result/bivas_WTCCC_",pheno[i],"_sp.RData",sep=""))
  info <- get(paste("info_",pheno[i],sep=""))

  pos_igess[[i]] <- get(paste("igess_WTCCC_",pheno[i],sep = ""))$gammas

  pos_bivas[[i]] <- getPos(get(paste("fit_WTCCC_",pheno[i],sep = "")))

  thr <- max(pos_bivas[[i]]$var_pos,pos_bivas[[i]]$group_pos,pos_igess[[i]])
  thr <- ifelse((1-thr)<1e-8,1-1e-8,thr)

  pos_igess[[i]] <- ifelse(pos_igess[[i]]>thr,thr,pos_igess[[i]])
  pos_bivas[[i]]$var_pos <- ifelse(pos_bivas[[i]]$var_pos>thr,thr,pos_bivas[[i]]$var_pos)
  pos_bivas[[i]]$group_pos <- ifelse(pos_bivas[[i]]$group_pos>thr,thr,pos_bivas[[i]]$group_pos)


  manhattan(data.frame(pos=1-pos_igess[[i]],info),chr="CHR",bp="BP",p="pos",snp="SNP",
            genomewideline = -log10(0.1),suggestiveline = -log10(0.05),ylim=c(0,8.5),
            ylab=expression("-log"[10]*"(fdr)"),main=paste(pheno[i],": varbvs-SNP",sep = ""))

  manhattan(data.frame(pos=1-pos_bivas[[i]]$var_pos,info),chr="CHR",bp="BP",p="pos",snp="SNP",
            genomewideline = -log10(0.1),suggestiveline = -log10(0.05),ylim=c(0,8.5),
            ylab=expression("-log"[10]*"(fdr)"),main=paste(pheno[i],": BIVAS-SNP",sep = ""))

  idx <- (with(info,match(levels(as.factor(group)),group)))
  info_g <- info[idx,]

  manhattan(data.frame(pos=1-pos_bivas[[i]]$group_pos,info_g),chr="CHR",bp="BP",p="pos",snp="SNP",
            genomewideline = -log10(0.1),suggestiveline = -log10(0.05),ylim=c(0,8.5),
            ylab=expression("-log"[10]*"(fdr)"),main=paste(pheno[i],": BIVAS-gene",sep = ""))

}
dev.off()








#Timing of multithread NFBC-HDL 33277.267 18416.335 11319.268  7077.545
# time <- data.frame(time=c(21628.874, 11673.835, 5828.43, 4847.200, 4189.476),thread=c(2,4,8,12,16))
# time <- data.frame(time=c(33277.267, 18416.335, 11319.268, 7077.545),thread=c(1,2,4,8))
time <- data.frame(time=c(33322.206, 18414.793, 11380.385, 8399.517, 7017.151),thread=c(1,2,4,6,8))
pdf(file="time_HDL.pdf",width = 5,height=5.5)
ggplot(time, aes(x=thread, y=time)) + geom_line() + geom_point() + xlab("Thread") + ylab("Time (Seconds)") + theme_grey(base_size = 15)
dev.off()


#Lq at EM each procedure in NFBC-HDL
maxL <- max(sapply(fit_NFBC_HDL$Lq_List,length))
Lq <- as.data.frame(lapply(fit_NFBC_HDL$Lq_List, FUN=function(x) c(x,rep(max(x),maxL-length(x)))))
names(Lq) <- fit_NFBC_HDL$logodds
Lq <- cbind(Lq,id=1:maxL)
Lq <- melt(Lq,id="id")
names(Lq) <- c("iteration","logodds","Lq")
pdf(file = "Lq_trace.pdf",width = 5,5,height=5.5)
ggplot(Lq,aes(x=iteration,y=Lq,colour=as.numeric(as.character(logodds)),group=logodds)) + geom_line() + theme_grey(base_size = 15) + theme(legend.position = c(0.85, 0.2),legend.background = element_rect(color = "black",fill = "grey90", size = 1, linetype = "solid"))+
  scale_colour_gradient(name = "logodds", low = "red", high = "yellow") + coord_cartesian(ylim = c(-163000, max(Lq$Lq))) + ylab("Lower Bound")
dev.off()


#Lq trace NFBC-HDL
pdf(file = "Lq_HDL.pdf",width = 5,5,height=5.5)
ggplot(data.frame(Lq=fit_NFBC_HDL$Lq,pi=fit_NFBC_HDL$logodds),aes(x=pi,y=Lq)) +# geom_line() +
  geom_point() + xlab(expression(paste("log"[10]*"odds(",pi,")"))) + ylab("Lower Bound") + theme_grey(base_size = 15)
dev.off()

#Alpha trace NFBC-HDL
pdf(file = "alpha_HDL.pdf",width = 5,5,height=5.5)
ggplot(data.frame(Lq=fit_NFBC_HDL$alpha,pi=fit_NFBC_HDL$logodds),aes(x=pi,y=Lq)) +# geom_line() +
  geom_point() + xlab(expression(paste("log"[10]*"odds(",pi,")"))) + ylab(expression(alpha)) + theme_grey(base_size = 15)
dev.off()











# Supp simulation
########################################################################### large rho ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/BF_largeRho.RData")
names(out) <- c("Time","FDR","AUC","groupAUC","Power","Error","Method","CORR","SNR","Sparsity")
levels(out$Method) <- c("BIVAS","varbvs","cMCP","GEL")
x_grid <- c("SNR=0.5","SNR=1","SNR=2")
names(x_grid) <- levels(out$SNR)
y_grid <- c("CORR=0.9","CORR=0.95","CORR=0.99")
names(y_grid) <- levels(droplevels(out$CORR))


P_rho_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BF_AUC_bigRho.pdf",width = 12,height = 9)
P_rho_AUC
dev.off()


P_rho_gAUC <- ggplot(out[out$Method!="varbvs",],aes(x=Sparsity,y=groupAUC,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylab("group-AUC") +
  scale_fill_manual(values=c("#CC6666", "#E69F00", "#56B4E9"),drop=TRUE) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BF_gAUC_bigRho.pdf",width = 12,height = 9)
P_rho_gAUC
dev.off()


P_rho_error <- ggplot(out,aes(x=Sparsity,y=Error,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,0.02)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BF_error_bigRho.pdf",width = 12,height = 9)
P_rho_error
dev.off()


P_rho_time <- ggplot(out,aes(x=Sparsity,y=Time,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,75)) +
  scale_fill_manual(values=c("#CC6666", "#66CC99", "#E69F00", "#56B4E9")) +
  ylab("Time(second)") + theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BF_time_bigRho.pdf",width = 12,height = 9)
P_rho_time
dev.off()


P_rho_FDR <- ggplot(out[out$Method%in%c("BIVAS","varbvs"),],aes(x=Sparsity,y=FDR,fill=Method)) +
  geom_boxplot() + facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_B1_FDR_bigRho.pdf",width = 12,height = 9)
P_rho_FDR
dev.off()


P_rho_power <- ggplot(out[out$Method%in%c("BIVAS","varbvs"),],aes(x=Sparsity,y=Power,fill=Method)) +
  geom_boxplot() + facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_B1_power_bigRho.pdf",width = 12,height = 9)
P_rho_power
dev.off()



########################################################################### Fixed Effects -- SNR ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/Bivas_FE_rho0_pi005_SNR.RData")
names(out) <- c("Time","FDR","AUC","groupAUC","Power","Error","Intercept","C1","C2","C3","C4","C5","CORR","SNR","Sparsity")
x_grid <- c("SNR=0.5","SNR=1","SNR=2")
names(x_grid) <- levels(out$SNR)
y_grid <- "CORR=0"#c("CORR=-0.5","CORR=0","CORR=0.5")
names(y_grid) <- levels(droplevels(out$CORR))


P_Bivas_AUC   <- ggplot(out,aes(x=SNR,y=AUC)) + geom_boxplot() +
  # facet_grid(.~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title="AUC")
# pdf(file = "Supp_AUC_FE.pdf",width = 12,height = 9)
# P_Bivas_AUC
# dev.off()


P_Bivas_gAUC <- ggplot(out,aes(x=SNR,y=groupAUC)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylab("group-AUC") +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title="group AUC") + scale_fill_discrete(drop=FALSE)
# pdf(file = "Supp_gAUC_FE.pdf",width = 12,height = 9)
# P_Bivas_gAUC
# dev.off()


P_Bivas_error <- ggplot(out,aes(x=SNR,y=Error)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,0.02)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title="Mean Squared Error")
# pdf(file = "Supp_error_FE.pdf",width = 12,height = 9)
# P_Bivas_error
# dev.off()


P_Bivas_FDR <- ggplot(out,aes(x=SNR,y=FDR)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title="FDR")
# pdf(file = "Supp_FDR_FE.pdf",width = 12,height = 9)
# P_Bivas_FDR
# dev.off()


P_Bivas_power <- ggplot(out,aes(x=SNR,y=Power)) +geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title="power")
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_power
# dev.off()


P_Bivas_C1 <- ggplot(out,aes(x=SNR,y=C1)) +geom_boxplot() + geom_hline(aes(yintercept=1),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title=expression(omega[1]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C1
# dev.off()

P_Bivas_C2 <- ggplot(out,aes(x=SNR,y=C2)) +geom_boxplot() + geom_hline(aes(yintercept=2),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title=expression(omega[2]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C2
# dev.off()

P_Bivas_C3 <- ggplot(out,aes(x=SNR,y=C3)) +geom_boxplot() + geom_hline(aes(yintercept=3),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title=expression(omega[3]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C3
# dev.off()

P_Bivas_C4 <- ggplot(out,aes(x=SNR,y=C4)) +geom_boxplot() + geom_hline(aes(yintercept=4),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title=expression(omega[4]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C4
# dev.off()

P_Bivas_C5 <- ggplot(out,aes(x=SNR,y=C5)) +geom_boxplot() + geom_hline(aes(yintercept=5),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,title=expression(omega[5]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C5
# dev.off()

Fig1 <- grid.arrange(P_Bivas_AUC,P_Bivas_gAUC,P_Bivas_error,P_Bivas_FDR,P_Bivas_power,
                     P_Bivas_C1,P_Bivas_C2,P_Bivas_C3,P_Bivas_C4,P_Bivas_C5,ncol=5,nrow=2)
ggsave("/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_FE_SNR.pdf",Fig1,width=30,height=14)



########################################################################### Fixed Effects -- rho ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/Bivas_FE_SNR1_pi005_rho.RData")
names(out) <- c("Time","FDR","AUC","groupAUC","Power","Error","Intercept","C1","C2","C3","C4","C5","CORR","SNR","Sparsity")
out <- out[out$CORR%in%c("-0.5","0","0.5"),]
x_grid <- "SNR=1"#c("SNR=0.5","SNR=1","SNR=2")
names(x_grid) <- levels(out$SNR)
y_grid <- c("CORR=-0.5","CORR=0","CORR=0.5")
names(y_grid) <- levels(droplevels(out$CORR))


P_Bivas_AUC   <- ggplot(out,aes(x=CORR,y=AUC)) + geom_boxplot() +
  # facet_grid(.~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title="AUC")
# pdf(file = "Supp_AUC_FE.pdf",width = 12,height = 9)
# P_Bivas_AUC
# dev.off()


P_Bivas_gAUC <- ggplot(out,aes(x=CORR,y=groupAUC)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylab("group-AUC") +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title="group AUC") + scale_fill_discrete(drop=FALSE)
# pdf(file = "Supp_gAUC_FE.pdf",width = 12,height = 9)
# P_Bivas_gAUC
# dev.off()


P_Bivas_error <- ggplot(out,aes(x=CORR,y=Error)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,0.02)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title="Mean Squared Error")
# pdf(file = "Supp_error_FE.pdf",width = 12,height = 9)
# P_Bivas_error
# dev.off()


P_Bivas_FDR <- ggplot(out,aes(x=CORR,y=FDR)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title="FDR")
# pdf(file = "Supp_FDR_FE.pdf",width = 12,height = 9)
# P_Bivas_FDR
# dev.off()


P_Bivas_power <- ggplot(out,aes(x=CORR,y=Power)) +geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title="power")
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_power
# dev.off()


P_Bivas_C1 <- ggplot(out,aes(x=CORR,y=C1)) +geom_boxplot() + geom_hline(aes(yintercept=1),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title=expression(omega[1]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C1
# dev.off()

P_Bivas_C2 <- ggplot(out,aes(x=CORR,y=C2)) +geom_boxplot() + geom_hline(aes(yintercept=2),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title=expression(omega[2]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C2
# dev.off()

P_Bivas_C3 <- ggplot(out,aes(x=CORR,y=C3)) +geom_boxplot() + geom_hline(aes(yintercept=3),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title=expression(omega[3]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C3
# dev.off()

P_Bivas_C4 <- ggplot(out,aes(x=CORR,y=C4)) +geom_boxplot() + geom_hline(aes(yintercept=4),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title=expression(omega[4]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C4
# dev.off()

P_Bivas_C5 <- ggplot(out,aes(x=CORR,y=C5)) +geom_boxplot() + geom_hline(aes(yintercept=5),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25)) +
  labs(y=NULL,x=expression(rho),title=expression(omega[5]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C5
# dev.off()

Fig1 <- grid.arrange(P_Bivas_AUC,P_Bivas_gAUC,P_Bivas_error,P_Bivas_FDR,P_Bivas_power,
                     P_Bivas_C1,P_Bivas_C2,P_Bivas_C3,P_Bivas_C4,P_Bivas_C5,ncol=5,nrow=2)
ggsave("/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_FE_rho.pdf",Fig1,width=30,height=14)



########################################################################### Fixed Effects -- sparsity ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/Bivas_FE_SNR1_rho0_Sparsity.RData")
names(out) <- c("Time","FDR","AUC","groupAUC","Power","Error","Intercept","C1","C2","C3","C4","C5","CORR","SNR","Sparsity")
# out <- out[out$CORR%in%c("-0.5","0","0.5"),]
x_grid <- "SNR=1"#c("SNR=0.5","SNR=1","SNR=2")
names(x_grid) <- levels(out$SNR)
y_grid <- "CORR=0"#c("CORR=-0.5","CORR=0","CORR=0.5")
names(y_grid) <- levels(droplevels(out$CORR))


P_Bivas_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC)) + geom_boxplot() +
  # facet_grid(.~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title="AUC")
# pdf(file = "Supp_AUC_FE.pdf",width = 12,height = 9)
# P_Bivas_AUC
# dev.off()


P_Bivas_gAUC <- ggplot(out,aes(x=Sparsity,y=groupAUC)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylab("group-AUC") +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title="group AUC") + scale_fill_discrete(drop=FALSE)
# pdf(file = "Supp_gAUC_FE.pdf",width = 12,height = 9)
# P_Bivas_gAUC
# dev.off()


P_Bivas_error <- ggplot(out,aes(x=Sparsity,y=Error)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,0.02)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title="Mean Squared Error")
# pdf(file = "Supp_error_FE.pdf",width = 12,height = 9)
# P_Bivas_error
# dev.off()


P_Bivas_FDR <- ggplot(out,aes(x=Sparsity,y=FDR)) + geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title="FDR")
# pdf(file = "Supp_FDR_FE.pdf",width = 12,height = 9)
# P_Bivas_FDR
# dev.off()


P_Bivas_power <- ggplot(out,aes(x=Sparsity,y=Power)) +geom_boxplot() +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title="power")
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_power
# dev.off()


P_Bivas_C1 <- ggplot(out,aes(x=Sparsity,y=C1)) +geom_boxplot() + geom_hline(aes(yintercept=1),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title=expression(omega[1]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C1
# dev.off()

P_Bivas_C2 <- ggplot(out,aes(x=Sparsity,y=C2)) +geom_boxplot() + geom_hline(aes(yintercept=2),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title=expression(omega[2]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C2
# dev.off()

P_Bivas_C3 <- ggplot(out,aes(x=Sparsity,y=C3)) +geom_boxplot() + geom_hline(aes(yintercept=3),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title=expression(omega[3]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C3
# dev.off()

P_Bivas_C4 <- ggplot(out,aes(x=Sparsity,y=C4)) +geom_boxplot() + geom_hline(aes(yintercept=4),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title=expression(omega[4]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C4
# dev.off()

P_Bivas_C5 <- ggplot(out,aes(x=Sparsity,y=C5)) +geom_boxplot() + geom_hline(aes(yintercept=5),color="red",linetype="dashed") +
  # facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),text = element_text(size=25),axis.text.x = element_text(size=15)) +
  labs(y=NULL,title=expression(omega[5]))
# pdf(file = "Supp_power_FE.pdf",width = 12,height = 9)
# P_Bivas_C5
# dev.off()

Fig1 <- grid.arrange(P_Bivas_AUC,P_Bivas_gAUC,P_Bivas_error,P_Bivas_FDR,P_Bivas_power,
                     P_Bivas_C1,P_Bivas_C2,P_Bivas_C3,P_Bivas_C4,P_Bivas_C5,ncol=5,nrow=2)
ggsave("/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_FE_Sparsity.pdf",Fig1,width=30,height=14)



########################################################################### group selection weak signal ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/Bivas_weak.RData")
fdr_bivas <- fdrControl(out_bivas,0.01,"global")
fdr_varbvs <- fdr2FDR(1-with(out_varbvs,alpha%*%normalizelogweights(logw)))
fdr_i <- fdr_bivas$FDR[as.logical(eta)]
fdr_g <- fdr_bivas$groupFDR[as.logical(eta0)]

out1 <- data.frame(FDR=c(fdr_i,fdr_g),Vid=c(paste("V",rep(1:20,12),sep=""),rep("Group",12)),Gid=paste("G",c(t(replicate(20,1:12)),1:12),sep=""))
out1$Vid <- factor(out1$Vid,levels = c(paste("V",1:20,sep=""),"Group"))
out1$Gid <- factor(out1$Gid,levels = paste("G",1:12,sep=""))
out1$star <- ""
out1$star[out1$FDR<0.05] <- "*"

P1 <- ggplot(out1, aes(x=Gid, y=Vid)) +
  geom_tile(aes(fill = -log10(FDR)), colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue") +
  geom_text(aes(label = star),size=10, vjust = 0.8) +
  labs(fill=expression(-log[10](FDR)),y="FDR",x="group") +
  theme(text = element_text(size=25), axis.text.y =  element_text(size=17), axis.text.x = element_text(size=15),legend.position = "right")

out2 <- data.frame(FDR=fdr_varbvs[as.logical(eta)],Vid=paste("V",rep(1:20,12),sep=""),Gid=paste("G",c(t(replicate(20,1:12))),sep=""))
out2$Vid <- factor(out2$Vid,levels = c(paste("V",1:20,sep="")))
out2$Gid <- factor(out2$Gid,levels = paste("G",1:12,sep=""))
out2$star <- ""
out2$star[out2$FDR<0.05] <- "*"

P2 <- ggplot(out2, aes(x=Gid, y=Vid)) +
  geom_tile(aes(fill = -log10(FDR)), colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue") +
  geom_text(aes(label = star),size=10, vjust = 0.8) +
  labs(fill=expression(-log[10](FDR)),y="FDR",x="group") +
  theme(text = element_text(size=25), axis.text.y =  element_text(size=17), axis.text.x = element_text(size=15),legend.position = "right")


pdf("/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_group_weak_bivas.pdf",width=15,height=10)
P1
dev.off()
pdf("/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_group_weak_varbvs.pdf",width=15,height=10)
P2
dev.off()


########################################################################### EMVS vs varbvs: SNR ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/EMVS_varbvs_weak.RData")
names(out) <- c("Time","FDR","AUC","Power","Error","Method","CORR","SNR","Sparsity")
# levels(out$Method) <- c("EMVS","varbvs")
x_grid <- c("SNR=0.5","SNR=1","SNR=2","SNR=10")
names(x_grid) <- levels(out$SNR)
# y_grid <- c("CORR=0.9","CORR=0.95","CORR=0.99")
# names(y_grid) <- levels(droplevels(out$CORR))

P_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC,fill=Method)) + geom_boxplot() +
  facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_EMVS_AUC.pdf",width = 12,height = 9)
P_AUC
dev.off()


P_error <- ggplot(out,aes(x=Sparsity,y=Error,fill=Method)) + geom_boxplot() +
  facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_EMVS_error.pdf",width = 12,height = 9)
P_error
dev.off()


# P_time <- ggplot(out,aes(x=Sparsity,y=Time,fill=Method)) + geom_boxplot() +
#   facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
#   ylab("Time(second)") + theme_grey(base_size = 15) + theme(legend.position="bottom")
# pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_EMVS_time.pdf",width = 12,height = 9)
# P_time
# dev.off()Carbonetto P, Stephens M. Scalable variational inference for Bayesian variable selection in regression, and its accuracy in genetic association studies[J]. Bayesian analysis, 2012, 7(1): 73-108.


P_FDR <- ggplot(out,aes(x=Sparsity,y=FDR,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_EMVS_FDR.pdf",width = 12,height = 9)
P_FDR
dev.off()


P_power <- ggplot(out,aes(x=Sparsity,y=Power,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_EMVS_power.pdf",width = 12,height = 9)
P_power
dev.off()

########################################################################### BIVAS vs BSGS ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/BIVAS_BSGS_beta.RData")
# out <- data.frame(method=rep(c("BSGS_simple","BSGS_sample","bivas"),each=length(coef_bivas$beta)),beta=c(PE_simple$beta.est,PE_sample$beta.est,coef_bivas$beta))
out <- data.frame(BSGS=PE_simple$beta.est,bivas=coef_bivas$beta)

P_beta <- ggplot(out,aes(x=BSGS,y=bivas))+
  geom_smooth(method=lm,  linetype="dashed",color="darkred") + geom_abline(slope = 1) + geom_point(size=2)
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/BIVAS_BSGS_beta.pdf",width = 8,height = 6)
P_beta
dev.off()

# load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/BIVAS_BSGS3.RData")
# names(out) <- c("Time","FDR","group-FDR","AUC","group-AUC","Power","group-Power","Error","Method","CORR","SNR","Sparsity")
# # levels(out$Method) <- c("EMVS","varbvs")
# x_grid <- "1"#c("SNR=0.5","SNR=1","SNR=2","SNR=10")
# names(x_grid) <- levels(out$SNR)
# # y_grid <- c("CORR=0.9","CORR=0.95","CORR=0.99")
# # names(y_grid) <- levels(droplevels(out$CORR))
#
# P_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC,fill=Method)) + geom_boxplot() +
#   facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
#   theme_grey(base_size = 15) + theme(legend.position="bottom")
# pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BSGS_AUC.pdf",width = 12,height = 9)
# P_AUC
# dev.off()
#
#
# P_error <- ggplot(out,aes(x=Sparsity,y=Error,fill=Method)) + geom_boxplot() +
#   facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
#   ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom")
# pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BSGS_error.pdf",width = 12,height = 9)
# P_error
# dev.off()
#
#
# P_time <- ggplot(out,aes(x=Sparsity,y=Time,fill=Method)) + geom_boxplot() +
#   facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
#   ylab("Time(second)") + theme_grey(base_size = 15) + theme(legend.position="bottom")
# pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BSGS_time.pdf",width = 12,height = 9)
# P_time
# dev.off()
#
#
# P_FDR <- ggplot(out,aes(x=Sparsity,y=FDR,fill=Method)) +
#   geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
#   geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
#   theme_grey(base_size = 15) + theme(legend.position="bottom")
# pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BSGS_FDR.pdf",width = 12,height = 9)
# P_FDR
# dev.off()
#
#
# P_power <- ggplot(out,aes(x=Sparsity,y=Power,fill=Method)) +
#   geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
#   theme_grey(base_size = 15) + theme(legend.position="bottom")
# pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_BSGS_power.pdf",width = 12,height = 9)
# P_power
# dev.off()


########################################################################### BIVAS initialization ###########################################################################
load("/Users/cmx/Desktop/Bivas/R-package/bivas_simulation/Bivas_init_test2.RData")
names(out) <- c("Time","FDR","AUC","group-AUC","Power","Error","pi","alpha","sb2","se2","Method","CORR","SNR","Sparsity")
# levels(out$Method) <- c("EMVS","varbvs")
x_grid <- "1"#c("SNR=0.5","SNR=1","SNR=2","SNR=10")
names(x_grid) <- levels(out$SNR)
# y_grid <- c("CORR=0.9","CORR=0.95","CORR=0.99")
# names(y_grid) <- levels(droplevels(out$CORR))

P_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC,fill=Method)) + geom_boxplot() +
  facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_init_AUC.pdf",width = 12,height = 9)
P_AUC
dev.off()


P_error <- ggplot(out,aes(x=Sparsity,y=Error,fill=Method)) + geom_boxplot() +
  facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_init_error.pdf",width = 12,height = 9)
P_error
dev.off()


P_FDR <- ggplot(out,aes(x=Sparsity,y=FDR,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  geom_hline(aes(yintercept=0.1),color="red",linetype="dashed") +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_init_FDR.pdf",width = 12,height = 9)
P_FDR
dev.off()


P_power <- ggplot(out,aes(x=Sparsity,y=Power,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_init_power.pdf",width = 12,height = 9)
P_power
dev.off()


P_pi <- ggplot(out,aes(x=Sparsity,y=pi,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) + geom_hline(yintercept = 0.8,linetype="dotted", color = "red")+
  theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_init_pi.pdf",width = 12,height = 9)
P_pi
dev.off()


P_alpha <- ggplot(out,aes(x=Sparsity,y=alpha,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) + geom_hline(yintercept = 0.05,linetype="dotted", color = "red")+
theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_initS_alpha.pdf",width = 12,height = 9)
P_alpha
dev.off()


P_sb2 <- ggplot(out,aes(x=Sparsity,y=sb2,fill=Method)) +
  geom_boxplot() + facet_grid(.~SNR,labeller = labeller(SNR=x_grid)) + geom_hline(yintercept = 0.1,linetype="dotted", color = "red")+
theme_grey(base_size = 15) + theme(legend.position="bottom")
pdf(file = "/Users/cmx/Desktop/Bivas/Bivas-paper/image/Supp_init_sb2.pdf",width = 12,height = 9)
P_sb2
dev.off()

