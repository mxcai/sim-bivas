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

pdf(file = "BF_AUC.pdf",width = 12,height = 9)
P_BF_AUC   <- ggplot(out,aes(x=Sparsity,y=AUC,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) +
  theme_grey(base_size = 15) + theme(legend.position="bottom")
P_BF_AUC
dev.off()

pdf(file = "BF_gAUC.pdf",width = 12,height = 9)
P_BF_gAUC <- ggplot(out[out$Method!="varbvs",],aes(x=Sparsity,y=groupAUC,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylab("group-AUC") +
  theme_grey(base_size = 15) + theme(legend.position="bottom") + scale_fill_discrete(drop=FALSE)
P_BF_gAUC
dev.off()

pdf(file = "BF_error.pdf",width = 12,height = 9)
P_BF_error <- ggplot(out,aes(x=Sparsity,y=Error,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,0.02)) +
  ylab("Mean Squared Error") + theme_grey(base_size = 15) + theme(legend.position="bottom")
P_BF_error
dev.off()

pdf(file = "BF_time.pdf",width = 12,height = 9)
P_BF_time <- ggplot(out,aes(x=Sparsity,y=Time,fill=Method)) + geom_boxplot() +
  facet_grid(CORR~SNR,labeller = labeller(SNR=x_grid,CORR=y_grid)) + ylim(c(0,75)) +
  ylab("Time(second)") + theme_grey(base_size = 15) + theme(legend.position="bottom")
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
  geom_point() + xlab(expression(paste("logodds(",pi,")"))) + ylab("Lower Bound") + theme_grey(base_size = 15)
dev.off()

#Alpha trace NFBC-HDL
pdf(file = "alpha_HDL.pdf",width = 5,5,height=5.5)
ggplot(data.frame(Lq=fit_NFBC_HDL$alpha,pi=fit_NFBC_HDL$logodds),aes(x=pi,y=Lq)) +# geom_line() +
  geom_point() + xlab(expression(paste("logodds(",pi,")"))) + ylab(expression(alpha)) + theme_grey(base_size = 15)
dev.off()


