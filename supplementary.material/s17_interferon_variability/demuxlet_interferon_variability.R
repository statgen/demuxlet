## load in the cd4
library(cellrangerRkit);
load("10x.singlets.RData");
library(pheatmap);

## load in the dc data
immvar.dc.stim <- read.table("IFNb.6.5.repeatability.txt",header=T, sep="\t")
immvar.cd4.stim <- read.table("cd4.ifn.txt",header=T, sep="\t",fill=T);

immvar.dc.ctrl <- read.table("unstim.0.repeatability.txt",header=T,sep="\t")
immvar.cd4.ctrl <- read.table("cd4.base.txt",header=T, sep="\t",fill=T);

load("cd4.stim.RData");
cd4.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd4.stim.cor.mean) <- fData(total)$symbol[match(names(cd4.stim.cor.mean),rownames(fData(total)))]

load("cd4.ctrl.RData");
cd4.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd4.ctrl.cor.mean) <- fData(total)$symbol[match(names(cd4.ctrl.cor.mean),rownames(fData(total)))]


load("cd14.stim.RData");
cd14.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd14.stim.cor.mean) <- fData(total)$symbol[match(names(cd14.stim.cor.mean),rownames(fData(total)))]

load("cd14.ctrl.RData");
cd14.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd14.ctrl.cor.mean) <- fData(total)$symbol[match(names(cd14.ctrl.cor.mean),rownames(fData(total)))]


load("dc.stim.RData");
dc.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(dc.stim.cor.mean) <- fData(total)$symbol[match(names(dc.stim.cor.mean),rownames(fData(total)))]

load("dc.ctrl.RData");
dc.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(dc.ctrl.cor.mean) <- fData(total)$symbol[match(names(dc.ctrl.cor.mean),rownames(fData(total)))]

load("cd16.stim.RData");
cd16.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd16.stim.cor.mean) <- fData(total)$symbol[match(names(cd16.stim.cor.mean),rownames(fData(total)))]

load("cd16.ctrl.RData");
cd16.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd16.ctrl.cor.mean) <- fData(total)$symbol[match(names(cd16.ctrl.cor.mean),rownames(fData(total)))]
##matched <- unique(union(union(immvar.dc.stim$X, cd4.stim$X), names(cd4.stim.cor.mean)));

load("cd19.stim.RData");
cd19.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd19.stim.cor.mean) <- fData(total)$symbol[match(names(cd19.stim.cor.mean),rownames(fData(total)))]

load("cd19.ctrl.RData");
cd19.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd19.ctrl.cor.mean) <- fData(total)$symbol[match(names(cd19.ctrl.cor.mean),rownames(fData(total)))]
##matched <- unique(union(union(immvar.dc.stim$X, cd4.stim$X), names(cd4.stim.cor.mean)));

load("cd8.stim.RData");
cd8.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd8.stim.cor.mean) <- fData(total)$symbol[match(names(cd8.stim.cor.mean),rownames(fData(total)))]

load("cd8.ctrl.RData");
cd8.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd8.ctrl.cor.mean) <- fData(total)$symbol[match(names(cd8.ctrl.cor.mean),rownames(fData(total)))]
##matched <- unique(union(union(immvar.dc.stim$X, cd4.stim$X), names(cd4.stim.cor.mean)));

load("cd56.stim.RData");
cd56.stim.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd56.stim.cor.mean) <- fData(total)$symbol[match(names(cd56.stim.cor.mean),rownames(fData(total)))]

load("cd56.ctrl.RData");
cd56.ctrl.cor.mean <- apply(mean.cor,1,mean,na.rm=T);
names(cd56.ctrl.cor.mean) <- fData(total)$symbol[match(names(cd56.ctrl.cor.mean),rownames(fData(total)))]
##matched <- unique(union(union(immvar.dc.stim$X, cd4.stim$X), names(cd4.stim.cor.mean)));

immvar.dc.stim$X <- as.character(immvar.dc.stim$X);
immvar.dc.stim$X[which(immvar.dc.stim$X=="C17orf48")] <- "ADPRM";
immvar.dc.stim$X[which(immvar.dc.stim$X=="C1orf150")] <- "GCSAML";
immvar.dc.stim$X[which(immvar.dc.stim$X=="IL1F9")] <- "IL36G";
immvar.dc.stim$X[which(immvar.dc.stim$X=="IL28A")] <- "IFNL2";
immvar.dc.stim$X[which(immvar.dc.stim$X=="IL28B")] <- "IFNL3";
immvar.dc.stim$X[which(immvar.dc.stim$X=="IL29")] <- "IFNL1";
immvar.dc.stim$X[which(immvar.dc.stim$X=="TMEM178")] <- "TMEM178A";
immvar.dc.stim$X[which(immvar.dc.stim$X=="C3orf59")] <- "MB21D2";
immvar.dc.stim$X[which(immvar.dc.stim$X=="C9orf30")] <- "MSANTD3";
immvar.dc.stim$X[which(immvar.dc.stim$X=="C6orf192")] <- "SLC18B1";
immvar.dc.stim$X[which(immvar.dc.stim$X=="NCRNA00152")] <- "LINC00152";
immvar.dc.stim$X[which(immvar.dc.stim$X=="MGC57346")] <- "CRHR1";


matched <- intersect(as.character(immvar.dc.stim$X), names(cd4.stim.cor.mean));
#matched <- intersect(as.character(immvar.cd4.ctrl$X), intersect(as.character(immvar.dc.stim$X), names(cd4.stim.cor.mean)));

df.out <- data.frame(immvar.modc.ctrl.re=immvar.dc.ctrl[match(matched, immvar.dc.ctrl$X),2],
immvar.modc.stim.re=immvar.dc.stim[match(matched, immvar.dc.stim$X),2],
##immvar.dc.ctrl.fe=immvar.dc.ctrl[match(matched, immvar.dc.ctrl$X),4],
##immvar.dc.stim.fe=immvar.dc.stim[match(matched, immvar.dc.stim$X),4],
##immvar.dc.ctrl.herit=immvar.dc.ctrl[match(matched, immvar.dc.ctrl$X),6],
##immvar.cd4.ctrl.re=immvar.cd4.ctrl[match(matched,immvar.cd4.ctrl$X),4], ##immvar.cd4.stim.re=immvar.cd4.stim[match(matched,immvar.cd4.stim$X),4],
cd14.ctrl=cd14.ctrl.cor.mean[match(matched, names(cd14.ctrl.cor.mean))]^2,
cd14.stim=cd14.stim.cor.mean[match(matched, names(cd14.stim.cor.mean))]^2,
dc.ctrl=dc.ctrl.cor.mean[match(matched, names(dc.ctrl.cor.mean))]^2,
dc.stim=dc.stim.cor.mean[match(matched, names(dc.stim.cor.mean))]^2,
cd16.ctrl=cd16.ctrl.cor.mean[match(matched, names(cd16.ctrl.cor.mean))]^2,
cd16.stim=cd16.stim.cor.mean[match(matched, names(cd16.stim.cor.mean))]^2,
cd4.ctrl=cd4.ctrl.cor.mean[match(matched, names(cd4.ctrl.cor.mean))]^2,
cd4.stim=cd4.stim.cor.mean[match(matched, names(cd4.stim.cor.mean))]^2,
cd8.ctrl=cd8.ctrl.cor.mean[match(matched, names(cd8.ctrl.cor.mean))]^2,
cd8.stim=cd8.stim.cor.mean[match(matched, names(cd8.stim.cor.mean))]^2,
cd56.ctrl=cd56.ctrl.cor.mean[match(matched, names(cd56.ctrl.cor.mean))]^2,
cd56.stim=cd56.stim.cor.mean[match(matched, names(cd56.stim.cor.mean))]^2
)

df.cor <- cor(df.out, use='complete.obs')
diag(df.cor) <- NA;
pheatmap(df.cor, cluster_cols=F, cluster_rows=F)
pheatmap(df.cor[c(1,2),-c(1,2)], cluster_cols=F, cluster_rows=F)
