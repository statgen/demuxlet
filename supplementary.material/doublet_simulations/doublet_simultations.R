library(data.table);
library(ggplot2);

header <- c("HG00096", "HG00101", "HG00104", "HG00106", "HG00109", "HG00116", "HG00126", "HG00130", "HG00139", "HG00141", "HG00146",
"HG00173", "HG00183", "HG00185", "HG00238", "HG00242", "HG00245", "HG00246", "HG00255", "HG00260", "HG00263", "HG00266", "HG00272",
"HG00275", "HG00280", "HG00281", "HG00311", "HG0089", "HG00334", "HG00338", "HG00341", "HG00343", "HG00350", "HG00361", "HG00362",
"HG00369", "HG00378", "HG00383", "HG00384", "HG01334", "NA06986", "NA07347", "NA11831", "NA11893", "NA11894", "NA11920", "NA11931",
"NA12043", "NA12348", "NA12763", "NA12814", "NA12842", "NA20502", "NA20513", "NA20534", "NA20588", "NA20759", "NA20790", "NA20799",
"NA20801", "NA20807", "NA20809", "NA20812", "NA20826");

probs <- NULL;
snps <- NULL;
inds <- NULL;
reads <- NULL;
errs <- NULL;

for (i in seq(5000,100000,5000)) {
  for (j in c(2,4,8,16,32,64)) {
    if(file.exists(paste(j,"inds",i,"snps.best",sep="."))) {
      cat("reading: ", paste(j,"inds",i,"snps.best",sep="."),"\n");

      rst <- fread(paste(j,"inds",i,"snps.best",sep="."),sep="\t");
      rst <- rst[grep("DBL",rst$BEST),];
      truth <- fread(paste(j,".inds.",i,".doublets",sep=""));
      ##browser();

      ##matched1.1 <- match(rst$DBL.1ST, header)==truth$V2[match(rst$BARCODE,truth$V1)];
      ##matched2.1 <- match(rst$DBL.2ND, header)==truth$V3[match(rst$BARCODE,truth$V1)];

      ##matched1.2 <- match(rst$DBL.1ST, header)==truth$V3[match(rst$BARCODE,truth$V1)];
      ##matched2.2 <- match(rst$DBL.2ND, header)==truth$V2[match(rst$BARCODE,truth$V1)];

      ##hits <- length(c(intersect(which(matched1.1),which(matched2.1)),
      ##intersect(which(matched1.2),which(matched2.2))));

      hits <- length(na.omit(match(rst$BARCODE, truth$new.b)));

      probs <- c(probs, hits/308);
      errs <- c(errs, (nrow(rst)-hits)/nrow(rst))
      snps <- c(snps, mean(rst$N.SNP));
      reads <- c(reads, mean(rst$RD.TOTL)/(1-0.95))
      inds <- c(inds, j);
    }
  }
}

theoretical <- data.frame(theory=(1-1/c(2,4,8,16,32,64)),inds=c(2,4,8,16,32,64));

df <- data.frame(probs=probs, errs=errs, snps=snps, reads=reads, inds=inds);
a <- ggplot(df, aes(snps, probs, color=as.factor(inds)))+geom_point()+geom_line()+theme_bw()+geom_hline(aes(yintercept=theory, color=as.factor(inds)),linetype="dotted", data=theoretical);
ggsave(a, file="snps.vs.prob.pdf");

df <- data.frame(probs=probs, errs=errs, snps=snps, reads=reads, inds=inds);
a <- ggplot(df, aes(snps, errs, color=as.factor(inds)))+geom_point()+geom_line()+theme_bw()+scale_y_continuous(limits = c(0, 1));
ggsave(a, file="snps.vs.errs.pdf");

b <- ggplot(df, aes(reads, probs, color=as.factor(inds)))+geom_point()+geom_line()+theme_bw();
ggsave(b, file="reads.vs.prob.pdf");

#
# ##pdf("snps.vs.probs.pdf"); plot(snps, probs); dev.off();
#
# ## let's make some bins of SNPs to see the probability of correctly identifying
#
# ##matched <- c(trio.matched, snps.75000.matched, snps.50000.matched, snps.25000.matched, snps.10000.matched, snps.5000.matched)
# ##snps <- c(trio.rst$N.SNP, snps.75000.rst$N.SNP, snps.50000.rst$N.SNP, snps.25000.rst$N.SNP, snps.10000.rst$N.SNP, snps.5000.rst$N.SNP);
#
# ##output <- NULL;
#
# n.snps <- seq(1,100,1);
#
# for(n.snp in n.snps) {
#     indices <- which(trio.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
#     df <- rbind(df, data.frame(probs=length(which(trio.matched[indices]))/length(indices), snps=n.snp, tot.snps=90000, type="all"));
#
#     indices <- which(snps.75000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
#     df <- rbind(df, data.frame(probs=length(which(snps.75000.matched[indices]))/length(indices), snps=n.snp, tot.snps=75000, type="all"));
#
#     indices <- which(snps.50000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
#     df <- rbind(df, data.frame(probs=length(which(snps.50000.matched[indices]))/length(indices), snps=n.snp, tot.snps=50000, type="all"));
#
#     indices <- which(snps.25000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
#     df <- rbind(df, data.frame(probs=length(which(snps.25000.matched[indices]))/length(indices), snps=n.snp, tot.snps=25000, type="all"));
#
#     indices <- which(snps.10000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
#     df <- rbind(df, data.frame(probs=length(which(snps.10000.matched[indices]))/length(indices), snps=n.snp, tot.snps=10000, type="all"));
#
#     indices <- which(snps.5000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
#     df <- rbind(df, data.frame(probs=length(which(snps.5000.matched[indices]))/length(indices), snps=n.snp, tot.snps=5000, type="all"));
#
#     ##browser()
# }
#
# ggplot(df, aes(snps, probs, color=as.factor(tot.snps), pch=type))+geom_point();
#
# pdf("snps.vs.probs.10000.pdf");
# plot(n.snps, output);
# dev.off();
