library(data.table);
library(ggplot2);

truth <- fread("family.batch1.barcodes.txt")

header <- c("HG00146", "HG00147", "HG00500", "HG00501", "HG00502", "HG00512", "HG00513", "HG00514", "HG00524")

trio.rst <- fread("test.plp.demux.best");
trio.matched <- match(trio.rst$SNG.1ST, header)==truth$V2+1
trio.hits <- length(trio.matched);

## look at the downsmapling results
snps.75000.rst <- fread("snps.75000.best");
snps.75000.matched <- match(snps.75000.rst$SNG.1ST, header)==truth$V2[match(snps.75000.rst$BARCODE,truth$V1)]+1;
snps.75000.hits <- length(which(snps.75000.matched));

snps.50000.rst <- fread("snps.50000.best");
snps.50000.matched <- match(snps.50000.rst$SNG.1ST, header)==truth$V2[match(snps.50000.rst$BARCODE,truth$V1)]+1;
snps.50000.hits <- length(which(snps.50000.matched));

snps.25000.rst <- fread("snps.25000.best");
snps.25000.matched <- match(snps.25000.rst$SNG.1ST, header)==truth$V2[match(snps.25000.rst$BARCODE,truth$V1)]+1;
snps.25000.hits <- length(which(snps.25000.matched));

snps.10000.rst <- fread("snps.10000.best");
snps.10000.matched <- match(snps.10000.rst$SNG.1ST, header)==truth$V2[match(snps.10000.rst$BARCODE,truth$V1)]+1;
snps.10000.hits <- length(which(snps.10000.matched));

snps.5000.rst <- fread("snps.5000.best");
snps.5000.matched <- match(snps.5000.rst$SNG.1ST, header)==truth$V2[match(snps.5000.rst$BARCODE,truth$V1)]+1;
snps.5000.hits <- length(which(snps.5000.matched));

## let's plot average number of snps per cell by the probability of assigning

probs <- c(snps.5000.hits/6145,
snps.10000.hits/6145,
snps.25000.hits/6145,
snps.50000.hits/6145,
snps.75000.hits/6145,
trio.hits/6145);

snps <- c(mean(snps.5000.rst$N.SNP),
mean(snps.10000.rst$N.SNP),
mean(snps.25000.rst$N.SNP),
mean(snps.50000.rst$N.SNP),
mean(snps.75000.rst$N.SNP),
mean(trio.rst$N.SNP))

df <- data.frame(probs=probs, snps=snps, tot.snps=c(5000,10000,25000,50000,75000,90000), type="averages");

##pdf("snps.vs.probs.pdf"); plot(snps, probs); dev.off();

## let's make some bins of SNPs to see the probability of correctly identifying

##matched <- c(trio.matched, snps.75000.matched, snps.50000.matched, snps.25000.matched, snps.10000.matched, snps.5000.matched)
##snps <- c(trio.rst$N.SNP, snps.75000.rst$N.SNP, snps.50000.rst$N.SNP, snps.25000.rst$N.SNP, snps.10000.rst$N.SNP, snps.5000.rst$N.SNP);

##output <- NULL;

n.snps <- seq(1,100,1);

for(n.snp in n.snps) {
    indices <- which(trio.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
    df <- rbind(df, data.frame(probs=length(which(trio.matched[indices]))/length(indices), snps=n.snp, tot.snps=90000, type="all"));

    indices <- which(snps.75000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
    df <- rbind(df, data.frame(probs=length(which(snps.75000.matched[indices]))/length(indices), snps=n.snp, tot.snps=75000, type="all"));

    indices <- which(snps.50000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
    df <- rbind(df, data.frame(probs=length(which(snps.50000.matched[indices]))/length(indices), snps=n.snp, tot.snps=50000, type="all"));

    indices <- which(snps.25000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
    df <- rbind(df, data.frame(probs=length(which(snps.25000.matched[indices]))/length(indices), snps=n.snp, tot.snps=25000, type="all"));

    indices <- which(snps.10000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
    df <- rbind(df, data.frame(probs=length(which(snps.10000.matched[indices]))/length(indices), snps=n.snp, tot.snps=10000, type="all"));

    indices <- which(snps.5000.rst$N.SNP == n.snp); ## <= n.snp & trio.snps > n.snp-5);
    df <- rbind(df, data.frame(probs=length(which(snps.5000.matched[indices]))/length(indices), snps=n.snp, tot.snps=5000, type="all"));

    ##browser()
}

a <- ggplot(df, aes(snps, probs, color=as.factor(tot.snps), pch=type))+geom_point();
ggsave(a, file="snps.vs.probs.pdf")

# pdf("snps.vs.probs.10000.pdf");
# plot(n.snps, output);
# dev.off();
