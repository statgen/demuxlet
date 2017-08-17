library(reshape2);
library(ggplot2);

#m= # snps
#n= # inds
## assuming all SNPs have allele frequency of 50%
sim <- function(m, n, rep = 1000) {
  nuniq <- 0
  nuniq2 <- 0
  first.match <- NULL;
  for(i in 1:rep) {
    mat1 <- matrix(rbinom(m*n, 1, 0.5), n, m);
    mat2 <- matrix(rbinom(m*n, 1, 0.5), n, m);

    reads <- NULL;

    for(j in 1:m) {
      if(sample(1,1:2)==1) {
        reads <- c(reads, mat1[1,j]);
      } else {
        reads <- c(reads, mat2[1,j]);
      }
    }

    match <- t(apply(mat1,1,"==",reads)) | t(apply(mat2,1,"==",reads));

    if ( nrow(unique(matrix(rbinom(m*n, 2, 0.5), n, m))) == n ) { nuniq <- nuniq + 1}

    if ( length(which(rowSums(match)==m)) == 1 ) {nuniq2 <- nuniq2+1;}

    first.match <- c(first.match, which(rowSums(match)==m)[[1]]);

    ##browser();
  }
  return (list(nuniq/rep, nuniq2/rep, first.match))
}

M <- seq(5,30,5);
N <- seq(8,128,8);

out <- data.frame(); ##matrix(NA, nrow=length(M), ncol=length(N))

for(m.i in 1:length(M)) {
  print(m.i);
  m <- M[m.i];
  for(n.i in 1:length(N)) {
    for(rep in 1:10) {
      n <- N[n.i];
      ##out[m.i,n.i] <- sim(m,n)[[2]];
      out <- rbind(out, c(m,n,sim(m,n)[[2]]));
    }
  }
}

colnames(out) <- c("N.SNPs", "N.inds", "Prob")
out.ggplot <- ggplot(data=out, aes(N.inds, Prob, color=N.SNPs))+geom_boxplot()+theme_bw();

# colnames(out) <- as.factor(N);
# rownames(out) <- as.factor(M);
# out.melt <- melt(out);
# colnames(out.melt) <- c("N.SNPs","N.inds","P")
# sig.indices <- which(out.melt[,3]>.98);
# out.sig <- out.melt[sig.indices,]
# out.sig.unique <- out.sig[!duplicated(out.sig$N.inds),]
#

##out.ggplot <- ggplot(out.melt, aes(N.inds, N.SNPs))+geom_raster(aes(fill=P))+scale_fill_distiller(palette="Spectral")+geom_line(aes(N.inds, N.SNPs), data=out.melt[rownames(out.sig.unique),])+theme_bw()

ggsave(out.ggplot, file="nsnps.vs.ninds.pdf")
