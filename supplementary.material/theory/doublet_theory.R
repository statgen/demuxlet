library(ggplot2);
library(reshape2);

## some theoretical basis for demuxlet

## figure s1
## number of samples
Ns <- seq(1,99,1);

## d0 is the observed doublet rate when u cells are loaded
d0 <- 0.01;
u <- 1000;

#################
## MODEL 1
#################
## lambda1 is the poisson rate parameter of observing a doublet
## probability of P(singlet) = P(0 multiplet) = exp(-lambda1)*lambda1^0/0! = exp(-lambda1);
## probability of multiplet is P(multiplet) = 1 - P(singlet) = 1 - exp(-lambda1)
## but we need to estimate lambda1 from d0, u and x (number of cells)
## the doublet rate grows linearly as the number of cells (x) by 1/Z
Z <- -u/log(1-d0);
##lambda1 <- x/Z;

##################
## MODEL 2
##################

## 750k beads
## capturing u cells gives us some estimate of lambda
## sum_750k(1-P(0 cells)) = u
## P(0 cells) = exp(-lambda2)
## sum_750k(1-exp(-lambda2)) = u
## lambda2 = -log(1-(u/7.5e5))
lambda2 = -log(1-(u/7.5e5));

## number of cells
##Nc <- seq(1000,100000,1000);

## doublet rates
##rs <- Nc*lambda^2/Nc

## undetectable doublet rate
ds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10)

total.cells <- matrix(NA, nrow=length(Ns), ncol=length(ds));
singlet.cells <- matrix(NA, nrow=length(Ns), ncol=length(ds));

for(N.i in 1:length(Ns)) {
  ## detectable doublets
  N <- Ns[N.i];

  ##detectable.doublet.perc <- 1/N;
  ##undetectable.doublet.perc <- 1-1/N;

  for(d.i in 1:length(ds)) {
    d <- ds[d.i];
    ##total.cells[N.i, r.i] <- r/(1e-5*undetectable.doublet.perc);
    total.cells[N.i, d.i] <- -Z*log(1-N*d);
    singlet.cells[N.i, d.i] <- -Z*(1-N*d)*log(1-N*d);
  }
}

colnames(total.cells) <- colnames(singlet.cells) <- as.factor(ds);
rownames(total.cells) <- rownames(singlet.cells) <- as.factor(Ns);

total.cells.melt <- melt(total.cells);
singlet.cells.melt <- melt(singlet.cells);
cells.melt <- rbind(data.frame(total.cells.melt,cell.type="total"), data.frame(singlet.cells.melt,cell.type="singlets"))

colnames(cells.melt) <- c("N.inds", "doublet.rate","N.cells","cell.type")

cells.ggplot <- ggplot(cells.melt, aes(N.inds,N.cells,color=factor(doublet.rate),shape=cell.type))+geom_point()+theme_bw()+scale_x_log10()+scale_y_log10();

ggsave(cells.ggplot, file="doublet.cells.pdf")

## now let's factor in the cost
total.cost <- matrix(NA, nrow=length(Ns), ncol=length(ds));
singlet.cost <- matrix(NA, nrow=length(Ns), ncol=length(ds));
##library.prep.cost <- matrix(NA, nrow=length(Ns), ncol=length(ds));

for(N.i in 1:length(Ns)) {
  ## detectable doublets
  N <- Ns[N.i];

  ##detectable.doublet.perc <- 1/N;
  ##undetectable.doublet.perc <- 1-1/N;

  for(d.i in 1:length(ds)) {
    d <- ds[d.i];
    ##total.cells[N.i, r.i] <- r/(1e-5*undetectable.doublet.perc);
    total.cost[N.i, d.i] <- (1750/total.cells[N.i, d.i])+25000*total.cells[N.i, d.i]/350e6*2000/total.cells[N.i, d.i];
    singlet.cost[N.i, d.i] <- (1750/singlet.cells[N.i, d.i])+25000*total.cells[N.i, d.i]/350e6*2000/singlet.cells[N.i, d.i];
  }
}


total.cost.melt <- melt(total.cost);
singlet.cost.melt <- melt(singlet.cost);
cost.melt <- rbind(data.frame(total.cost.melt,cell.type="total"), data.frame(singlet.cost.melt,cell.type="singlets"))

colnames(cost.melt) <- c("N.inds", "doublet.rate","Cost","cell.type")

cost.ggplot <- ggplot(cost.melt, aes(N.inds,Cost,color=factor(doublet.rate),shape=cell.type))+geom_point()+theme_bw()+scale_x_log10()+scale_y_log10();

ggsave(cost.ggplot, file="doublet.cost.pdf")
