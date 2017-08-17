library(data.table);

##################################
## BATCH 1 ANALYSIS
##################################

batch1.w1 <- fread("cramore/jy-a.seqaggr.eagle.10.sm.best");
batch1.w2 <- fread("cramore/jy-b.seqaggr.eagle.10.sm.best");
batch1.w3 <- fread("cramore/jy-c.seqaggr.eagle.10.sm.best");

batch1.w1.llk.singlets <- which(batch1.w1$SNG.LLK1-batch1.w1$LLK12 > 1);
batch1.w2.llk.singlets <- which(batch1.w2$SNG.LLK1-batch1.w2$LLK12 > 1);
batch1.w3.llk.singlets <- which(batch1.w3$SNG.LLK1-batch1.w3$LLK12 > 1);

batch1.w1.llk.doublets <- which(batch1.w1$LLK12-batch1.w1$SNG.LLK1 > 1);
batch1.w2.llk.doublets <- which(batch1.w2$LLK12-batch1.w2$SNG.LLK1 > 1);
batch1.w3.llk.doublets <- which(batch1.w3$LLK12-batch1.w3$SNG.LLK1 > 1);

batch1.w1.llk.amb <- intersect(which(batch1.w1$LLK12-batch1.w1$SNG.LLK1 < 1),which(batch1.w1$SNG.LLK1-batch1.w1$LLK12 < 1));
batch1.w2.llk.amb <- intersect(which(batch1.w2$LLK12-batch1.w2$SNG.LLK1 < 1),which(batch1.w2$SNG.LLK1-batch1.w2$LLK12 < 1));
batch1.w3.llk.amb <- intersect(which(batch1.w3$LLK12-batch1.w3$SNG.LLK1 < 1),which(batch1.w3$SNG.LLK1-batch1.w3$LLK12 < 1));

batch1.w1.post.singlets <- which(batch1.w1$PRB.DBL < 0.1)
batch1.w2.post.singlets <- which(batch1.w2$PRB.DBL < 0.1)
batch1.w3.post.singlets <- which(batch1.w3$PRB.DBL < 0.1)

batch1.w1.post.doublets <- which(batch1.w1$PRB.DBL > 0.9)
batch1.w2.post.doublets <- which(batch1.w2$PRB.DBL > 0.9)
batch1.w3.post.doublets <- which(batch1.w3$PRB.DBL > 0.9)

batch1.w1.post.amb <- intersect(which(batch1.w1$PRB.DBL < 0.9), which(batch1.w1$PRB.DBL > 0.1))
batch1.w2.post.amb <- intersect(which(batch1.w2$PRB.DBL < 0.9), which(batch1.w2$PRB.DBL > 0.1))
batch1.w3.post.amb <- intersect(which(batch1.w3$PRB.DBL < 0.9), which(batch1.w3$PRB.DBL > 0.1))

pdf("batch1.w1.llk.vs.post.pdf");
plot(batch1.w1$LLK12-batch1.w1$SNG.LLK1, batch1.w1$PRB.DBL, xlim=c(-10,10))
dev.off();

pdf("batch1.w2.llk.vs.post.pdf");
plot(batch1.w2$LLK12-batch1.w2$SNG.LLK1, batch1.w2$PRB.DBL, xlim=c(-10,10))
dev.off();

pdf("batch1.w3.llk.vs.post.pdf");
plot(batch1.w3$LLK12-batch1.w3$SNG.LLK1, batch1.w3$PRB.DBL, xlim=c(-10,10))
dev.off();

##################################
## BATCH 4 ANALYSIS
##################################


batch4 <- fread("cramore/batch4.cramore.test0.1.sm.sorted.2.filtered.best");

batch4.llk.singlets <- which(batch4$SNG.LLK1-batch4$LLK12 > 1);

batch4.llk.doublets <- which(batch4$LLK12-batch4$SNG.LLK1 > 1);

batch4.llk.amb <- intersect(which(batch4$LLK12-batch4$SNG.LLK1 < 1),which(batch4$SNG.LLK1-batch4$LLK12 < 1));

batch4.post.singlets <- which(batch4$PRB.DBL < 0.1)

batch4.post.doublets <- which(batch4$PRB.DBL > 0.9)

batch4.post.amb <- intersect(which(batch4$PRB.DBL < 0.9), which(batch4$PRB.DBL > 0.1))

pdf("batch4.llk.vs.post.pdf");
plot(batch4$LLK12-batch4$SNG.LLK1, batch4$PRB.DBL, xlim=c(-10,10))
dev.off();

pdf("batch4.llk.vs.post.pdf");
plot(batch4$LLK12-batch4$SNG.LLK1, batch4$PRB.DBL, xlim=c(-10,10))
dev.off();

pdf("batch4.w3.llk.vs.post.pdf");
plot(batch4$LLK12-batch4$SNG.LLK1, batch4$PRB.DBL, xlim=c(-10,10))
dev.off();


##################################
## BATCH 5 ANALYSIS
##################################

batch5.w1 <- fread("cramore/Ye10.1.sm.sorted.2.filtered.best");
batch5.w2 <- fread("cramore/Ye20.1.sm.sorted.2.filtered.best");

batch5.w1.llk.singlets <- which(batch5.w1$SNG.LLK1-batch5.w1$LLK12 > 1);
batch5.w2.llk.singlets <- which(batch5.w2$SNG.LLK1-batch5.w2$LLK12 > 1);

batch5.w1.llk.doublets <- which(batch5.w1$LLK12-batch5.w1$SNG.LLK1 > 1);
batch5.w2.llk.doublets <- which(batch5.w2$LLK12-batch5.w2$SNG.LLK1 > 1);

batch5.w1.llk.amb <- intersect(which(batch5.w1$LLK12-batch5.w1$SNG.LLK1 < 1),which(batch5.w1$SNG.LLK1-batch5.w1$LLK12 < 1));
batch5.w2.llk.amb <- intersect(which(batch5.w2$LLK12-batch5.w2$SNG.LLK1 < 1),which(batch5.w2$SNG.LLK1-batch5.w2$LLK12 < 1));

batch5.w1.post.singlets <- which(batch5.w1$PRB.DBL < 0.1)
batch5.w2.post.singlets <- which(batch5.w2$PRB.DBL < 0.1)

batch5.w1.post.doublets <- which(batch5.w1$PRB.DBL > 0.9)
batch5.w2.post.doublets <- which(batch5.w2$PRB.DBL > 0.9)

batch5.w1.post.amb <- intersect(which(batch5.w1$PRB.DBL < 0.9), which(batch5.w1$PRB.DBL > 0.1))
batch5.w2.post.amb <- intersect(which(batch5.w2$PRB.DBL < 0.9), which(batch5.w2$PRB.DBL > 0.1))

pdf("batch5.w1.llk.vs.post.pdf");
plot(batch5.w1$LLK12-batch5.w1$SNG.LLK1, batch5.w1$PRB.DBL, xlim=c(-10,10))
dev.off();

pdf("batch5.w2.llk.vs.post.pdf");
plot(batch5.w2$LLK12-batch5.w2$SNG.LLK1, batch5.w2$PRB.DBL, xlim=c(-10,10))
dev.off();
