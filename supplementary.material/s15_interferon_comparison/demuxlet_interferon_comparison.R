library(ggplot2);
library(data.table);

t <- fread("t.ifn.txt");
##t.timecourse <- fread("t.ifn.timecourse.txt")
##t.timecourse <- cbind(t.timecourse[,c(1,2)],t.timecourse[,c(14)]-t.timecourse[,c(13)]);
dc <- fread("dc.ifn.txt");

dc.small <- dc[match(t$gene,dc$"#gene"),];

cors <- NULL;

cd4.dscrnaseq <- fread("cd4.dscrnaseq.txt");
cd4.matched <- intersect(cd4.dscrnaseq$featureData.symbol, dc.small$"#gene")

cors[1,1] <- cor(cd4.dscrnaseq$log2FoldChange[match(cd4.matched,cd4.dscrnaseq$featureData.symbol)], dc.small$log2FoldChange[match(cd4.matched, dc.small$"#gene")])

dc.dscrnaseq <- fread("dc.dscrnaseq.txt");
dc.matched <- intersect(dc.dscrnaseq$featureData.symbol, dc.small$"#gene")

cors[1,2] <- cor(dc.dscrnaseq$log2FoldChange[match(dc.matched,dc.dscrnaseq$featureData.symbol)], dc.small$log2FoldChange[match(dc.matched, dc.small$"#gene")])
