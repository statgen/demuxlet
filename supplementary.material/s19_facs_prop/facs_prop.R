library(ggplot2);

load("lupus.facs.props.rda");
load("lupus.scRNAseq.prop.rda");

out.df <- data.frame(scRNA.seq=df.props.complete[,1], FACS=bcells.facs.ctrl, cell="CD19+");
out.df <- rbind(out.df,
  data.frame(scRNA.seq=df.props.complete[,2], FACS=monocyte.facs.ctrl, cell="CD14+"));
out.df <- rbind(out.df,
  data.frame(scRNA.seq=df.props.complete[,3], FACS=cd4.facs.ctrl, cell="CD4+"));
out.df <- rbind(out.df,
  data.frame(scRNA.seq=df.props.complete[,4], FACS=cd8.facs.ctrl, cell="CD8+"));
##out.df <- rbind(out.df,
##  data.frame(scRNA.seq=df.props.complete[,8], FACS=nk.facs.ctrl, cell="CD56+"));

ctrl.ggplot <- ggplot(data=out.df, aes(scRNA.seq, FACS, color=cell))+geom_point()
ggsave(ctrl.ggplot, file="ctrl.ggplot.pdf");

out.stim.df <- data.frame(scRNA.seq=df.props.complete[,1], FACS=bcells.facs.stim, cell="CD19+");
out.stim.df <- rbind(out.stim.df,
  data.frame(scRNA.seq=df.props.complete[,2]+df.props.complete[,5]+df.props.complete[,6], FACS=monocyte.facs.stim, cell="CD14+"));
out.stim.df <- rbind(out.stim.df,
  data.frame(scRNA.seq=df.props.complete[,3], FACS=cd4.facs.stim, cell="CD4+"));
out.stim.df <- rbind(out.stim.df,
  data.frame(scRNA.seq=df.props.complete[,4], FACS=cd8.facs.stim, cell="CD8+"));
##out.stim.df <- rbind(out.stim.df,
##  data.frame(scRNA.seq=df.props.complete[,8], FACS=nk.facs.stim, cell="CD56+"));

stim.ggplot <- ggplot(data=out.stim.df, aes(scRNA.seq, FACS, color=cell))+geom_point()
ggsave(stim.ggplot, file="stim.ggplot.pdf");
