draw.alphabeta <- function(alphas, betas, p0, C0, ns, fname) {
    df <- data.frame()
    for(alpha in alphas) {
        for(beta in betas) {
            for(n in ns) {
                sol <- uniroot(function(x) { (1-(n-1)/n*beta)*(1-(1-p0)^(x/C0)) - p0 }, c(C0-1,C0*n+1))
                C <- sol$root
                ED <- C0*( 1-(1-p0)^(C/C0) )/(0-log(1-p0))        
                ndoub <- (1-(n-1)/n*beta)*(1-(1-p0)^(C/C0))*ED
                nsing <- (1-alpha)*C0*(1-p0)^(C/C0)*(1-(1-p0)^(C/C0))/(-log(1-p0))
                df <- rbind(df, data.frame(N.SAMPLEs=n, TOTAL.CELLs=C, TOTAL.SINGLETs=nsing, TOTAL.DOUBLETs=ndoub, TOTAL.DROPLETs=ED, CELLs.PER.SAMPLE = C/n, SINGLETs.PER.SAMPLE = nsing/n, NON.SINGLETs = ED-nsing, TypeIErr = alpha, Power = beta));
            }
        }
    }
    return(df)


}

##plot.singlets <- function(p.thres, beta = 0.95, alpha = 0.01, C0 = 1000, p0 = 0.01) {
df0 <- data.frame()
df0 <- draw.alphabeta(c(0.005,0.01,0.02,0.05,0.1), c(0.95), 0.01, 1000, 1:100)
pdf("v2.alphabeta.pdf",width=6,height=4)
print( ggplot(df0,aes(N.SAMPLEs,TOTAL.SINGLETs,colour=as.factor(TypeIErr*100)))+geom_line(size=0.7) + labs(colour="% singlets\nclassified as\nnon-singlets\n(95% power)") )
#dev.off()

df0 <- data.frame()
df0 <- draw.alphabeta(c(0.02), c(0.85, 0.90, 0.95, 0.98, 0.99), 0.01, 1000, 1:100)
#pdf("v2.alpha.pdf",width=6,height=4)
print( ggplot(df0,aes(N.SAMPLEs,TOTAL.SINGLETs,colour=as.factor(Power)))+geom_line(size=0.7) + labs(colour="% true multiplets\nclassified as\nnon-singlets\n(2% Type.I.Err)") )
dev.off()
#print(df)
#}

