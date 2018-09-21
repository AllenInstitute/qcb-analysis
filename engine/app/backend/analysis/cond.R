rm(list=ls())

library(nlme)
library(rjson)
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
#args = c("fix","cell_seg_vol","dna_seg_vol" )

if (args[1]=='fix') {
  origin = -1
  print('Fixed origin')
} else {
  origin = +1
  print('Non fixed origin')
}

x_lab = paste('x =',args[2])
y_lab = paste('x =',args[3])

features <- fromJSON(file="../../static/results/features.json")
#features <- fromJSON(file="/Users/viana/projects/qcb/qcb-analysis/engine/app/static/results/../../static/results/features.json")

features <- features['DATA'][[1]]

Table <- NULL
for (data in features) {
  Table <- rbind(Table,data.frame(x=data$x,y=data$y))
}

Table <- within(Table,x<-x*(nrow(Table)/sum(Table$x)))
Table <- within(Table,y<-y*(nrow(Table)/sum(Table$y)))

nptsx <- 128

ordered_x_id <- order(Table$x)

ExpandedTable <- NULL
for (i in seq(1,nrow(Table)-nptsx,nptsx)) {
  subTable <- Table[ordered_x_id[i:(i+nptsx-1)],]
  subTable$bin <- i
  ExpandedTable <- rbind(ExpandedTable,subTable)  
}

bin_levels <- unique(ExpandedTable$bin)
bin_labels <- seq(1,length(bin_levels),1)
ExpandedTable$bin <- factor(ExpandedTable$bin,levels=bin_levels,labels=bin_labels)

STest <- data.frame(bin=levels(ExpandedTable$bin), pvalue=rep( NA,length(levels(ExpandedTable$bin))))
STest <- within(STest, sigma<-pvalue)
STest$bin <- factor(STest$bin,levels=bin_labels,labels=bin_labels)

meanTable <- NULL

bin_id = 1
for (b in bin_labels) {
  x <- ExpandedTable$x[which(ExpandedTable$bin==b)]
  y <- ExpandedTable$y[which(ExpandedTable$bin==b)]
  stest <- shapiro.test(y)
  STest$sigma[bin_id] = sd(y)
  STest$pvalue[bin_id] = stest$p.value
  bin_id = bin_id + 1
  
  meanTable <- rbind(meanTable,data.frame(x=mean(x),y=mean(y),sdx=sd(x),sdy=sd(y)))
}

ExpandedTable$normal <- rep(ifelse(STest$pvalue>0.05,"normal","non-normal"),as.numeric(table(ExpandedTable$bin)))

ExpandedTable$sigma <- rep(STest$sigma,as.numeric(table(ExpandedTable$bin)))

if (args[1]=='free') {
  
  lrr <- lm(formula = y~x, data = ExpandedTable)
  lrw <- lm(formula = y~x, data = ExpandedTable, weights=sigma)
  mean_lrr <- lm(formula = y~x, data=meanTable)
  mean_lrw <- lm(formula = y~x, data=meanTable, weights=sdy)
  mean_nlr_pow <- nls(formula = y~a*(x**b), data=meanTable, start=list("a"=0.5,"b"=1))
  
  mean_nlr_log <- tryCatch({
    nls(formula = y~a*log(b*x+1), data=meanTable, start=list("a"=0.6,"b"=1))},
    error=function(condition) {
      print('Log and pow models are similar')
      return(mean_nlr_pow)
    })
  
} else {
  
  lrr <- lm(formula = y~x-1, data = ExpandedTable)
  lrw <- lm(formula = y~x-1, data = ExpandedTable, weights=sigma)
  mean_lrr <- lm(formula = y~x-1, data=meanTable)
  mean_lrw <- lm(formula = y~x-1, data=meanTable, weights=sdy)
  mean_nlr_pow <- nls(formula = y~a*(x**b)-1, data=meanTable, start=list("a"=0.5,"b"=1))
  
  mean_nlr_log <- tryCatch({
    nls(formula = y~a*log(b*x+1)-1, data=meanTable, start=list("a"=0.6,"b"=1))},
    error=function(condition) {
      print('Log and pow models are similar')
      return(mean_nlr_pow)
    })
  
}  

AICres <- AIC(lrr,lrw)

mean_AICres <- AIC(mean_lrr,mean_lrw,mean_nlr_pow,mean_nlr_log)

ExpandedTable <- within(ExpandedTable, y_pred_lrr <- predict(lrr))
ExpandedTable <- within(ExpandedTable, y_pred_lrw <- predict(lrw))

meanTable <- within(meanTable, y_pred_lrr<-predict(mean_lrr))
meanTable <- within(meanTable, residual_lrr<-y-y_pred_lrr)
meanTable <- within(meanTable, y_pred_lrw<-predict(mean_lrw))
meanTable <- within(meanTable, residual_lrw<-y-y_pred_lrw)
meanTable <- within(meanTable, y_pred_nl_pow<-predict(mean_nlr_pow))
meanTable <- within(meanTable, residual_nl_pow<-y-y_pred_nl_pow)
meanTable <- within(meanTable, y_pred_nl_log<-predict(mean_nlr_log))
meanTable <- within(meanTable, residual_nl_log<-y-y_pred_nl_log)

f1 <- ggplot(ExpandedTable) + geom_density(aes(y,fill=normal)) +
  facet_wrap(~bin, ncol=1) + xlab("p(y|x)") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text.x=element_blank(), legend.position="none")

f2 <- ggplot(STest) + geom_point(aes(bin,pvalue)) + scale_y_log10() + geom_hline(yintercept=0.05) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle('Shapiro Normality Test')

f3 <- ggplot(STest) + geom_point(aes(bin,sigma)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

f4 <- ggplot(ExpandedTable) + geom_point(aes(x,y), alpha=0.05, col='purple') +
  geom_line(aes(x,y_pred_lrr), col='black') +
  geom_line(aes(x,y_pred_lrw), col='green') +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  theme_bw() + xlab(x_lab) + ylab(y_lab) +
  ggtitle(toString(paste(row.names(AICres),format(AICres$AIC,digits=1))))

f5 <- ggplot(meanTable,aes(x,y)) +
  geom_point(data=ExpandedTable,aes(x,y),alpha=0.05,col='purple') +
  geom_errorbar(aes(ymin=y-sdy,ymax=y+sdy, width=0.25*sdx)) +
  geom_errorbarh(aes(xmin=x-sdx,xmax=x+sdx, height=0.25*sdy)) +
  geom_point(col='red') +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  theme_bw() + xlab(x_lab) + ylab(y_lab)

f6 <- ggplot(meanTable,aes(x,y)) +
  geom_errorbar(aes(ymin=y-sdy,ymax=y+sdy, width=0.25*sdx)) +
  geom_errorbarh(aes(xmin=x-sdx,xmax=x+sdx, height=0.25*sdy)) +
  geom_point(col='red') +
  theme_bw() + xlab(x_lab) + ylab(y_lab) +
  geom_line(aes(x,y_pred_lrr), col='black') +
  geom_line(aes(x,y_pred_lrw), col='green') +
  geom_line(aes(x,y_pred_nl_pow), col='red') +
  geom_line(aes(x,y_pred_nl_log), col='blue') +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  ggtitle(toString(paste(row.names(mean_AICres),format(mean_AICres$AIC,digits=1)))) +
  theme(plot.title=element_text(size=8))

f7 <- ggplot(meanTable) +
  geom_point(aes(y,y_pred_lrr),col='black') +
  geom_point(aes(y,y_pred_lrw),col='green') +
  geom_point(aes(y,y_pred_nl_pow),col='red') +
  geom_point(aes(y,y_pred_nl_log),col='blue') +
  theme_bw() + xlab('y') + ylab('y_pred') +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  geom_abline(intercept=0, slope=1, linetype="dashed")

f8 <- ggplot(meanTable) +
  geom_point(aes(x,residual_lrr),col='black') +
  geom_point(aes(x,residual_lrw),col='green') +
  geom_point(aes(x,residual_nl_pow),col='red') +
  geom_point(aes(x,residual_nl_log),col='blue') +
  theme_bw() + xlab('x') + ylab('residual')

jpeg(file = '../../static/results/analysis.jpg', width=1200, height=800)
grid.arrange(arrangeGrob(f4,f5),f1,arrangeGrob(f2,f3),arrangeGrob(f6,f7,f8),ncol=4,widths=c(2,1,2,2))
dev.off()
