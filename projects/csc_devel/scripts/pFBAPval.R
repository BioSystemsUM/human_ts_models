library(ggplot2)
data(warpbreaks)
ggplot(warpbreaks, aes(x = tension, y = breaks)) +
  geom_bar(stat = "identity") + facet_wrap(.~wool) +
  ggtitle("Breaks for wool A and B")
ggplot(warpbreaks, aes(x = tension, y = breaks, fill = wool)) +
  geom_bar(stat = "identity", position = "dodge") + 
  ggtitle("Breaks for wool A and B") + ylab("Mean breaks")
warpbreaks
#######################
library(limma)
setwd('/home/tania/CSCs/human_ts_models/projects/csc_devel/results')
dfN = 'pfbaOrdAbsVal.tab'
sbN = 'subsytemInfo.tab'
df = read.csv(dfN, sep='\t', row.names = 1)
sb = read.csv(sbN, sep='\t', row.names = 1)
sbL = lapply(split(sb, sb$Subsystem), row.names)
sbs = unique(sb$Subsystem)
el = 'Metabolism_of_other_amino_acids'
ts = 'AML'
Lstup = list()
Lstdown = list()
for (ts in colnames(df)){
  tsvup = vector()
  tsvdown = vector()
  for(el in sbs){
    f = df[ts]
    ff = f[order(f[,1], decreasing = T), ,drop=F]
    ind = rownames(ff)%in%sbL[[el]]
    resup = geneSetTest(ind, ff[,1], alternative = 'up', ranks.only = F, nsim = 1000)
    resdown = geneSetTest(ind, ff[,1], alternative = 'down', ranks.only = F, nsim = 1000)
    tsvup[el] = resup
    tsvdown[el] = resdown
  }
  Lstup[[ts]] = tsvup
  Lstdown[[ts]] = tsvdown
}
upDf = do.call(cbind, Lstup)
downDf = do.call(cbind, Lstdown)
# replace pvals not signifficant by 1 (cause log1 = 0) so that they do not appear in bar graph bellow:
thr = 0.01
upDf[upDf>=thr] = 1
downDf[downDf>=thr] = 1
# apply -log and +log (respectively) to up and down enriched:
upDf[] = sapply(upDf, function(x) {-log2(x)}) # [] is to keep row and column names
downDf[] = sapply(downDf, log2)
newDf <- data.frame(upDf) # does a copy of upDf (can be checked with tracemem(newDf) == tracemem(upDf) cause results in False)
newDf[downDf < 0] = downDf[downDf < 0] # put up and down values in same dataframe
newDf = newDf[rowSums(newDf[,1:ncol(newDf)]) != 0,]
library(reshape2)
newDf = cbind('Subsystem' = rownames(newDf), data.frame(newDf, row.names = NULL))
mf = melt(newDf, id.vars='Subsystem', measure.vars=colnames(newDf)[2:ncol(newDf)])
ggplot(mf, aes(x = variable, y = value, fill = Subsystem)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  ggtitle('Enrichment pFBA') + ylab('(+/-)log(p-value)') +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))