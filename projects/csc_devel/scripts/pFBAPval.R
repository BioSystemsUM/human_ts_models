### Packages:
library(limma)
source('~/CSCs/human_ts_models/projects/csc_devel/src/paths.py') # path to paths.py file. It is hardcoded.

setwd(paste(projFld, BaseDir, sep= '/'))

rkGenSetEnr = function (dfN, outp) {
  ### Ranked reaction set enrichment analysis:
  sbN = 'support/Pathways/subsytemInfo.tab'
  # load dataframe with difference between absolute values of pfba fluxes of CSCs and CCs:
  df = read.csv(dfN, sep='\t', row.names = 1, check.names = F)
  # load dataframe with generic model reaction ids and corresponding subsystems:
  sb = read.csv(sbN, sep='\t', row.names = 1)
  # do gene-set enrichment analysis for all tissues and subsytems:
  sbL = lapply(split(sb, sb$Subsystem), row.names) # list of subsystems, each has a vector with corresponding reaction ids
  sbs = unique(sb$Subsystem) # vector with all unique subsytems
  Lstup = list()
  for (ts in colnames(df)){ # for each tissue
    tsvup = vector()
    tsvdown = vector()
    f = df[ts]
    ff = f[order(f[,1], decreasing = T), ,drop=F] # order dataframe
    for(el in sbs){ # for each subsystem
      # bollean vector where 'Trues' exist when in the ordered vector of pfba diff. (from a tissue)
      # a reaction belongs to that subsystem:
      ind = rownames(ff)%in%sbL[[el]]
      # ranked gene-set enrichment test:
      resup = geneSetTest(ind, ff[,1], alternative = 'up', ranks.only = F, nsim = 1000) # pvalue for subsytem enrichment in top of list (CSCs)
      tsvup[el] = resup # pvalue is saved to a vector and labeled with corresponding substem
    }
    Lstup[[ts]] = tsvup # vector is saved into a list where each element is a tissue
  }
  upDf = do.call(cbind, Lstup) # binds list of all tissues into dataframe
  
  applylog = function(up, pathOut) {
    # apply -log to up-enriched:
    up[] = sapply(up, function(x) {-log2(x)}) # [] is to keep row and column names
    new <- data.frame(up) # does a copy of upDf (can be checked with tracemem(newDf) == tracemem(upDf) cause results in False)
    # exclude rows with all 0's (no enrichment):
    new = new[rowSums(new[,1:ncol(new)]) != 0,]
    write.table(new, pathOut, sep='\t', col.names = T, row.names = T, quote = F)
  }
  ## now including just significant subsytems:
  # replace pvals not signifficant by 1 (cause log1 = 0) so that they do not appear in dataframe:
  thr = 0.05
  upDf[upDf >= thr] = 1
  applylog(up = upDf, pathOut = outp)
}
# do ranked reaction set enrichment analysis for difference between pfba values of CSCs vs CCs of different tissues:
dfN = 'support/Simulations/pfbaAbsValDiffTss.tab'
outp = 'support/Simulations/pfbaSubsDiffPvalSign.tab'
rkGenSetEnr(dfN, outp)

