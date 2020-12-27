### Modules:
'''
import sys
sys.path.extend(['/home/tbarata/CSCs/human_ts_models/projects/csc_devel/src', '/home/tbarata/cobamp/src', '/home/tbarata/troppo/src'])
'''
import os
import matplotlib.pyplot as plt
import seaborn as sns
from urllib.request import urlretrieve
from scaleGenExpFunct import scaleGenExp, GetThresAndRcScores, GetDictOfReacScoresDF, EucDstWithinSampleGroups, SimEucDistWithinSampleGroups, GetPvalueDfs, BestHighestNbGroups


### Paths
#BaseDir = '/home/tbarata/CSCs/human_ts_models/projects/csc_devel'
BaseDir = '/home/tania/CSCs/human_ts_models/projects/csc_devel'
modelDir = os.path.join(BaseDir, 'data/MetModels/iHuman')
modelPathSBML, _ = urlretrieve('https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/master/model/Human-GEM.xml')
modelProcPath = os.path.join(modelDir, 'iHuman_processedEMEM.xml')
modelProcInfoPath = os.path.join(modelDir, 'iHuman_processedInfoEMEM.xlsx')
modelMedAdapPath = os.path.join(modelDir, 'iHuman_processed_MedAdapEMEM.xml')
modelInfoPath = os.path.join(modelDir, 'iHumanInfo.xlsx')
HamsPathName = 'HAMs_iHuman.xlsx'  # ham's medium. when univ. model doesn't grow on chosen medium, metabolites from ham's media are added to chosen medium for model to grow
HamsPath = os.path.join(modelDir, HamsPathName)
BaseMediumPathName = 'EMEM_iHuman.xlsx' # chosen medium
BaseMediumPath = os.path.join(modelDir, BaseMediumPathName)
MediumPathName = 'EMEM_iHuman_plus.xls' # chosen medium plus metabolites from ham's medium that allow universal model to grow
MediumPath = os.path.join(modelDir, MediumPathName)
BiomassID = 'biomass_human'

### Apply min-max normalization to gene expression of each study:
# also return generic model and model adapted for medium composition
model, modelMedAdap = scaleGenExp(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, HamsPath, BaseMediumPath, MediumPath, BiomassID, BaseDir)

### Get and save reaction scores for all combinations of thresholds + and/or rules for all studies' samples in dictionary
# dict with studies samples as keys and values are dictionaries,
# inner dictionaries have threshold + And/Or rules as keys and other dicts as values.
# lower dictionaries have reactions as keys and reaction scores as values:
GetThresAndRcScores(BaseDir, model, GenExpAllName='MinMaxNormExp_AllGenes.csv', GenExpMetbName='MinMaxNormExp_MetbGenes.csv')

### Get/save dictionary with threshold+and/or rule combinations as keys and dataframe as value,
### dataframe has react ids as rows and samples as columns, where values are reactions scores:
GetDictOfReacScoresDF(BaseDir)

### Get/save dataframe with thresholds as columns and groups as rows, where values are average of real euclidean distance between reactions scores of donors of same group:
# group is composed of dif donors/cell lines of same cell type (CSCs/CCs,etc) and same study
# Also get dictionary with group name and list of corresponding donors/cell lines' ids
GrpDct = EucDstWithinSampleGroups(BaseDir)

### Get/save list of 1000 simulated dataframes where each dataframe is similar to above dataframe but distances are between donors/cell lines of randomly selected groups of same size as original ones.
### and get/save list of 1000 simulated dataframes where each dataframe has 'Trues' when distance in a simulated group is lower than real distance observed for corresponding real group:
SimEucDistWithinSampleGroups(BaseDir, GrpDct, maxSim = 1000)

### Get:
### - dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group.
###   a.k.a. dataframe of p-values
### - dataframe with 'Trues' when a group pvalue is below a user defined threshold
PvalSimDstBellowObs, IsBellow, usrVal = GetPvalueDfs(BaseDir, usrVal=1E-3)

### - Save dataframe with threshold combinations that have highest number of groups with p-value bellow user defined value.
MGrps = list(map(lambda x: 'G'+str(x),range(1,11))) # sample groups with microarray samples
RGrps = list(map(lambda x: 'G'+str(x),range(11,23))) # sample groups with rna-seq samples
BestHighestNbGroups(IsBellow, BaseDir)

### Barplot of pvalues of each group of samples for each all genes vs only metabolic genes:
# create dataframes with average p-values by groups for all genes vs only metabolic genes:
# allgenes:
allGenesNames = [thr for thr in PvalSimDstBellowObs.columns if 'AllGenes' in thr[0]] # select all genes
allGenesV = PvalSimDstBellowObs[allGenesNames].mean(axis=1)
allGenesStd = PvalSimDstBellowObs[allGenesNames].std(axis=1)
# metbGenesNames:
metbGenesNames = [thr for thr in PvalSimDstBellowObs.columns if 'MetbGenes' in thr[0]] # select only metabolic genes
metbGenesV = PvalSimDstBellowObs[metbGenesNames].mean(axis=1)
metbGenesStd = PvalSimDstBellowObs[metbGenesNames].std(axis=1)
# create graph:
path = os.path.join(BaseDir, 'results/TestGeneScores', 'AllvsMetbGenes.png')
labels = list(PvalSimDstBellowObs.index)
x = np.arange(len(labels)) # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots(figsize=(20,12))
mmp = ax.bar(x - width/2, allGenesV, width, label = 'All Genes')
msp = ax.bar(x + width/2, metbGenesV, width, label = 'Metabolic Genes')
ax.set_ylabel('# simulations with dist bellow observed dist / # simulations', fontsize= 14)
ax.set_xticks(x)
ax.tick_params(axis='x', which='minor', labelsize=16)
ax.tick_params(axis='y', which='minor', labelsize=13)
ax.set_xticklabels(labels)
ax.legend(prop={'size': 16})
#autolabel(mmp, ax)
#autolabel(msp, ax)
fig.tight_layout()
plt.show()
plt.savefig(path)

### Barplot of pvalues of each group of samples for local1 vs local1B:
# create dataframes with average p-values by groups for all genes vs only metabolic genes:
# local1:
allGenesNames = [thr for thr in PvalSimDstBellowObs.columns if 'local1_' in thr[0]] # select all genes
allGenesV = PvalSimDstBellowObs[allGenesNames].mean(axis=1)
allGenesStd = PvalSimDstBellowObs[allGenesNames].std(axis=1)
# local1B:
metbGenesNames = [thr for thr in PvalSimDstBellowObs.columns if 'local1B' in thr[0]] # select only metabolic genes
metbGenesV = PvalSimDstBellowObs[metbGenesNames].mean(axis=1)
metbGenesStd = PvalSimDstBellowObs[metbGenesNames].std(axis=1)
# create graph:
path = os.path.join(BaseDir, 'results/TestGeneScores', 'AllvsMetbGenes.png')
labels = list(PvalSimDstBellowObs.index)
x = np.arange(len(labels)) # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots(figsize=(20,12))
mmp = ax.bar(x - width/2, allGenesV, width, label = 'local1')
msp = ax.bar(x + width/2, metbGenesV, width, label = 'local1B')
ax.set_ylabel('# simulations with dist bellow observed dist / # simulations', fontsize= 14)
ax.set_xticks(x)
ax.tick_params(axis='x', which='minor', labelsize=16)
ax.tick_params(axis='y', which='minor', labelsize=13)
ax.set_xticklabels(labels)
ax.legend(prop={'size': 16})
#autolabel(mmp, ax)
#autolabel(msp, ax)
fig.tight_layout()
plt.show()
plt.savefig(path)

### Barplot of pvalues of each group of samples for each and/or rule:
# create dataframes with average p-values by groups for each and/or rule:
# minmax:
minmaxNames = [thr for thr in PvalSimDstBellowObs.columns if thr[1]=='minmax'] # select local 1 threshold combinations
minmaxV = PvalSimDstBellowObs[minmaxNames].mean(axis=1)
minmaxStd = PvalSimDstBellowObs[minmaxNames].std(axis=1)
# minsum:
minsumNames = [thr for thr in PvalSimDstBellowObs.columns if thr[1]=='minsum'] # select local 1 threshold combinations
minsumV = PvalSimDstBellowObs[minsumNames].mean(axis=1)
minsumStd = PvalSimDstBellowObs[minsumNames].std(axis=1)
# create graph:
path = os.path.join(BaseDir, 'results/TestGeneScores', 'AndOrStrategies.png')
labels = list(PvalSimDstBellowObs.index)
x = np.arange(len(labels)) # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots(figsize=(20,12))
mmp = ax.bar(x - width/2, minmaxV, width, label = 'min-max', yerr=minmaxStd)
msp = ax.bar(x + width/2, minsumV, width, label = 'min-sum', yerr=minsumStd)
ax.set_ylabel('# simulations with dist bellow observed dist / # simulations', fontsize= 14)
ax.set_xticks(x)
ax.tick_params(axis='x', which='minor', labelsize=16)
ax.tick_params(axis='y', which='minor', labelsize=13)
ax.set_xticklabels(labels)
ax.legend(prop={'size': 16})
#autolabel(mmp, ax)
#autolabel(msp, ax)
fig.tight_layout()
plt.show()
plt.savefig(path)
## boxplot:
fig, ax = plt.subplots(figsize=(20,12))
BoxDf = PvalSimDstBellowObs.copy()
BoxDf['MetbGenes_local2_2_3.0_2.0', 'Group'] = BoxDf.index
BoxDf = pd.melt(BoxDf, col_level=1, id_vars=['Group'], var_name=['And/Or Rule'])
BoxDf.rename({'value':'# simulations with dist bellow observed dist / # simulations'}, axis=1, inplace=True)
ax = sns.boxplot(x='Group', y='# simulations with dist bellow observed dist / # simulations', hue='And/Or Rule', data=BoxDf)
fig.tight_layout()
plt.show()
plt.savefig(path)






--
fig, ax = plt.subplots(figsize=(20,12))
BoxDf = PvalSimDstBellowObs.copy()
BoxDf['Group'] = BoxDf.index
BoxDf = pd.melt(BoxDf, id_vars=['Group'], var_name=['Threshold Strategy'])
BoxDf['Threshold Strategy'] = BoxDf['Threshold Strategy'].apply(lambda x: x.split('_')[0])
BoxDf.rename({'value':'# simulations with dist bellow observed dist / # simulations'}, axis=1, inplace=True)
ax = sns.boxplot(x='Group', y='# simulations with dist bellow observed dist / # simulations', hue='Threshold Strategy', data=BoxDf)
fig.tight_layout()
plt.show()
plt.savefig(path)




### Barplot of pvalues of each group of samples for each threshold strategy:
# create dataframes with average p-values by groups for each threshold strategy:
# local 1:
local1Names = [thr for thr in PvalSimDstBellowObs.columns if 'local1_' in thr[0]] # select local 1 threshold combinations
local1V = PvalSimDstBellowObs[local1Names].mean(axis=1)
local1Std = PvalSimDstBellowObs[local1Names].std(axis=1)
# local 1B:
local1BNames = [thr for thr in PvalSimDstBellowObs.columns if 'local1B' in thr[0]] # select local 1B threshold combinations
local1BV = PvalSimDstBellowObs[local1BNames].mean(axis=1)
local1BStd = PvalSimDstBellowObs[local1BNames].std(axis=1)
# local2:
local2Names = [thr for thr in PvalSimDstBellowObs.columns if 'local2_' in thr[0]] # select local 2 threshold combinations
local2V = PvalSimDstBellowObs[local2Names].mean(axis=1)
local2Std = PvalSimDstBellowObs[local2Names].std(axis=1)
# local2:
local2BNames = [thr for thr in PvalSimDstBellowObs.columns if 'local2B' in thr[0]] # select local 2B threshold combinations
local2BV = PvalSimDstBellowObs[local2BNames].mean(axis=1)
local2BStd = PvalSimDstBellowObs[local2BNames].std(axis=1)
# global:
globalNames = [thr for thr in PvalSimDstBellowObs.columns if 'global' in thr[0]] # select local 1 threshold combinations
globalV = PvalSimDstBellowObs[globalNames].mean(axis=1)
globalStd = PvalSimDstBellowObs[globalNames].std(axis=1)
# create barplot:
path = os.path.join(BaseDir, 'results/TestGeneScores', 'thresholdStrategies.png')
labels = list(PvalSimDstBellowObs.index)
x = np.arange(len(labels))*2 # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots(figsize=(20,12))
gp = ax.bar(x - 2*width, globalV, width, label = 'global')
l1p = ax.bar(x - width, local1V, width, label = 'local1')
l1Bp = ax.bar(x, local1BV, width, label = 'local1B')
l2p = ax.bar(x + width, local2V, width, label = 'local2')
l2Bp = ax.bar(x + 2*width, local2BV, width, label = 'local2B')
ax.set_ylabel('# simulations with dist bellow observed dist / # simulations', fontsize= 14)
ax.set_xticks(x)
ax.tick_params(axis='x', which='minor', labelsize=16)
ax.tick_params(axis='y', which='minor', labelsize=13)
ax.set_xticklabels(labels)
ax.legend(prop={'size': 16})
#autolabel(gp, ax)
#autolabel(l1p, ax)
#autolabel(l2p, ax)
fig.tight_layout()
plt.show()
plt.savefig(path)










### PCA of threshold strategies global, local1 and local2):
# create dataframe with all samples + threshold combinations on rows and reactions on column, with reaction scores as values - where NAs are replaced by median:
CreateAllSmpThrCmbDf(BaseDir)
path = os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'AllStudiesAllThresReacScoreDf.pkl')
infile = open(path,'rb')
AlThrSmp = pickle.load(infile)
infile.close()
FAlThrSmp = AlThrSmp.copy()
FAlThrSmp['ThresholdStrategy'] = FAlThrSmp['Threshold'].apply(lambda y: y[:6]) # add column with threshold strategy
# from previous dataframe get only reaction scores as numpy array:
FAlThrSmpValues = FAlThrSmp.iloc[:,:(FAlThrSmp.shape[1] - 3)].values
title = 'PCA for all threshold strategies'
thresholds = list(FAlThrSmp['Threshold'].apply(lambda y: y[:6]).unique())
path=os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAallThrStrategies_PC1PC2.png')
createPCAgraph(valuesDf=FAlThrSmpValues, completeDf=FAlThrSmp, index=FAlThrSmp.index, title=title, classes= thresholds, classesName='ThresholdStrategy', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
# now just for rnaseq samples:
FAlThrSmpR = FAlThrSmp[pd.Series(FAlThrSmp.index, index = FAlThrSmp.index).apply(lambda y: y.startswith('R'))]
FAlThrSmpRValues = FAlThrSmpR.iloc[:,:(FAlThrSmpR.shape[1] - 3)].values
title = 'PCA for all threshold strategies - RNAseq'
path=os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAallThrStrategies_rnaseq_PC1PC2.png')
createPCAgraph(valuesDf=FAlThrSmpRValues, completeDf=FAlThrSmpR, index=FAlThrSmpR.index, title=title, classes= thresholds, classesName='ThresholdStrategy', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
path=os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAallThrStrategies_rnaseq_PC3PC4.png')
createPCAgraph(valuesDf=FAlThrSmpRValues, completeDf=FAlThrSmpR, index=FAlThrSmpR.index, title=title, classes= thresholds, classesName='ThresholdStrategy', where2Save=path, legend=True, Ncomp=4, comp2plot=['PC3', 'PC4'])
# now just for microarray samples:
FAlThrSmpM = FAlThrSmp[pd.Series(FAlThrSmp.index, index = FAlThrSmp.index).apply(lambda y: y.startswith('M'))]
FAlThrSmpMValues = FAlThrSmpM.iloc[:,:(FAlThrSmpM.shape[1] - 3)].values
title = 'PCA for all threshold strategies - microarray'
path=os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAallThrStrategies_microarray_PC1PC2.png')
createPCAgraph(valuesDf=FAlThrSmpMValues, completeDf=FAlThrSmpM, index=FAlThrSmpM.index, title=title, classes= thresholds, classesName='ThresholdStrategy', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])

### PCA of AND/OR rules
# create dataframe with all samples + threshold combination on rows and reactions on column, with reaction scores as values - where NAs are replaced by median:
path = os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'AllStudiesAllThresReacScoreDf.pkl')
infile = open(path,'rb')
AlThrSmp = pickle.load(infile)
infile.close()
# from previous dataframe get only reaction scores as numpy array:
AlThrSmpValues = AlThrSmp.iloc[:,:(AlThrSmp.shape[1] - 2)].values
title = 'PCA for all AND/OR rules'
aor = list(AlThrSmp['And_OR_rule'].unique())
path = os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAandorRules_PC1PC2.png')
createPCAgraph(valuesDf=AlThrSmpValues, completeDf=AlThrSmp, index=AlThrSmp.index, title=title, classes= aor, classesName='And_OR_rule', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
# now just for rnaseq samples:
AlThrSmpR = AlThrSmp[pd.Series(AlThrSmp.index, index = AlThrSmp.index).apply(lambda y: y.startswith('R'))]
AlThrSmpRValues = AlThrSmpR.iloc[:,:(AlThrSmpR.shape[1] - 2)].values
title = 'PCA for all AND/OR rules - RNAseq'
path = os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAandorRules_rnaseq_PC1PC2.png')
createPCAgraph(valuesDf=AlThrSmpRValues, completeDf=AlThrSmpR, index=AlThrSmpR.index, title=title, classes= aor, classesName='And_OR_rule', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
# now just for microarray samples:
AlThrSmpM = AlThrSmp[pd.Series(AlThrSmp.index, index = AlThrSmp.index).apply(lambda y: y.startswith('M'))]
AlThrSmpMValues = AlThrSmpM.iloc[:,:(AlThrSmpM.shape[1] - 2)].values
title = 'PCA for all AND/OR rules - microarray'
path = os.path.join(BaseDir, 'Models2TestGeneScores', 'ThresholdByStudy', 'PCAandorRules_microarray_PC1PC2.png')
createPCAgraph(valuesDf=AlThrSmpMValues, completeDf=AlThrSmpM, index=AlThrSmpM.index, title=title, classes= aor, classesName='And_OR_rule', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
