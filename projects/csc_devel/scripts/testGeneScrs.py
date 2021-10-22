### Modules:
import sys
import os
from testGeneScrsFunct import scaleGenExp, GetThresAndRcScores, GetDictOfReacScoresDF, EucDstWithinSampleGroups, SimEucDistWithinSampleGroups, GetPvalueDfs
from BuildPaths import modelProcPath, modelMedAdapPath
from paths import BaseDir, projFld
from testGeneScrsFunct import BoxplotTestThres, PCATestThres, getThrsLwstPval

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Apply min-max normalization to gene expression of each study:
# also return generic model and model adapted for medium composition
model, modelMedAdap = scaleGenExp(modelProcPath, modelMedAdapPath, BaseDir)

### Get and save reaction scores for all combinations of thresholds + and/or rules for all studies' samples in dictionary
# dict with studies samples as keys and values are dictionaries,
# inner dictionaries have threshold + And/Or rules as keys and other dicts as values.
# lower dictionaries have reactions as keys and reaction scores as values:
GetThresAndRcScores(BaseDir, model, GenExpAllName='MinMaxNormExp_AllGenes.csv', GenExpMetbName='MinMaxNormExp_MetbGenes.csv')

### Produce a dictionary with threshold+and/or rule combinations as keys and dataframe as value,
### dataframe has react ids as rows and samples as columns, where values are reactions scores:
GetDictOfReacScoresDF(BaseDir)

### Produce a dataframe with thresholds as columns and groups as rows, where values are average of real euclidean distance between reactions scores of donors of same group:
# group is composed of dif donors/cell lines of same cell type (CSCs/CCs,etc) and same study
GrpDct = EucDstWithinSampleGroups(BaseDir)

### Produce list of 1000 simulated dataframes where each dataframe is similar to abovemetioned dataframe but distances are between donors/cell lines of randomly selected groups of same size as original ones.
### and produce list of 1000 simulated dataframes where each dataframe has 'Trues' when distance in a simulated group is lower than real distance observed for corresponding real group:
MGrps = list(map(lambda x: 'G'+str(x),range(1,11))) # sample groups with microarray samples
RGrps = list(map(lambda x: 'G'+str(x),range(11,23))) # sample groups with rna-seq samples
SimEucDistWithinSampleGroups(BaseDir, GrpDct, RGrps, MGrps, maxSim = 1000)

### Get dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group.
### a.k.a. dataframe of p-values
PvalSimDstBellowObs = GetPvalueDfs(BaseDir)

### Save for sample groups G1 (microarray study) and G17 (rna-seq study) the thresholds with lowest values of abovementioned p-values:
## Chosen sample groups (G1 and G17) are those with cell lines for which we have lethality scores data
getThrsLwstPval(PvalSimDstBellowObs, BaseDir)

### Boxplots of reaction scores to compare dif. threshold strategies:
BoxplotTestThres(BaseDir, PvalSimDstBellowObs, MGrps, RGrps)

### PCA graphs of reaction to compare dif. threshold strategies:
PCATestThres(BaseDir)
