### Modules:
import sys
import os
import cobra
from troppo.tasks.task_io import JSONTaskIO
from statsAndSimFunct import MdCmpPlot, MdCmpHeatMap, pFBA, prepPfba2PvalCalc
from BuildPaths import modelProcPath, modelMedAdapPath, TasksJsonPath, BaseDir, modelPathMat
from paths import projFld

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get human1 metabolic model - generic and the one adapted for medium composition:
model = cobra.io.read_sbml_model(modelProcPath)
modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)

### Get essential tasks' list:
task_list = JSONTaskIO().read_task(TasksJsonPath)

# dict with tissue name and corresponding models:
# note: models of donnors/cell lines that although sucessfully reconstructed/gapfilled did not match between conditions were excluded
#      i.e. imagine study 1 has in CCs two sucessful models (donors A and B) and a another sucessful model for donor B in CSCs, so donor A is excluded
StdD = {'M2.2.13': ['M2.2.13.CCs_Hep3B', 'M2.2.13.CCs_Huh7', 'M2.2.13.CSCs_Hep3B', 'M2.2.13.CSCs_Huh7'],
        'M6.3.14': ['M6.3.14.CCs_1', 'M6.3.14.CSCs_1'],
        'M8.2.3': ['M8.2.3.CCs_1', 'M8.2.3.CCs_2', 'M8.2.3.CCs_3', 'M8.2.3.CSCs_1', 'M8.2.3.CSCs_2', 'M8.2.3.CSCs_3'],
        'M9.1.4':['M9.1.4.CCs_1', 'M9.1.4.CCs_2','M9.1.4.CCs_3', 'M9.1.4.CSCs_1_CXCR4+MET+CD44+', 'M9.1.4.CSCs_2_CXCR4+MET+CD44+', 'M9.1.4.CSCs_3_CXCR4+MET+CD44+'],
        'M13.2.11': ['M13.2.11.CCs_UT33', 'M13.2.11.CSCs_UT33'],
        'R1.1.1': ['R1.1.1.CCs_3', 'R1.1.1.CSCs_3'],
        'R4.10.2': ['R4.10.2.CCs_3', 'R4.10.2.CSCs_3'],
        'R5.8.2.TN': ['R5.8.2.TN.CCs_HCC1937', 'R5.8.2.TN.CCs_SUM149PT', 'R5.8.2.TN.CSCs_HCC1937', 'R5.8.2.TN.CSCs_SUM149PT'],
        'R11.1.4': ['R11.1.4.CCs_NCIH1703', 'R11.1.4.CSCs_NCIH1703'],
        'R12.5.4': ['R12.5.4.CCs_2', 'R12.5.4.CCs_5', 'R12.5.4.CSCs_2', 'R12.5.4.CSCs_5']
        }
# dict with tissue code vs tissue name:
TissueCode = {'1': 'prostate', '2': 'liver', '4': 'glioblastoma', '5': 'breast', '6':'ovary', '8': 'AML', '9': 'kidney', '11': 'lung', '12': 'pancreas', '13': 'head and neck'}

### Create plots with statistics for reconstructed models' composition (reactions/genes/metabolites):
AllMdDf = MdCmpPlot(BaseDir, StdD, modelMedAdap, TissueCode)
MdCmpHeatMap(model, modelPathMat, BaseDir, StdD, TissueCode, AllMdDf)

### Do pFBA simulation on reconstructed models:
pFBA(BaseDir, modelMedAdap)

### - Prepare for calculation of pvalue showing metabolic subsystem enrichment in different tissues:
prepPfba2PvalCalc(BaseDir, StdD, TissueCode)
