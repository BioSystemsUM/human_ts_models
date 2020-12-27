### Modules:
'''
import sys
sys.path.extend(['/home/tbarata/CSCs/human_ts_models/projects/csc_devel/src', '/home/tbarata/cobamp/src', '/home/tbarata/troppo/src'])
'''
import os
from urllib.request import urlretrieve
from ModelInfo import preprocessiHuman
from tasksFunct import TestTasksUniversalModel
from reconstAllStudFunct import ReconstGapfill_AllMdsTasks, MdCmpPlot, MdCmpHeatMap, pFBA, ConvertTasksRecon2Human1

### Paths
#BaseDir = '/home/tbarata/CSCs/human_ts_models/projects/csc_devel'
BaseDir = '/home/tania/CSCs/human_ts_models/projects/csc_devel'
modelDir = os.path.join(BaseDir, 'data/MetModels/iHuman')
modelPathSBML, _ = urlretrieve('https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/master/model/Human-GEM.xml')
modelPathMat = os.path.join(modelDir, 'Human-GEM.mat')
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
TasksExcelPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/56metTasks4iHuman.xlsx')
TasksJsonPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/56metTasks4iHuman.json')
AllTasksJsonPath = os.path.join('/home/tania/CSCs/human_ts_models/shared/task_sets/nl2019_tasks_r3d_compact.json')
ProcAllTasksJsonPath = os.path.join(BaseDir, 'support/MetbTasks/allTasks/processedTasks.json')
AssocReconHumanPath = os.path.join(BaseDir, 'support/MetbTasks/allTasks/humanGEMMetAssoc.JSON')
BiomassMetbExtID = 'temp001s' # id of biomass metabolite in external compartment
StudyNumberR = 'R5.8.2.TN'
StudyNumberM = 'M2.1.13'

### Pre-process iHuman metabolic model:
# remove blocked reactions and adapt boundaries for medium composition
model, modelMedAdap = preprocessiHuman(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, HamsPath, BaseMediumPath, MediumPath, BiomassID)

### Pre-process essential tasks and test them on universal model:
task_list = TestTasksUniversalModel(TasksExcelPath, TasksJsonPath, modelMedAdap)

### Reconstruct models with best threshold and algorithm and gapfill for model growth and for essential tasks:
ReconstGapfill_AllMdsTasks(BaseDir, StudyNumberR, StudyNumberM, model, modelMedAdap, task_list, BiomassMetbExtID)

# dict with tissue name and corresponding models:
# note: models of donnors/cell lines that although sucessfully reconstructed/gapfilled did not match between conditions were excluded
#      i.e. imagine study 1 has in CCs two sucessful models (donors A and B) and a another sucessful model for donor B in CSCs, so donor A is excluded
StdD = {'M2.2.13': ['M2.2.13.CCs_Hep3B', 'M2.2.13.CCs_Huh7', 'M2.2.13.CSCs_Hep3B', 'M2.2.13.CSCs_Huh7'],
        'M8.2.3': ['M8.2.3.CCs_1', 'M8.2.3.CCs_2', 'M8.2.3.CCs_3', 'M8.2.3.CSCs_1', 'M8.2.3.CSCs_2', 'M8.2.3.CSCs_3'],
        'M9.1.4':['M9.1.4.CCs_2', 'M9.1.4.CSCs_2_CXCR4+MET+CD44+'],
        'M13.2.11': ['M13.2.11.CCs_UT16A', 'M13.2.11.CCs_UT24A', 'M13.2.11.CSCs_UT16A', 'M13.2.11.CSCs_UT24A'],
        'R1.1.1': ['R1.1.1.CCs_3', 'R1.1.1.CSCs_3'],
        'R5.8.2.TN': ['R5.8.2.TN.CCs_HCC1937', 'R5.8.2.TN.CSCs_HCC1937'],
        'R11.1.4': ['R11.1.4.CCs_ChaGoK1', 'R11.1.4.CSCs_ChaGoK1']
        }
# dict with tissue code vs tissue name:
TissueCode = {'1': 'prostate', '2': 'liver', '5': 'breast', '8': 'AML', '9': 'kidney', '11': 'lung', '13': 'head and neck'}

### Create plots with statistics for reconstructed models' composition (reactions/genes/metabolites):
AllMdDf = MdCmpPlot(BaseDir, StdD, modelMedAdap, TissueCode)
MdCmpHeatMap(model, modelPathMat, BaseDir, StdD, TissueCode, AllMdDf)

### Do pFBA simulation on reconstructed models:
pFBA(BaseDir, modelMedAdap, AllMdDf)

### Create plots for test of all tasks:
# convert tasks used in recon3d (consensus list) to human1:
ConvertTasksRecon2Human1(TJsonPath=AllTasksJsonPath, MetAssPath=AssocReconHumanPath, finalpath=ProcAllTasksJsonPath)

