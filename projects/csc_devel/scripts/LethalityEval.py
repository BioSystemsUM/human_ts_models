'''
import sys
sys.path.extend(['/home/tbarata/CSCs/human_ts_models/projects/csc_devel/src', '/home/tbarata/cobamp/src', '/home/tbarata/troppo/src'])
'''
### Modules:
import os
from urllib.request import urlretrieve
from ModelInfo import preprocessiHuman
from tasksFunct import TestTasksUniversalModel, EssTaskRcFct
from LethalityEvalFunct import RcnstGapFill4ThrsEval, CorrExpSimulatedEssentialGenes

### Paths
#BaseDir = '/home/tbarata/CSCs/human_ts_models/projects/csc_devel'
BaseDir = '/home/tania/CSCs/human_ts_models/projects/csc_devel'
sampleM = 'M2.2.13.CCs_Huh7'  # select a cell line to test in metabolic simulation - of a microarray dataset
StudyNumberM = 'M2.1.13'
LethalityCellLineM = 'ACH-000480' # cell line (HUH7) from lethality experimental data. see sample_info.csv at 'support' directory
sampleR = 'R5.8.2.TN.CCs_HCC1937'  # select a cell line to test in metabolic simulation - of a rnaseq dataset
StudyNumberR = 'R5.8.2.TN'
LethalityCellLineR = 'ACH-000223' # cell line (HCC1937) from lethality experimental data. see sample_info.csv at 'support' directory
modelDir = os.path.join(BaseDir, 'data/MetModels/iHuman') # directory from which to retrieve model info
StudyDirM = os.path.join(BaseDir,'data/studies/' + str(StudyNumberM)) # directory with data of microarray study
GeneExpressionPathM = os.path.join(StudyDirM, 'NormData', 'NormGeneExp_Groups.tab')
StudyDirR = os.path.join(BaseDir,'data/studies/' + str(StudyNumberR)) # directory with data of rnaseq study
GeneExpressionPathR = os.path.join(StudyDirR, 'NormData', 'NormGeneExp_Groups.tab')
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
BiomassID = 'biomass_human' # model reaction id for biomass reaction
BiomassMetbExtID = 'temp001s' # id of biomass metabolite in external compartment
TasksExcelPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/56metTasks4iHuman.xlsx')
TasksJsonPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/56metTasks4iHuman.json')
LethalCellLinePath = os.path.join(BaseDir, 'data/Achilles_gene_effect_v2.csv')

### Pre-process iHuman metabolic model:
# remove blocked reactions and adapt boundaries for medium composition
model, modelMedAdap = preprocessiHuman(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, HamsPath, BaseMediumPath, MediumPath, BiomassID)

### Pre-process essential tasks and test them on universal model:
task_list = TestTasksUniversalModel(TasksExcelPath, TasksJsonPath, modelMedAdap)

### Identify reactions in generic model that if excluded cause one or more of the essential tasks to fail:
path = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/EssentialReactions4Tasks.tab')
EssTaskRcFct(modelMedAdap, task_list, path)

### Reconstruct models for thresholds+and/or rules combinations to test and gapfill for essential tasks and model growth:
protected = ['biomass_human', 'HMR_10023', 'HMR_10024']
RcnstGapFill4ThrsEval(BaseDir, sampleM, StudyNumberM, model, modelMedAdap, task_list, protected=protected, BiomassMetbExtID=BiomassMetbExtID)
RcnstGapFill4ThrsEval(BaseDir, sampleR, StudyNumberR, model, modelMedAdap, task_list, protected=protected, BiomassMetbExtID=BiomassMetbExtID)

### Determine correlation between experimental and simulated lethal genes for eah combination of threshold strategy + and/or rule:
# get dataframe with threshold comb and correlation scores
# best threshold comb is printed
finaldfM = CorrExpSimulatedEssentialGenes(LethalCellLinePath, modelMedAdap, BaseDir, StudyNumberM, LethalityCellLineM)
finaldfR = CorrExpSimulatedEssentialGenes(LethalCellLinePath, modelMedAdap, BaseDir, StudyNumberR, LethalityCellLineR)







