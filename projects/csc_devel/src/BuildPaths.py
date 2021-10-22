### Modules:
import os
from urllib.request import urlretrieve
from paths import BaseDir, mdDir, mdPthSBML, mdPathMat, mdPrcPath, mdPrcInfPath, mdMdAdpPath, mdInfPath, \
    MediumPathName, StudyNumberM, StudyNumberR, projFld

### Paths:
modelDir=os.path.join(BaseDir, mdDir) # directory where to retrieve model info
modelPathSBML, _ = urlretrieve(mdPthSBML)
modelPathMat = os.path.join(modelDir, mdPathMat) # path to generic .mat model
modelMedAdapPath = os.path.join(modelDir, mdMdAdpPath) # path to medium adapted model
modelInfoPath = os.path.join(modelDir, mdInfPath) # path to info on medium adapted model
modelProcPath = os.path.join(modelDir, mdPrcPath) # path to processed generic model
modelProcInfoPath = os.path.join(modelDir, mdPrcInfPath) # path to info on processed generic model
MediumPath = os.path.join(modelDir, MediumPathName) # path to file with chosen medium plus metabolites from ham's medium that allow generic model to grow
StudyDirM = os.path.join(BaseDir,'data/studies/' + str(StudyNumberM)) # directory with data of microarray study
GeneExpressionPathM = os.path.join(StudyDirM, 'NormData', 'NormGeneExp_Groups.tab')
StudyDirR = os.path.join(BaseDir,'data/studies/' + str(StudyNumberR)) # directory with data of rnaseq study
GeneExpressionPathR = os.path.join(StudyDirR, 'NormData', 'NormGeneExp_Groups.tab')
TasksExcelPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/56metTasks4iHuman.xlsx')
TasksJsonPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/56metTasks4iHuman.json')
LethalCellLinePath = "https://ndownloader.figshare.com/files/22543691" # link to deepmap achilles dataset v2
AllTasksJsonPath = os.path.join(os.path.expanduser(projFld), 'shared/task_sets/nl2019_tasks_r3d_compact.json')
ProcAllTasksJsonPath = os.path.join(BaseDir, 'support/MetbTasks/allTasks/processedTasks.json')
AssocReconHumanPath = os.path.join(BaseDir, 'support/MetbTasks/allTasks/humanGEMMetAssoc.JSON')
BiomassMetbExtID = 'temp001s' # id of biomass metabolite in external compartment










