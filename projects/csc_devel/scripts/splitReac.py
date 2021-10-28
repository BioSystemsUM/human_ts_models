### Modules:
import sys
import os
from ModelInfo import preprocessiHuman
from BuildPaths import BaseDir, modelInfoPath, MediumPath, modelProcPath, modelProcInfoPath, modelMedAdapPath, TasksExcelPath, TasksJsonPath, modelPathSBML
from paths import BiomassID, projFld
from LethalityEvalFunct import createLstRct2Test
from tasksFunct import TestTasksUniversalModel

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Pre-process Human1 metabolic model:
# adapt boundaries for medium composition
model, modelMedAdap = preprocessiHuman(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, MediumPath, BiomassID)

### Split model reaction ids (except boundaries/exchanges) into different files
# to be used for assessing reactions that are essential for essential tasks in generic model
createLstRct2Test(model, BaseDir)

### Pre-process essential tasks and test them on generic model:
TestTasksUniversalModel(TasksExcelPath, TasksJsonPath, modelMedAdap)

