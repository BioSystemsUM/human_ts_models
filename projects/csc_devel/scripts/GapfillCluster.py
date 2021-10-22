### Modules:
import sys
import os
import cobra
from LethalityEvalFunct import RcnstGapFill4ThrsEval
from BuildPaths import BaseDir, modelMedAdapPath, TasksJsonPath
from paths import projFld, BiomassMetbExtID
from troppo.tasks.task_io import JSONTaskIO

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

StudyNumber=sys.argv[1] # study number
fname=sys.argv[2] # name of file with reconst. alg. output - list of reactions to keep

### Get human1 metabolic model generic model adapted for medium composition:
modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)

### Get task list:
task_list = JSONTaskIO().read_task(TasksJsonPath)

### Reconstruct models for thresholds+and/or rules combinations to test and gapfill for essential tasks and model growth:
protected = ['biomass_human', 'HMR_10023', 'HMR_10024'] # biomass reactions to protect
Finalpath = os.path.join(BaseDir, 'support/LethalityEval/GapFillRes')
if not os.path.exists(Finalpath):
    os.makedirs(Finalpath)
Intpath = os.path.join(BaseDir, 'support/LethalityEval/Reac4Reconst/ReconstOut', StudyNumber, fname)
RcnstGapFill4ThrsEval(BaseDir, BiomassMetbExtID, StudyNumber=StudyNumber, modelMedAdap=modelMedAdap, task_list=task_list, protected=protected, Intpath=Intpath, Finalpath=Finalpath, testThres=True, fname=fname)