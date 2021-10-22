### Modules:
import sys
import os
import cobra
from LethalityEvalFunct import RcnstGapFill4ThrsEval
from BuildPaths import modelMedAdapPath, TasksJsonPath, BaseDir
from troppo.tasks.task_io import JSONTaskIO
from paths import projFld, BiomassMetbExtID

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project
fname=sys.argv[1] # name of file with reconst. alg. output - list of reactions to keep

### Get human1 metabolic model adapted for medium composition:
modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)

### Get task list:
task_list = JSONTaskIO().read_task(TasksJsonPath)

### Reconstruct models for sample studies and gapfill for essential tasks and model growth:
protected = ['biomass_human', 'HMR_10023', 'HMR_10024'] # biomass reactions to protect
Finalpath = os.path.join(BaseDir, 'results/Reconst')
if not os.path.exists(Finalpath):
    os.makedirs(Finalpath)
Intpath = os.path.join(BaseDir, 'support/Reconst/ReconstOut', fname)
RcnstGapFill4ThrsEval(BaseDir, BiomassMetbExtID, StudyNumber='', modelMedAdap=modelMedAdap, task_list=task_list, protected=protected, Intpath=Intpath, Finalpath=Finalpath, testThres=False, fname=fname)


