### Modules:
import sys
sys.path.extend(['/home/tbarata/CSCs/human_ts_models/projects/csc_devel/src', '/home/tbarata/cobamp/src', '/home/tbarata/troppo/src'])
import os
import pandas as pd
import cobra
from tasksFunct import EvalTasksReacKo
from BuildPaths import modelProcPath, modelMedAdapPath, TasksJsonPath, BaseDir
from paths import projFld
from troppo.tasks.task_io import JSONTaskIO

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get human1 metabolic generic model and generic model adapted for medium composition:
model = cobra.io.read_sbml_model(modelProcPath)
modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)

### Get task list:
task_list = JSONTaskIO().read_task(TasksJsonPath)

### Identify reactions essential for essential tasks in generic model:
# close model exchange/boundary reactions (in "human1" all boundaries are exchange reactions):
taskModel = modelMedAdap.copy()
for r in taskModel.boundary:
    r.knock_out()
# get group of reactions of generic model (except boundaries/exchange reactions):
fname = sys.argv[1] # name of file with a group of reactions (access bash variable)
frids = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/ModelRc2Test/' + fname)).iloc[:,0]
frids = frids.tolist()
# save reactions essential to at least one task into a list:
TasksReacKoLst = list()
for rid in frids:
    rid, bool = EvalTasksReacKo(rid, taskModel, task_list)
    if bool:
        TasksReacKoLst.append(rid)
i = fname[6:]
pd.Series(TasksReacKoLst).to_csv(os.path.join(BaseDir, 'support/LethalityEval/EssRcKo/EssRcKo' +str(i)), index=False)




