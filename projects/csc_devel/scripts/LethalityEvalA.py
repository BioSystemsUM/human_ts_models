### Modules:
import sys
import os
import cobra
from tasksFunct import EssTaskRcFct
from LethalityEvalFunct import ThrsRcScores4aSample
from BuildPaths import BaseDir, modelProcPath
from paths import projFld

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get Human1 generic metabolic model:
model = cobra.io.read_sbml_model(modelProcPath)

### Save all reactions in generic model that if excluded cause one or more of the essential tasks to fail:
dir = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks')
path = os.path.join(dir, 'EssentialReactions4Tasks.tab')
EssTaskRcFct(path, BaseDir)

### get a list of thresholds to test and corresponding dicts with (reaction:reaction scores) for best thresholds of a rnaseq and of a microarray study:
ThrsRcScores4aSample(BaseDir, model)







