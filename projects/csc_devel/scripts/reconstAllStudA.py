### Modules:
import sys
import cobra
import os
from reconstAllStudFunct import StudRcSc
from BuildPaths import modelProcPath, BaseDir
from paths import projFld, StudyNumberR, StudyNumberM

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get human1 generic metabolic model:
model = cobra.io.read_sbml_model(modelProcPath)

### get a list of models to reconst and corresponding dicts with (reaction:reaction scores) using the best threshold for rnaseq and microarray studies:
StudRcSc(BaseDir, model, StudyNumberR, StudyNumberM)