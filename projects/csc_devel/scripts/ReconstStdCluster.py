### Modules:
import sys
import os
import cobra
from LethalityEvalFunct import RcnstThrsEval
from BuildPaths import modelProcPath, BaseDir
from paths import projFld, StudyNumberR, StudyNumberM

i = int(sys.argv[1]) # index to select one model to test from list of all models
projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get human1 generic metabolic model:
model = cobra.io.read_sbml_model(modelProcPath)

### Get list of reactions to keep during model reconstruction, for all studies:
path2MdNames = os.path.join(BaseDir, 'support/Reconst/reconstMdNames.pkl') # list of model names
path2MdScores = os.path.join(BaseDir, 'support/Reconst/reconstScores.pkl')  # list of corresponding reaction scores dictionaries
pathDir = os.path.join(BaseDir, 'support/Reconst/ReconstOut')
RcnstThrsEval(BaseDir, StudyNumberR, StudyNumberM, model, i, path2MdNames, path2MdScores, pathDir, testThres=False, alg='')

