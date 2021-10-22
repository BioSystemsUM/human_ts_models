### Modules:
import sys
import os
import cobra
from LethalityEvalFunct import RcnstThrsEval
from BuildPaths import modelProcPath, modelMedAdapPath, BaseDir, StudyNumberM, StudyNumberR
from paths import StudyNumberR, StudyNumberM, projFld

i = int(sys.argv[1]) # index to select one thresholds to test from list of all those to test
sample = sys.argv[2] # defines input file name/sample name, eg. M2.2.13.CCs_Huh7

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get iHuman metabolic model -generic and adapted for medium composition:
model = cobra.io.read_sbml_model(modelProcPath)
modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)

### Get list of reactions to keep during model reconstruction for thresholds+and/or rules combinations to test:
StudyNumber = '.'.join(sample.split('.')[:-1])
path2MdNames = os.path.join(BaseDir, 'support/TestGeneScores/models2Test', sample + '.pkl')  # list of threshold names to test
path2MdScores = os.path.join(BaseDir, 'support/TestGeneScores/rcScores2Test', sample + '.pkl')  # list of reaction scores dictionaries to test
pathDir = os.path.join(BaseDir, 'support/LethalityEval/Reac4Reconst/ReconstOut', StudyNumber)
RcnstThrsEval(BaseDir, StudyNumberR, StudyNumberM, model, i, path2MdNames, path2MdScores, pathDir, testThres=True, alg='both')

