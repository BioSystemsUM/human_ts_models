### Modules:
import sys
import os
import cobra
from LethalityEvalFunct import CorrExpSimulatedEssentialGenes
from BuildPaths import BaseDir, modelMedAdapPath, LethalCellLinePath
from paths import projFld

projFld = os.path.expanduser(projFld)
os.chdir(projFld) # set working directory as root of project

### Get human1 metabolic generic model:
modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)

### Determine correlation between experimental and simulated lethal genes for best to test thresholds in the microarray and in the rnaseq study:
# dictionary mapping study to cell line id to test lethality genes correlation:
stdCLineD = {'M2.2.13':'ACH-000480', 'R5.8.2.TN':'ACH-000223'}
# get dataframe with threshold comb and correlation scores
# best threshold comb is printed:
for StudyNumber, LethalityCellLine in stdCLineD.items():
    CorrExpSimulatedEssentialGenes(LethalCellLinePath, modelMedAdap, BaseDir, StudyNumber = StudyNumber, LethalityCellLine = LethalityCellLine)

