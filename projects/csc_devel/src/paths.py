BaseDir='projects/csc_devel' # basic directory
projFld='~/CSCs/human_ts_models' # path to git root project directory
mdPthSBML='https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/master/model/Human-GEM.xml' # path to generic SBML model
# following paths are relative to BaseDir:
mdDir='data/MetModels/iHuman' # directory where to retrieve model info
# following paths are relative to mdDir:
mdPathMat='Human-GEM.mat' # file name of generic .mat model
mdMdAdpPath='iHuman_MedAdap.xml' # file name of medium adapted model
mdInfPath='iHumanInfo.xlsx' # name of file with info on medium adapted model
MediumPathName='HAMs_iHuman.xlsx' # name of file with HAM's medium
mdPrcPath='iHuman_processed.xml' # file name of processed generic model
mdPrcInfPath='iHuman_processedInfo.xlsx' # name of file with info on processed generic model
# following are file names:
BiomassID='biomass_human' # model reaction id for biomass
BiomassMetbExtID='temp001s' # id of biomass metabolite in external compartment
sampleM='M2.2.13.CCs_Huh7' # name of sample of a microarray study to be used for best threshold simulation
StudyNumberM='M2.2.13' # name of microarray study which sample is to be used for best threshold simulation
LethalityCellLineM='ACH-000480' # deepmap cell line id corresponding to above mentioned microarray sample. for cell line ids go to https://depmap.org/portal/. cell line (HUH7) from lethality experimental data.
sampleR='R5.8.2.TN.CCs_HCC1937'  # select a cell line to test in metabolic simulation - of a rnaseq dataset
StudyNumberR='R5.8.2.TN'
LethalityCellLineR='ACH-000223' # cell line (HCC1937) from lethality experimental data. deepmap cell line id at https://depmap.org/portal/

