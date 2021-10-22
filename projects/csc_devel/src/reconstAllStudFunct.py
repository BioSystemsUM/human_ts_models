### Modules:
import os
import pickle
import pandas as pd
from LethalityEvalFunct import saveResult
from tasksFunct import EvalAllTasksAtSameTime


### Functions:
def StudRcSc(BaseDir, model, StudyNumberR, StudyNumberM):
    '''
    - get a list of models to reconst and corresponding dicts with (reaction:reaction scores) using the best thresholds for rnaseq and microarray
    :param BaseDir: basic directory
    :param model: generic model Not adapted for medium composition
    '''
    # open dict with threshold combinations as keys and dataframes as values, where columns are samples and values are reaction scores:
    pathAllStdAllThr = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct.pkl')
    infile = open(pathAllStdAllThr,'rb')
    AllStdAllThr = pickle.load(infile)
    infile.close()
    Datapath = os.path.join(BaseDir, 'data/studies')
    Dirs = [files for files in os.listdir(Datapath) if os.path.isdir(os.path.join(Datapath, files))]
    namesLst = list()
    RcScLst = list()
    for StudyNumber in Dirs:
        nmLst = list()
        rcLst = list()
        # identify best threshold from dataframe for each study/sample group:
        if StudyNumber.startswith('R'):
            thrsDf = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumberR, 'CorrScores.tab'), sep='\t', index_col=0)
        else:
            thrsDf = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumberM, 'CorrScores.tab'), sep='\t', index_col=0)
        besthres = thrsDf.loc[thrsDf['value'].idxmax()].drop(['value', 'biomass_cutoff', 'essential_cutoff'])
        thr = ('_'.join([besthres.Genes, besthres.thres, str(besthres.globalmin), str(besthres.globalmax), str(besthres.local)]), besthres.int_function)
        # add to list (for all studies) of model names and list of corresponding dictionaries with reaction scores:
        for smp in AllStdAllThr[thr]:
            if '.'.join(smp.split('.')[:-1]) == StudyNumber:
                nmLst.append(smp)
                rcLst.append(AllStdAllThr[thr][smp].fillna(0).to_dict()) # replace 'NaN' (corresponding to genes present in one microarray/rnaseq platform that are absent in another) by 0 (it's the score that is given to 'None' values bellow)
        # add scores 'None' to reactions that do not have gene-reaction-rule (cause those had been excluded from 'AllStudiesAllthresholdsDct.pkl' in function GetDictOfReacScoresDF of testGeneScrs.py):
        mdRcIds = {r.id for r in model.reactions}
        WithGeneReacRule = set(rcLst[0].keys())
        NoGeneReacRule = mdRcIds - WithGeneReacRule # get ids of reactions that do not have gene-reaction rule
        NoGeneReacRuleDct = {rid: None for rid in NoGeneReacRule}
        for d in rcLst:
            d.update(NoGeneReacRuleDct)
        namesLst = namesLst + nmLst
        RcScLst = RcScLst + rcLst
    if not os.path.exists(os.path.join(BaseDir, 'support/Reconst')):
        os.makedirs(os.path.join(BaseDir, 'support/Reconst'))
    saveResult(res=namesLst, path=os.path.join(BaseDir, 'support/Reconst/reconstMdNames.pkl'))  # save list of model names to test
    saveResult(res=RcScLst, path=os.path.join(BaseDir, 'support/Reconst/reconstScores.pkl'))  # save list of corresponding reaction scores to test
    pd.DataFrame([len(namesLst)]).to_csv(os.path.join(BaseDir, 'support/Reconst/NumbModels.tab'), sep='\t', header=False, index=False)
