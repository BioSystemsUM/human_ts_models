### Modules:
import os
import numpy as np
import pickle
import ast
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import cobra
import json
import copy
from functools import reduce
from LethalityEvalFunct import Reconst_taskEval, ReconstructionReac_ReacScores
from sklearn.preprocessing import MinMaxScaler
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as dst
from cobra.flux_analysis import flux_variability_analysis, pfba
from tasksFunct import EvalAllTasksAtSameTime
from troppo.tasks.task_io import JSONTaskIO

### Functions:
def ReconstGapfill_AllMdsTasks(BaseDir, StudyNumberR, StudyNumberM, model, modelMedAdap, task_list, BiomassMetbExtID):
    '''
    - Reconstruct models with best threshold and algorithm and gapfill for model growth and for essential tasks
    :param BaseDir: basic directory
    :param StudyNumberR: study number of rna-seq study for which best threshold + algorithm was assessed
    :param StudyNumberM:study number of microarray study for which best threshold + algorithm was determined
    :param model: generic model Not adapted for medium composition
    :param modelMedAdap: generic model adapted for medium composition
    :param task_list: list of tasks to test and to gapfill for
    :param BiomassMetbExtID: id of biomass metabolite in external compartment
    '''
    # open dict with threshold combinations as keys and dataframes as values, where columns are samples and values are reaction scores:
    pathAllStdAllThr = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct.pkl')
    infile = open(pathAllStdAllThr,'rb')
    AllStdAllThr = pickle.load(infile)
    infile.close()
    # identify best threshold from dataframe:
    thrsDfR = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumberR, 'CorrScores.tab'), sep='\t', index_col=0)
    besthresR = thrsDfR.loc[thrsDfR['value'].idxmax()].drop(['value', 'biomass_cutoff', 'essential_cutoff'])
    thrsDfM = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumberM, 'CorrScores.tab'), sep='\t', index_col=0)
    besthresM = thrsDfM.loc[thrsDfM['value'].idxmax()].drop(['value', 'biomass_cutoff', 'essential_cutoff'])
    # as the best threshold for microarray studies is also the best for rnaseq studies (see besthresR == besthresM) we will no seprarate microarray from rnaseq studies:
    thr = ('_'.join([besthresR.Genes, besthresR.thres, str(besthresR.globalmin), str(besthresR.globalmax), str(besthresR.local)]), besthresR.int_function)
    alg = besthresR.algorithm
    # get list of model names and list of corresponding dictionaries with reaction scores:
    namesLst = list()
    RcScLst = list()
    for c in AllStdAllThr[thr]:
        namesLst.append(c)
        RcScLst.append(AllStdAllThr[thr][c].to_dict())
    # get reactions to include/exclude from reconstructed models according to best threshold and algorithm:
    path = os.path.join(BaseDir, 'results/Reconst', 'reconstReac.csv')
    ReconstructionReac_ReacScores(model, namesLst, RcScLst, path, alg=alg, testThres=False)
    # gapfill reconstructed models with "EFM" method and gapfill for tasks using cobrapy gapfill:
    protected = ['biomass_human', 'HMR_10023', 'HMR_10024']
    Finalpath = os.path.join(BaseDir, 'results/Reconst')
    BoolDict = pd.read_csv(path, header=0, index_col=[0, 1]).T.to_dict()  # dict where key is ('algorithm', 'model') and value another dict with (reactID: True/False)
    Reconst_taskEval(modelMedAdap, task_list, BiomassMetbExtID, BaseDir, protected, Finalpath, BoolDict, StudyNumb='')

def GetMetbsGenesFromReact(model, ReacBoolDf):
    '''
    - returns dataframes with True/False for reaction/genes/metabolies in models/cell types (indicating if reaction/gene/metabolite is present or not)
    :param model: univeral model
    :param ReacBoolDf: dataframe with True/False for reactions in models/cell types (indicating if reaction is present or not)
    :return ReacBoolDf: input ReacBoolDf with extra row wih total number of reactions
    :return GMdDf: dataframe with True/False for genes in models/cell types (indicating if gene is present or not)
    :return MMdDf: dataframe with True/False for metabolites in models/cell types (indicating if metabolite is present or not)
    '''
    AllGenes = [g.id for g in model.genes] # list of genes in universal model
    AllMtbs = [m.id for m in model.metabolites] # list of metabolites in universal model
    genesD = dict()
    metbsD = dict()
    BoolDf = ReacBoolDf.copy()
    for rcModelName in BoolDf:  # for each model/cell type
        mdActGenesLst = list()
        mdActMetbsLst = list()
        rcModelBool = BoolDf[rcModelName]  # get reaction True/False of a model/cell type
        for rcId, rcVl in zip(rcModelBool.index, rcModelBool):  # for each reaction
            if rcVl:
                Arc = model.reactions.get_by_id(rcId)
                for g in Arc.genes:
                    mdActGenesLst.append(g.id)
                for m in Arc.metabolites:
                    mdActMetbsLst.append(m.id)
        mdActGenesLst = np.unique(mdActGenesLst) # all active genes in that model/cell type
        mdActMetbsLst = np.unique(mdActMetbsLst) # all active metabolites in that model/cell type
        mdGenesD = {gid: True if gid in mdActGenesLst else False for gid in AllGenes} # dict with gene id vs 'True' when gene is in that model and 'False' when it is not
        mdMtbsD = {mid: True if mid in mdActMetbsLst else False for mid in AllMtbs} # dict with metabolite id vs 'True' when metabolite is in that model and 'False' when it is not
        genesD[rcModelName] = mdGenesD # dict with model id vs inner dict with gene id and True/False
        metbsD[rcModelName] = mdMtbsD # dict with model id vs inner dict with metabolite id and True/False
    GMdDf = pd.DataFrame.from_dict(genesD) # dataframe with genes on rows and models on columns where values are True/False
    MMdDf = pd.DataFrame.from_dict(metbsD) # dataframe with metabolites on rows and models on columns where values are True/False
    # add row to dataframes with sum of all active reactions/metabolites/genes of each model/cell type:
    BoolDf.loc['Reactions'] = BoolDf.sum() # 1 means 'True', 0 means 'False'
    GMdDf.loc['Genes'] = GMdDf.sum()
    MMdDf.loc['Metabolites'] = MMdDf.sum()
    return BoolDf, GMdDf, MMdDf

def medianMdSameStd(MdTotalDf, StdD, TissueCode, axis, datype):
    '''
    - given a series/dataframe with values for each model, obtain same dict with with each tissue and cell type as key and corresponding median value number
    :param MdTotalDf: series/dataframe with values for different models
    :param StdD: dict with study number and corresponding models
    :param TissueCode: dict with tissue id tissue name correspondence
    :param axis: axis across which aplly medium. axis 0 is across rows and axis 1 across columns. for series use axis 0
    :param datype: string indicating MdTotalDf data type. either 'series' or 'dataframe'
    '''
    AvgD = dict()
    for std, mds in StdD.items(): # for each study and corresponding models
        CCsmp = [mdN for mdN in mds if 'CCs' in mdN] # get CCs model names for that study
        CSCsmp = [mdN for mdN in mds if 'CSCs' in mdN] # get CSCs model names for that study
        ts = TissueCode[std.split('.')[0][1:]]
        if datype == 'series':
            AvgD['.'.join([ts, 'CCs'])] = [MdTotalDf[CCsmp].median(axis=axis)]
            AvgD['.'.join([ts, 'CSCs'])] = [MdTotalDf[CSCsmp].median(axis=axis)]
        elif datype == 'dataframe':
            AvgD['.'.join([ts, 'CCs'])] = MdTotalDf[CCsmp].median(axis=axis)
            AvgD['.'.join([ts, 'CSCs'])] = MdTotalDf[CSCsmp].median(axis=axis)
    return AvgD

def MdCmpPlot(BaseDir, StdD, modelMedAdap, TissueCode):
    '''
    - create boxplot and barplots with statistics for reconstructed model composition (reactions/genes/metabolites)
    :param BaseDir: basic directory
    :param StdD: dict with tissue name and corresponding model
    :param modelMedAdap: generic metabolic model adapted for medium composition
    :param TissueCode: dictionary with number identifying a tissue vs corresponding tissue name
    :return AllMdDf: dataframe with reactions on rows and models on columns with 'True' values when reaction is present in model
    '''
    # get dataframe with model names in columns, reactions as rows and values (LB, UB):
    StdPath = os.path.join(BaseDir, 'results/Reconst2/RecModelBounds.tab')
    modelBds = pd.read_csv(StdPath, sep='\t', index_col=0)
    modelBds.columns = list(map(lambda x: '_'.join(x.split('_')[1:]),modelBds.columns))
    # select just matched donors that had models sucessfully reconstructed and gapfilled:
    modelBds = modelBds[reduce(lambda x,y: x+y, StdD.values())]
    # when bounds are different from (0,0) replace by 'True' otherwise 'False':
    AllMdDf = modelBds.applymap(lambda x: True if ast.literal_eval(x) != (0,0) else False)
    # get dataframes with reactions/genes/metabolites on rows and models on columns where values are 1 when reaction/gene/metabolite is active,
    # the last row has total number of active reactions/genes/metabolites:
    RctMdDf, GMdDf, MMdDf = GetMetbsGenesFromReact(model=modelMedAdap, ReacBoolDf=AllMdDf)
    # get total number of active reactions/genes/metabolites per model:
    MdTotalRctDf = pd.DataFrame(RctMdDf.loc['Reactions'])
    MdTotalRctDf.columns = ['Number']
    MdTotalGDf = pd.DataFrame(GMdDf.loc['Genes'])
    MdTotalGDf.columns = ['Number']
    MdTotalMDf = pd.DataFrame(MMdDf.loc['Metabolites'])
    MdTotalMDf.columns = ['Number']
    # concat reactions/genes/metabolites info on a dataframe:
    MdTotalRctDf['Property'] = pd.Series(['Reactions'] * MdTotalRctDf.shape[0]).values
    MdTotalGDf['Property'] = pd.Series(['Genes'] * MdTotalGDf.shape[0]).values
    MdTotalMDf['Property'] = pd.Series(['Metabolites'] * MdTotalMDf.shape[0]).values
    MedDf = pd.concat([MdTotalRctDf, MdTotalGDf, MdTotalMDf])
    MedDf.reset_index(inplace=True)
    MedDf.rename(columns={'index':'Donor'}, inplace=True)
    MedDf['Cell Type'] = MedDf['Donor'].apply(lambda x: x.split('.')[len(x.split('.'))-1].split('_')[0]) # add column indicating if is CCs or CSC
    MedDf['Tissue'] = MedDf['Donor'].apply(lambda x : TissueCode[x.split('.')[0][1:]]) # add column indicating study name
    MedDf['Number'] = MedDf['Number'].astype('float')
    JRct = MedDf[MedDf['Property']=='Reactions']
    JRct = JRct.rename(columns={'Number':'Number of Reactions'})
    # boxplot split by tissue:
    plt.figure(figsize=(12, 7))
    Medbplt = sns.boxplot(x='Tissue', y='Number of Reactions', data=JRct, palette='coolwarm', hue='Cell Type')
    handles, labels = Medbplt.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2])
    fig = Medbplt.get_figure()
    plt.tight_layout()
    fig.savefig(os.path.join(BaseDir, 'results/Pathways','Rc_boxplotissue.png'))
    plt.close(fig)
    # barplot split by tissue:
    plt.figure(figsize=(12,7))
    nplt = sns.barplot(x='Tissue', y='Number of Reactions', data=JRct, palette='coolwarm', hue='Cell Type')
    nplt.set_xlabel('Study',fontsize=15)
    nplt.set_ylabel('Number of Reactions',fontsize=15)
    nplt.tick_params(labelsize=10)
    nplt.set_xticklabels(nplt.get_xticklabels(), rotation=90)
    plt.tight_layout()
    nplt.legend(loc=1)
    fig = nplt.get_figure()
    fig.savefig(os.path.join(BaseDir, 'results/Pathways','Rc_barplotissue.png'))
    plt.close(fig)
    # boxplot:
    plt.figure(figsize=(12,7))
    Medbxplt = sns.boxplot(x='Property', y='Number', data=MedDf, palette='coolwarm', hue='Cell Type')
    handles, labels = Medbxplt.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2])
    fig = Medbxplt.get_figure()
    fig.savefig(os.path.join(BaseDir, 'results/Pathways','RcGeneMetb_boxplot.png'))
    plt.close(fig)
    return AllMdDf

def GetRcSbs(model, modelPathMat):
    '''
    - create dict with reaction id vs reaction subsystem
    :param model: generic metabolic modelin .sbml format
    :param modelPathMat: path to generic model in .mat format, which contains model reactions' subsystems
    '''
    #create dict with reaction id vs reaction subsytem
    modelMat = cobra.io.load_matlab_model(modelPathMat) # retrieve reaction subsystem from .mat model
    RcSbsD = dict()
    for r in model.reactions:
        rSubs = modelMat.reactions.get_by_id(r.id).subsystem
        r_sub = re.sub(" ","_",re.sub(" / ", "/", re.sub("\[array\(\[|\], dtype=<.*\)\]|\],\n.*", "", re.sub("'", "", rSubs))))
        RcSbsD[r.id] = [r_sub]
    return RcSbsD

def applyclustermap(data, filePath, csize, rsize, cclust, rclust, fcl, btSp, fts):
    '''
    - does a clean clustermap graph
    :param data: data for heatmap
    :param filePath: path to file where to save graph
    :param csize: width of a column
    :param rsize: width of a row
    :param cclust: bool indicating whether to cluster (show dendrogram) columns or not
    :param rclust: bool indicating whether to cluster rows or not
    :param fcl: factor controlling overall font scale
    :param btSp: factor controlling bottom white space for column names
    :param fts: factor controlling row font size
    :return: saves the graph
    '''
    nc = data.shape[1]
    nr = data.shape[0]
    Fcl = fcl / nc
    sns.set(font_scale=fcl)
    res = sns.clustermap(data, col_cluster=cclust, row_cluster=rclust, cmap='coolwarm', linecolor='black', xticklabels=True, yticklabels=True, linewidths=0.15, figsize=(csize * nc, rsize * nr), dendrogram_ratio=(0.2, 0.05), cbar=True, cbar_kws={"shrink": .82})
    res.fig.subplots_adjust(bottom=btSp * Fcl)
    res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_ymajorticklabels(), fontsize= fts / nr)
    res.savefig(filePath, dpi=800)
    plt.close('all')

def MdCmpHeatMap(model, modelPathMat, BaseDir, StdD, TissueCode, AllMdDf):
    '''
    - create heatmap with tissue vs subsystem with values:
      * 1, when on the top 10% subsytems of list of CSC-CC difference
      * -1, when on the bottom 10% subsytems of list of CSC-CC difference
      * 0, otherwise
    - count number of tissues that subsystems on top/bottom of list have and save in dataframe
    :param model: universal model
    :param modelPathMat: path to model in format .mat, which contains model reactions' subsystems
    :param BaseDir: base directory
    :param StdD: dictionary with study names as keys and corresponding model names as values
    :param TissueCode: dictionary with tissues ids and corresponding tissues' names
    :param AllMdDf: dataframe with booleans indicating whether a reaction is 'active' (True) or not (False) in a model. columns are models and rows are reactions
    '''
    # create dict with reaction id vs reaction subsystem:
    RcSbsD = GetRcSbs(model, modelPathMat)
    # add subsystem info to dataframe with 'True' values when reaction is active in model:
    RcSbsF = pd.DataFrame.from_dict(RcSbsD, orient='index', columns=['Subsystem'])
    mdDf = RcSbsF.join(AllMdDf, how='outer')
    RcSbsF.to_csv(os.path.join(BaseDir, 'results/subsytemInfo.tab'), sep='\t') # save subsystem info in a dataframe
    # percentage of reactions in each subsytem that are active in each model:
    AllMdDfSbs = mdDf.groupby('Subsystem').sum()
    NbRcInSbs = mdDf.pivot_table(index=['Subsystem'], aggfunc='size') # number of reactions in each subsytem in universal model
    PercSbsRct = AllMdDfSbs.apply(lambda x: x/NbRcInSbs*100)
    # get median percentage of each subsystem with active reactions, per tissue and cell type (CCs or CSCs):
    MedD = medianMdSameStd(PercSbsRct, StdD, TissueCode, axis=1, datype='dataframe')
    MedDf = pd.DataFrame.from_dict(MedD)
    # get difference beween CSCs and CCs in each tissue:
    MedDdiff = {ts: MedDf[ts + '.CSCs'] -  MedDf[ts + '.CCs'] for ts in TissueCode.values()}
    MedDfDiff = pd.DataFrame.from_dict(MedDdiff)
    # exclude subsystems with low number of reactions (3 or less):
    CSC_CC_df = MedDfDiff[NbRcInSbs >= 4]
    NbRcInSbs[NbRcInSbs < 4].to_csv(os.path.join(BaseDir, 'results/Pathways', 'SubsLessThan4Reac'), sep='\t', header=False) # save info on excluded subsytems with low number of reactions (less than 4)
    # exclude 'exchange/demand reactions':
    # (number of these reactions is dif. depending on .sbml and .mat model, besides in this analysis they do not change between reconstructed models)
    CSC_CC_df = CSC_CC_df.drop(index='Exchange/demand_reactions')
    # exclude subsystems where no tissue had active reactions (all row is 0):
    CSC_CC_df = CSC_CC_df.loc[(CSC_CC_df != 0).any(axis=1)]
    # exclude other weird subsytems:
    CSC_CC_df = CSC_CC_df.drop(index='Isolated')
    CSC_CC_df = CSC_CC_df.drop(index='Miscellaneous')
    # get top and bottom pathways in each study and then merge info:
    perc = int(CSC_CC_df.shape[0] * 0.10) # number of subsytems in top and bottom 10%
    tCSCL = list()
    tCCL = list()
    for ts in CSC_CC_df: # for each tissue get order subsytems by those with more existent reactions in CSCs to those with more in CCs
        df = CSC_CC_df[ts].sort_values(ascending=False)
        topCSC = pd.DataFrame(df.iloc[:perc+1])
        topCCs = pd.DataFrame(df.iloc[-perc:])
        tCSCL.append(topCSC)
        tCCL.append(topCCs)
    tCSCDf = reduce(lambda x,y: x.join(y, how='outer'), tCSCL)
    tCCDf = reduce(lambda x, y: x.join(y, how='outer'), tCCL)
    # get number of tissues where a subsystem is on top/bottom of ordered list (on top of those with more existing reactions in CSCs than CCs and vice-versa):
    tCSCDfCp = tCSCDf.copy()
    tCCDfCp = tCCDf.copy()
    tCSCDf[tCSCDf.notna()] = 1 # when subsytem has a value in top/bottom, counts as 1, otherwise (if nan) counts as 0
    tCSCDf = tCSCDf.fillna(0)
    tCSCDf.sum(axis=1).sort_values(ascending=False).to_csv(os.path.join(BaseDir, 'results/Pathways/MoreCSCRank.tab'), sep='\t', header = True)
    tCCDf[tCCDf.notna()] = 1
    tCCDf = tCCDf.fillna(0)
    tCCDf.sum(axis=1).sort_values(ascending=False).to_csv(os.path.join(BaseDir, 'results/Pathways/MoreCCRank.tab'), sep='\t', header=True)
    # heatmap with existence of subsytems on top/bottom of ordered CSC-CC difference list:
    tCSCDfCp[tCSCDfCp.notna()] = 1 # using the copied dataframes above, if subs. is "on top of CSC-CC list" counts as 1
    tCCDfCp[tCCDfCp.notna()] = - 1 # using the copied dataframes above, if subs. is "on bottom of CSC-CC list" counts as -1
    joined = tCSCDfCp.join(tCCDfCp, how='outer', lsuffix = '.CSC', rsuffix = '.CC') # allows for index of both dataframes to be the same, by merging them
    MCSC = joined[['.'.join([c, 'CSC']) for c in tCCDfCp.columns]]  # then we split/subset dataframes again
    MCSC.columns = [c.split('.')[0] for c in MCSC.columns]
    MCC = joined[['.'.join([c, 'CC']) for c in tCCDfCp.columns]] # split/subset dataframes again
    MCC.columns = [c.split('.')[0] for c in MCC.columns]
    new = MCSC.copy()
    new[(MCSC == 1) & (MCC.isna())] = 1 # if subsytem is "on top of CSC-CC list" and not in bottom counts as 1
    new[(MCSC.isna()) & (MCC == -1)] = -1 # if subsyem is "on bottom of CSC-CC list" and not in top counts as -1
    new = new.fillna(0) # all NAs (neither on top nor bottom of CSC-CC list) are replaced by 0
    path = os.path.join(BaseDir, 'results/Pathways', 'clusterMoreDiff.png')
    applyclustermap(data=new, filePath=path, csize=0.7, rsize=0.2, cclust=False, rclust=False, fcl=2, btSp=0.2, fts=800)

def pFBA(BaseDir, modelMedAdap, AllMdDf):
    '''
    - do pFBA in reconstructed models and saves info on dataframe
    :param BaseDir: basic directory
    :param modelMedAdap: generic metabolic modle adapted for medium composition
    :param AllMdDf: dataframe with booleans indicating whether a reaction is 'active' (True) or not (False) in a model. columns are models and rows are reactions
    '''
    excId = [exc.id for exc in modelMedAdap.exchanges] # get exchange reactions ids
    PfbaFluxD = dict()
    for mdN in AllMdDf: # for each model
        rcMd = modelMedAdap.copy()
        mdSeries = AllMdDf[mdN]
        for rcid in mdSeries.index: # for each reaction
            rcbool = mdSeries[rcid]
            if (not rcbool) and (rcid not in excId):
                rcMd.reactions.get_by_id(rcid).bounds = (0,0)
        opt = rcMd.optimize()
        if opt.status != 'optimal' or opt.objective_value < 1E-9:  # if model is 'infeasible' or doesn't grow, print warning and continue to next model
            print(mdN + ' infeasible or does not grow')
            continue
        # do pFBA on the model:
        pfbasol = pfba(rcMd)
        pfbaflx = pfbasol.fluxes
        PfbaFluxD[mdN] = pfbaflx # save model pFBA fluxes to dictionary with model name as key
    pd.DataFrame.from_dict(PfbaFluxD).to_csv(os.path.join(BaseDir, 'results/pFBA.tab'), sep='\t')

def prepPfba2PvalCalc (BaseDir, StdD, TissueCode):
    '''
    - apply module to pFBA fluxes, get median for each tissue and cell type and get difference between CSCs and CCs in each tisue; save results
    :param BaseDir: basic directory
    :param StdD: dictionary with study names as keys and corresponding model names as values
    :param TissueCode: dictionary with tissues ids and corresponding tissues' names
    '''
    # get pfba flux results:
    pfba_flux = pd.read_csv(os.path.join(BaseDir, 'results/pFBA.tab'), sep='\t', index_col=[0])
    # apply absolute value to each flux:
    AbsPfba = pfba_flux.applymap(abs)
    # cal. median for each tissue and cell type:
    tsAbsPfba = medianMdSameStd(AbsPfba, StdD, TissueCode, axis=1, datype='dataframe')
    tspfba = pd.DataFrame.from_dict(tsAbsPfba)
    # get difference beween CSCs and CCs in each tissue:
    pfbadifD = {ts: tspfba[ts + '.CSCs'] -  tspfba[ts + '.CCs'] for ts in TissueCode.values()}
    pfbadif = pd.DataFrame.from_dict(pfbadifD)
    pfbadif.to_csv(os.path.join(BaseDir, 'results/pfbaOrdAbsVal.tab'),sep='\t')

def Eval256TasksReconstMds (modelMedAdap, AllMdDf, AllTaskLst):
    # get reconstructed models:
    excId = [exc.id for exc in modelMedAdap.exchanges] # get exchange reactions ids
    for mdN in AllMdDf: # for each model
        rcMd = modelMedAdap.copy()
        mdSeries = AllMdDf[mdN]
        for rcid in mdSeries.index: # for each reaction
            rcbool = mdSeries[rcid]
            if (not rcbool) and (rcid not in excId):
                rcMd.reactions.get_by_id(rcid).bounds = (0,0)
        # eval tasks:
        for k in rcMd.boundary:  # close reconstructed model boundaries to test tasks
            k.knock_out()
        tasksName, batch_res_tasks = EvalAllTasksAtSameTime(rcMd, AllTaskLst)
        # save tasks eval result:
        TaskResD = dict()
        for Tname, res in zip(tasksName, batch_res_tasks):
            print(Tname, res)
            TaskResD[Tname] = res[0]


metbTasks = pd.read_excel(AllTasksExcelPath, sheet_name='TASKS', usecols=['ID', 'DESCRIPTION', 'SHOULD FAIL', 'IN', 'IN LB', 'IN UB', 'OUT', 'OUT LB', 'OUT UB', 'EQU', 'EQU LB', 'EQU UB'])
metbTasks['SHOULD FAIL'].fillna(0, inplace=True)
metbTasks['SHOULD FAIL'] = metbTasks['SHOULD FAIL'].apply(bool)
metbTasks['IN LB'].fillna(0, inplace=True) # when no value means LB = 0
metbTasks['IN UB'].fillna(1000, inplace=True) # when no value means UB = 1000
metbTasks['OUT LB'].fillna(0, inplace=True) # when no value means LB = 0
metbTasks['OUT UB'].fillna(1000, inplace=True) # when no value means UB = 1000
# when '=>' in 'EQU' and 'EQU LB' has no value means 'EQU LB' = 0:
metbTasks['EQU LB'][np.isnan(metbTasks['EQU LB']) & metbTasks['EQU'].apply(lambda x: '<=>' not in x if isinstance(x, str) else False)] = 0
# when '<=>' in 'EQU' and 'EQU LB' has no value means 'EQU LB' = -1000:
metbTasks['EQU LB'][np.isnan(metbTasks['EQU LB']) & metbTasks['EQU'].apply(lambda x: '<=>' in x if isinstance(x, str) else False)] = -1000
# when 'EQU' has a equation (is a string) and 'EQU UB' has no value means 'EQU UB' = 1000:
metbTasks['EQU UB'][metbTasks['EQU'].apply(lambda x: True if isinstance(x, str) else False) & metbTasks['EQU UB'].apply(lambda x: np.isnan(x))] = 1000
metbTasks['EQU'] = metbTasks['EQU'].apply(str) # for the column to avoid 'ffil' operation bellow
metbTasks['EQU LB'] = metbTasks['EQU LB'].apply(str) # for the column to avoid 'ffil' operation bellow
metbTasks['EQU UB'] = metbTasks['EQU UB'].apply(str) # for the column to avoid 'ffil' operation bellow
metbTasks = metbTasks.fillna(method='ffill') # for metabolites of same task in different rows to have same 'ID' and 'DESCRIPTION'
#metbTasks['ID'] = metbTasks['ID'].apply(lambda x: str(int(x)))
mt = metbTasks.copy()
mt['IN'] = metbTasks['IN'].apply(lambda x: x.split(';')) # split metabolites 'IN' into a list
mt['OUT'] = metbTasks['OUT'].apply(lambda x: x.split(';')) # split metabolites 'OUT' into a list
mt['IN LB'] = mt['IN LB'].apply(lambda x: [x]) * mt['IN'].apply(len) # repeat LB values for each 'IN' metabolite
mt['IN UB'] = mt['IN UB'].apply(lambda x: [x]) * mt['IN'].apply(len) # repeat UB values for each 'IN' metabolite
mt['OUT LB'] = mt['OUT LB'].apply(lambda x: [x]) * mt['OUT'].apply(len) # repeat LB values for each 'OUT' metabolite
mt['OUT UB'] = mt['OUT UB'].apply(lambda x: [x]) * mt['OUT'].apply(len) # repeat UB values for each 'OUT' metabolite


df = pd.DataFrame(mt.groupby('ID')['IN'].apply(lambda x: set([el for o in list(x) for el in o])))
df.reset_index(inplace=True)
df = df[['ID', 'IN']]
def removeDupBounds(clNm):
    mNm = clNm.split(' ')[0]
    a = mt.groupby('ID')[clNm].apply(lambda x: [el for o in list(x) for el in o])
    return pd.Series([v[:s] for v, s in zip(a,df[mNm].apply(len))])
df['IN LB'] = removeDupBounds(clNm = 'IN LB')
df['IN UB'] = removeDupBounds(clNm = 'IN UB')
df['OUT'] = mt.groupby('ID')['OUT'].apply(lambda x: set([el for o in list(x) for el in o])).reset_index(drop=True)
df['OUT LB'] = removeDupBounds(clNm = 'OUT LB')
df['OUT UB'] = removeDupBounds(clNm = 'OUT UB')

from urllib.request import urlretrieve
from troppo.tasks.task_io import ExcelTaskIO
URL = 'https://github.com/SysBioChalmers/Human-GEM/raw/master/data/metabolicTasks/metabolicTasks_Essential.xlsx'

path, _ = urlretrieve(URL)
path = '/home/tania/metabolicTasks_Essential.xlsx'
task_list = ExcelTaskIO().read_task(path)


def ConvertTasksRecon2Human1(TJsonPath, MetAssPath, finalpath):
    '''
    - convert tasks with metabolite ids of recon3d to tasks with metabolite ids of human1
    :param TJsonPath: path to .json file with tasks from recon3d (consensus) list
    :param MetAssPath: path to .json file with association between recon3d and human1 metabolites' ids
    :param finalpath: path to json file that contains the result, a list of tasks that can be tested on human1 model
    '''
    # get tasks with metabolite ids in recon3d format:
    tasks = JSONTaskIO().read_task(TJsonPath)
    # get datframe with association between metabolites ids in recon3d and ihuman and compartments:
    with open(MetAssPath) as json_file:
        metbass = json.load(json_file)
    df = pd.concat([pd.Series(metbass['mets']), pd.Series(metbass['metRecon3DID'])], axis=1)
    df.columns = ['Human1', 'Recon3D']
    df['Compartment'] = df['Human1'].apply(lambda x: x[-1])
    # exclude tasks with metabolites that do not exist in human1:
    t2Remove = list()
    for t in tasks:
        INB = [True for im, v in t.inflow_dict.items() if (df['Recon3D'] == im[:-3]).sum() == 0]  # get True if not all IN metabolites can be converted to iHuman ID
        OUTB = [True for im, v in t.outflow_dict.items() if (df['Recon3D'] == im[:-3]).sum() == 0] # get True if not all OUT metabolites can be converted to iHuman ID
        if len(INB) != 0 or len(OUTB) != 0:
            t2Remove.append(t.name)
    tasks = [t for t in tasks if t.name not in t2Remove] # remove tasks
    # replace ids of metabolites from recon3d to human1:
    Mt2Remove = list()
    tasksInD = dict()
    tasksOutD = dict()
    for t in tasks: # all tasks to test have reaction_dict empty: [t.name for t in tasks if len(t.reaction_dict) != 0] is []
        inDct = {''.join([im[:-2], 's]']) if (im[-2] == 'x' or  im[-2] == 'e') else im: v for im, v in t.inflow_dict.items()} # replace x and e compartments by s (both correspond to extracellular comp in human1)
        outDct = {''.join([im[:-2], 's]']) if (im[-2] == 'x' or im[-2] == 'e') else im: v for im, v in t.outflow_dict.items()} # replace x and e compartments by s (both correspond to extracellular comp in human1)
        inD = dict()
        outD = dict()
        for im, v in inDct.items(): # for inflow metabolites
            if len(df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]) == 0: # if combination of metabolite id + compartment do not exist in human1
                Mt2Remove.append(t.name) # append task, to latter on remove it
            else: # else, replace recon3d metabolite id by id in human1
                inD[df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]['Human1'].iloc[0]] = v
        for im, v in outDct.items(): # for outflow metabolites
            if len(df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]) == 0: # if combination of metabolite id + compartment do not exist in human1
                Mt2Remove.append(t.name) # append task, to latter on remove it
            else: # else, replace recon3d metabolite id by id in human1
                outD[df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]['Human1'].iloc[0]] = v
        tasksInD[t.name] = inD
        tasksOutD[t.name] = outD
    tasks = [t for t in tasks if t.name not in set(Mt2Remove)] # remove reactions with ombination of metabolite id + compartment that do not exist in human1
    tasksN = copy.copy(tasks)
    for t in tasksN:
        t.inflow_dict = tasksInD[t.name]
        t.outflow_dict = tasksOutD[t.name]
    JSONTaskIO().write_task(finalpath, tasksN)

## preprocess tasks:
    taskModel = model.copy()
    # close opened boundary reactions of the model to evaluate tasks:
    for k in taskModel.boundary:
        k.knock_out()
     = preprocessTasks(TasksJsonPath, taskModel)
    # evaluate tasks on generic model:
    tasksName, batch_res_tasks = EvalAllTasksAtSameTime(taskModel, task_list)
    # failed tasks on universal model will be removed from list of tasks to test in reconstructed model:
    task_list = RemoveFailedTasks(tasksName, batch_res_tasks, task_list, TasksJsonPath)












corr.keys()
dict_keys(['mets', 'metBiGGID', 'metHMDBID', 'metKEGGID', 'metLipidMapsID', 'metMetaNetXID', 'metsNoComp', 'metHMR2ID', 'metHepatoNET1ID', 'metRecon3DID', 'metEHMNID', 'metChEBIID', 'metPubChemID'])


### Pre-process essential tasks and test them on universal model:
task_list = TestTasksUniversalModel(TasksExcelPath, TasksJsonPath, modelMedAdap, BaseDir)
'''
'''




pvals = pd.read_csv(os.path.join(BaseDir, 'results/Pathways/pfbaPval.tab'), sep='\t', index_col=0)
path = os.path.join(BaseDir, 'results/Pathways/pfbaPvalHeatMap.png')
data=pvals
csize=0.6
rsize=0.4
cclust=True
rclust=True
nc = data.shape[1]
nr = data.shape[0]
fcl = 12 / nc
sns.set(font_scale=fcl)
res = sns.clustermap(data, col_cluster=cclust, row_cluster=rclust, cmap='coolwarm', linecolor='black', xticklabels=True, yticklabels=True, linewidths=0.1, figsize=(csize * nc, rsize * nr), dendrogram_ratio=(0.2, 0.05), cbar=True,cbar_kws={"shrink": .82})
res.fig.subplots_adjust(bottom=0.2 * fcl)
res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_ymajorticklabels(), fontsize=200 / nr)
res.savefig(path, dpi=800)
plt.close('all')
'''
