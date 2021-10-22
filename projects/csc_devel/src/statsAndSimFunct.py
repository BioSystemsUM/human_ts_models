### Modules:
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import cobra
import re
from functools import reduce
import pandas as pd
import ast

### Functions:
def subsetCellTypes(mdF, StdD, TissueCode):
    '''
    - given a dataframe with booleans for each model, finds features specific/or common for certain cell types (CSCs and/or CCs).
      a feature is considered to occur in CSCs/CCs of a tissue if >=50% of donnors/cell lines in that tissue contain that feature
    :param mdF: dataframe  with booleans where columns are donors/cell lines of different tissues and rows are tasks/genes/metabolites etc..
    :param StdD: dict with tissue name and corresponding models to use
    :param TissueCode: dict with tissue code vs tissue name
    :return inCSCnotCC: dataframe with dif tissues on columns with 'true' when feature (row) is in CSC but not in CC
    :return inCSCinCC: dataframe with dif tissues on columns with 'true' when feature (row) is in CSC and in CC
    :return inCCnotCSC: dataframe with dif tissues on columns with 'true' when feature (row) is in CC but not CSC
    '''
    CSCspecD = dict()
    commonD = dict()
    CCspecD = dict()
    for std, mds in StdD.items(): # for each study and corresponding models
        ts = TissueCode[std.split('.')[0][1:]]
        CCs = [mdN for mdN in mds if 'CCs' in mdN] # get CCs model names for a study
        CSCs = [mdN.replace('CCs', 'CSCs') + '_CXCR4+MET+CD44+' if mdN.split('.')[0]== 'M9' else  mdN.replace('CCs', 'CSCs') for mdN in CCs]
        CCbol = mdF[CCs].apply(lambda row: True if row.mean() >= 0.5 else False, axis=1) # boolean indicating if feature occurs in at least 50% of donors/cell lines of CCs of that tissue
        CSCbol = mdF[CSCs].apply(lambda row: True if row.mean() >= 0.5 else False, axis=1) # boolean indicating if feature occurs in at least 50% of donors/cell lines of CSCs of that tissue
        # 'True' for CSC and not CC of that tissue:
        CSCspec = set(CSCbol.index[CSCbol]) - set(CCbol.index[CCbol])
        # 'True' for both CSC and CC of that tissue:
        common = set(CSCbol.index[CSCbol]).intersection(set(CCbol.index[CCbol]))
        # 'True' for CC and not CSC of that tissue:
        CCspec = set(CCbol.index[CCbol]) - set(CSCbol.index[CSCbol])
        CSCspecD[ts] = CSCspec
        commonD[ts] = common
        CCspecD[ts] = CCspec
    CSCspecT = reduce(lambda x,y: x.union(y), CSCspecD.values()) # set with all features in at least one donor in a tissue
    commonT = reduce(lambda x, y: x.union(y), commonD.values())  # set with all features in at least one donor in a tissue
    CCspecT = reduce(lambda x, y: x.union(y), CCspecD.values())  # set with all features in at least one donor in a tissue
    inCSCnotCC = pd.DataFrame.from_dict({ts: {el: True if {el}.issubset(g) else False for el in CSCspecT} for ts, g in CSCspecD.items()})
    inCSCinCC = pd.DataFrame.from_dict({ts: {el: True if {el}.issubset(g) else False for el in commonT} for ts, g in commonD.items()})
    inCCnotCSC = pd.DataFrame.from_dict({ts: {el: True if {el}.issubset(g) else False for el in CCspecT} for ts, g in CCspecD.items()})
    return inCSCnotCC, inCSCinCC, inCCnotCSC

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
    dirpath = os.path.join(BaseDir, 'results/Reconst/RecModelBounds')
    model_dicts = dict()
    for file in os.listdir(dirpath):
        df = pd.read_csv(os.path.join(dirpath, file), sep='\t')
        if not df.empty:
            df = df.T
            df.columns = df.iloc[0]
            df = df.drop(df.index[0]).T
            model_dicts.update(df.to_dict())
    mdD = {'_'.join(k.split('_')[1:]):v for k,v in model_dicts.items()} # exclude label of reconstruction algorithm
    modelBds = pd.DataFrame.from_dict(mdD)
    # select just matched donors that had models sucessfully reconstructed and gapfilled:
    modelBds = modelBds[reduce(lambda x,y: x+y, StdD.values())]
    modelBds = modelBds.applymap(lambda x:ast.literal_eval(x))
    # save dataframe:
    StdPath = os.path.join(BaseDir, 'results/Reconst/RecModelBounds.tab')
    modelBds.to_csv(StdPath, sep='\t')
    # when bounds are different from (0,0) replace by 'True' otherwise 'False':
    AllMdDf = modelBds.applymap(lambda x: True if x != (0,0) else False)
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
    # barplot reactions split by tissue and cell type:
    fig, ax = plt.subplots(figsize=(18/2,7/1.5))
    Medbplt = sns.barplot(x='Tissue', y='Number of Reactions', data=JRct, hue='Cell Type')
    handles, labels = Medbplt.get_legend_handles_labels()
    plt.legend(handles, labels, bbox_to_anchor=(1.015, 1), loc=2, borderaxespad=0., fontsize=15/1.5)
    Medbplt.set_xlabel('', fontsize=20/1.5)
    Medbplt.set_ylabel('Number of Reactions', fontsize=20/1.5, labelpad=5)
    Medbplt.set(ylim=(7000, 11500))
    Medbplt.tick_params(labelsize=15/1.5)
    Medbplt.tick_params(axis='x', rotation=90, labelsize=15/1.15)
    fig.tight_layout()
    #plt.show()
    plt.savefig(os.path.join(BaseDir, 'results/Pathways','Rc_barplotissue.png'))
    plt.savefig(os.path.join(BaseDir, 'results/Pathways', 'Rc_barplotissue.svg'), dpi=300)
    plt.close(fig)
    # boxplot:
    fig, ax = plt.subplots(figsize=(12/2,7/2))
    Medbxplt = sns.boxplot(x='Property', y='Number', data=MedDf, hue='Cell Type')
    handles, labels = Medbxplt.get_legend_handles_labels()
    plt.legend(handles, labels, bbox_to_anchor=(1.015, 1), loc=2, borderaxespad=0., fontsize=15/1.5)
    Medbxplt.set_xlabel('', fontsize=5/1.5, labelpad=0)
    Medbxplt.set_ylabel('Number', fontsize=18/1.5, labelpad=2.5)
    Medbxplt.tick_params(axis='x', labelsize=17/1.5)
    Medbplt.tick_params(axis='y', labelsize=12/1.5)
    fig.tight_layout()
    #plt.show()
    plt.savefig(os.path.join(BaseDir, 'results/Pathways','RcGeneMetb_boxplot.png'))
    plt.savefig(os.path.join(BaseDir, 'results/Pathways', 'RcGeneMetb_boxplot.svg'), dpi=300)
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

def medianMdSameStd(MdTotalDf, StdD, TissueCode, axis, datype, calc = 'median'):
    '''
    - given a series/dataframe with values for each model, obtain dict with with each tissue and cell type as key and corresponding median/mean/std value number
    :param MdTotalDf: series/dataframe with values for different models
    :param StdD: dict with study number and corresponding models
    :param TissueCode: dict with tissue id tissue name correspondence
    :param axis: axis across which aplly medium. axis 0 is across rows and axis 1 across columns. for series use axis 0
    :param datype: string indicating MdTotalDf data type. either 'series' or 'dataframe'
    :param calc: calculation to apply: can be 'median' (default), 'mean' or std
    '''
    AvgD = dict()
    for std, mds in StdD.items(): # for each study and corresponding models
        CCsmp = [mdN for mdN in mds if 'CCs' in mdN] # get CCs model names for that study
        CSCsmp = [mdN for mdN in mds if 'CSCs' in mdN] # get CSCs model names for that study
        ts = TissueCode[std.split('.')[0][1:]]
        if calc == 'median':
            if datype == 'series':
                AvgD['.'.join([ts, 'CCs'])] = [MdTotalDf[CCsmp].median(axis=axis)]
                AvgD['.'.join([ts, 'CSCs'])] = [MdTotalDf[CSCsmp].median(axis=axis)]
            elif datype == 'dataframe':
                AvgD['.'.join([ts, 'CCs'])] = MdTotalDf[CCsmp].median(axis=axis)
                AvgD['.'.join([ts, 'CSCs'])] = MdTotalDf[CSCsmp].median(axis=axis)
        if calc == 'mean':
            if datype == 'series':
                AvgD['.'.join([ts, 'CCs'])] = [MdTotalDf[CCsmp].mean(axis=axis)]
                AvgD['.'.join([ts, 'CSCs'])] = [MdTotalDf[CSCsmp].mean(axis=axis)]
            elif datype == 'dataframe':
                AvgD['.'.join([ts, 'CCs'])] = MdTotalDf[CCsmp].mean(axis=axis)
                AvgD['.'.join([ts, 'CSCs'])] = MdTotalDf[CSCsmp].mean(axis=axis)
        if calc == 'std':
            if datype == 'series':
                AvgD['.'.join([ts, 'CCs'])] = [MdTotalDf[CCsmp].std(axis=axis)]
                AvgD['.'.join([ts, 'CSCs'])] = [MdTotalDf[CSCsmp].std(axis=axis)]
            elif datype == 'dataframe':
                AvgD['.'.join([ts, 'CCs'])] = MdTotalDf[CCsmp].std(axis=axis)
                AvgD['.'.join([ts, 'CSCs'])] = MdTotalDf[CSCsmp].std(axis=axis)
    return AvgD

def applyclustermap(data, filePath, csize, rsize):
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
    fcl = 12 / nc
    sns.set(font_scale=fcl)
    colrs = sns.color_palette("dark:white_r", as_cmap=True)
    res = sns.clustermap(data, cmap=colrs, linecolor='black', xticklabels=True, yticklabels=True,
                            linewidths=0.1, figsize=(csize * nc, rsize * nr), dendrogram_ratio=(0.2, 0.03), cbar_pos = (0.9, 0.92, 0.02, 0.05), rasterized=True)
    res.ax_cbar.tick_params(labelsize=16, pad= 10)
    xlabel = res.ax_heatmap.get_xmajorticklabels()
    res.ax_heatmap.set_xticklabels(xlabel, fontsize=15)
    #plt.show()
    res.savefig(filePath, dpi=300)
    plt.close('all')

def MdCmpHeatMap(model, modelPathMat, BaseDir, StdD, TissueCode, AllMdDf):
    '''
    - create heatmap with tissue vs subsystem with values:
      * 1, when on the top 10% subsytems of list of CSC-CC difference
      * 0, otherwise
    - count number of tissues that each subsystem is on top of list have and save in dataframe
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
    path = os.path.join(BaseDir, 'support/Pathways')
    if not os.path.exists(path):
        os.makedirs(path)
    RcSbsF.to_csv(os.path.join(BaseDir, 'support/Pathways/subsytemInfo.tab'), sep='\t') # save subsystem info in a dataframe
    # percentage of reactions in each subsytem that are active in each model:
    AllMdDfSbs = mdDf.groupby('Subsystem').sum()
    NbRcInSbs = mdDf.pivot_table(index=['Subsystem'], aggfunc='size') # number of reactions in each subsytem in universal model
    PercSbsRct = AllMdDfSbs.apply(lambda x: x/NbRcInSbs*100)
    # get median percentage of each subsystem with active reactions, per tissue and cell type (CCs or CSCs):
    MedD = medianMdSameStd(MdTotalDf=PercSbsRct, StdD=StdD, TissueCode=TissueCode, axis=1, datype='dataframe')
    MedDf = pd.DataFrame.from_dict(MedD)
    ## get difference beween CSCs and CCs in each tissue:
    MedDdiff = {ts: MedDf[ts + '.CSCs'] -  MedDf[ts + '.CCs'] for ts in TissueCode.values()}
    MedDfDiff = pd.DataFrame.from_dict(MedDdiff)
    # exclude subsystems with low number of reactions (3 or less):
    CSC_CC_df = MedDfDiff[NbRcInSbs >= 4]
    lnbRc = pd.DataFrame(NbRcInSbs[NbRcInSbs < 4])
    lnbRc.columns = ['Number of Reactions']
    lnbRc.index.name = 'Subsystem with 3 or less reactions'
    lnbRc.to_csv(os.path.join(BaseDir, 'support/Pathways', 'SubsLessThan4Reac'), sep='\t', header=True) # save info on excluded subsytems with low number of reactions (less than 4)
    # exclude 'exchange/demand reactions':
    # (number of these reactions is dif. depending on .sbml and .mat model, besides in this analysis they do not change between reconstructed models)
    CSC_CC_df = CSC_CC_df.drop(index='Exchange/demand_reactions')
    # exclude subsystems where no tissue had difference between CSC and CCs:
    CSC_CC_df = CSC_CC_df.loc[(CSC_CC_df != 0).any(axis=1)]
    # exclude other weird subsytems:
    CSC_CC_df = CSC_CC_df.drop(index='Isolated')
    CSC_CC_df = CSC_CC_df.drop(index='Miscellaneous')
    # get top pathways in each study and then merge info:
    perc = int(CSC_CC_df.shape[0] * 0.10) # number of subsytems in top 10% - top 12
    tCSCL = list()
    for ts in CSC_CC_df: # for each tissue get order subsytems by those with more existent reactions in CSCs to those with more in CCs
        df = CSC_CC_df[ts].sort_values(ascending=False)
        tp = df.unique()[:perc] # subsytems in top 12 value in CSCs are not necessarily 12 subsytems (as some may have the same value)
        tp = tp[tp > 0.5] # guarantee 0 is not among top 12 values in CSCs (as 0 means no difference..)
        topCSC = pd.DataFrame(df[df.isin(tp)])
        tCSCL.append(topCSC)
    tCSCDf = reduce(lambda x,y: x.join(y, how='outer'), tCSCL)
    # get number of tissues where a subsystem is on top of ordered list (on top of those with more existing reactions in CSCs than CCs):
    tCSCDfCp = tCSCDf.copy()
    tCSCDf[tCSCDf.notna()] = 1 # when subsytem has a value in top CSC-CC, counts as 1, otherwise (if nan) counts as 0
    tCSCDf = tCSCDf.fillna(0)
    moreCSCDf = pd.DataFrame(tCSCDf.sum(axis=1).sort_values(ascending=False))
    moreCSCDf.columns = ['Number of Tissues']
    moreCSCDf.to_csv(os.path.join(BaseDir, 'results/Pathways/MoreCSCRank.tab'), sep='\t', header = True)
    # heatmap with existence of subsytems on top of ordered CSC-CC difference list:
    tCSCDfCp[tCSCDfCp.notna()] = 1 # using the copied dataframes above, if subs. is "on top of CSC-CC list" counts as 1
    new = tCSCDfCp.copy()
    new = new.fillna(0) # all NAs (neither on top nor bottom of CSC-CC list) are replaced by 0
    path = os.path.join(BaseDir, 'results/Pathways', 'Top_CSCvsCCDiff.png')
    new.index = map(lambda x: x.replace('_', ' '), new.index)
    applyclustermap(data=new, filePath=path, csize=1.1, rsize=0.2)
    path = os.path.join(BaseDir, 'results/Pathways', 'Top_CSCvsCCDiff.svg')
    applyclustermap(data=new, filePath=path, csize=1.1, rsize=0.2)

def pFBA(BaseDir, modelMedAdap):
    '''
    - do pFBA in reconstructed models and saves info on dataframe
    :param BaseDir: basic directory
    :param modelMedAdap: generic metabolic modle adapted for medium composition
    '''
    path = os.path.join(BaseDir, 'results/Reconst/RecModelBounds.tab')
    modelBds = pd.read_csv(path, sep='\t', index_col=[0])
    modelBds = modelBds.applymap(lambda x: ast.literal_eval(x))
    # verify models grow and are feasible:
    for mdN in modelBds:
        mdBds = modelBds[mdN]
        with modelMedAdap as gmd:
            for rcid, bd in zip(mdBds.index, mdBds):
                if rcid in [r.id for r in gmd.reactions]:
                    gmd.reactions.get_by_id(rcid).bounds = bd
            opt = gmd.optimize()
            if opt.status != 'optimal' or opt.objective_value < 1E-9:
                print(mdN + ' infeasible or does not grow')
    # do pFBA:
    PfbaFluxD = dict()
    for mdN in modelBds:
        print(mdN)
        mdBds = modelBds[mdN]
        with modelMedAdap as gmd:
            for rcid, bd in zip(mdBds.index, mdBds):
                if rcid in [r.id for r in gmd.reactions]:
                    gmd.reactions.get_by_id(rcid).bounds = bd
            pfbasol = cobra.flux_analysis.pfba(gmd)
            pfbaflx = pfbasol.fluxes
            PfbaFluxD[mdN] = pfbaflx  # save model pFBA fluxes to dictionary with model name as key
    df = pd.DataFrame.from_dict(PfbaFluxD)
    pathD = os.path.join(BaseDir, 'results/Simulations')
    if not os.path.exists(pathD):
        os.makedirs(pathD)
    df.to_csv(os.path.join(pathD, 'pFBA.tab'), sep ='\t')

def prepPfba2PvalCalc (BaseDir, StdD, TissueCode):
    '''
    - apply module to pFBA fluxes, get median for each tissue and cell type and get difference between CSCs and CCs in each tisue; save results
    :param BaseDir: basic directory
    :param StdD: dictionary with study names as keys and corresponding model names as values
    :param TissueCode: dictionary with tissues ids and corresponding tissues' names
    '''
    # get pfba flux results:
    pfba_flux = pd.read_csv(os.path.join(BaseDir, 'results/Simulations/pFBA.tab'), sep='\t', index_col=[0])
    # apply absolute value to each flux:
    AbsPfba = pfba_flux.applymap(abs)
    # cal. median for each tissue and cell type:
    tsAbsPfba = medianMdSameStd(AbsPfba, StdD, TissueCode, axis=1, datype='dataframe')
    tspfba = pd.DataFrame.from_dict(tsAbsPfba)
    # get difference between CSCs and CCs in each tissue:
    pfbadifD = {ts: tspfba[ts + '.CSCs'] -  tspfba[ts + '.CCs'] for ts in TissueCode.values()}
    pfbadif = pd.DataFrame.from_dict(pfbadifD)
    pathD = os.path.join(BaseDir, 'support/Simulations')
    if not os.path.exists(pathD):
        os.makedirs(pathD)
    pfbadif.to_csv(os.path.join(pathD, 'pfbaAbsValDiffTss.tab'), sep='\t')

def heatMpPvpFBA(BaseDir):
    '''
    plots heatmap where:
     * positive values - subsystems enriched in top of list with difference between CSCs and CCs abs(pFBA fluxes)
                         i.e enriched for CSCs but not CCs
                         value is -log(p-value) (the higher the value, the more significant)
     * negative values - subsytems enriched in bottom of list with difference between CSCs and CCs abs(pFBA fluxes)
                         i.e. enriched for CCs but not CSCs
                         value is log(p-value) (the lower the value, the more significant)
     The above is for subsytems with significant p-value (p-value < 0.05)
     All remaining subsystems have value 0 in the plot.
    :param BaseDir: basic directory
    '''
    # load dataframe with subsystems on rows, tissue on columns where values are -/+log(p-value):
    pvalsSig = pd.read_csv(os.path.join(BaseDir, 'support/Simulations/pfbaSubsDiffPvalSign.tab'), sep='\t', index_col=0)
    pvalsSig.rename(columns={'head.and.neck':'head and neck'}, inplace=True)
    pathSig = os.path.join(BaseDir, 'results/Pathways/pfbaPvalHeatMapSig.svg') # where to save plot
    pvalsSig.index = map(lambda x: x.replace('_', ' '), pvalsSig.index)
    # exclude 'exchange/demand reactions':
    # (number of these reactions is dif. depending on .sbml and .mat model, besides in this analysis they do not change between reconstructed models)
    if 'Exchange/demand reactions' in pvalsSig.index:
        pvalsSig = pvalsSig.drop(index='Exchange/demand reactions')
    # exclude other weird subsytems:
    if 'Isolated' in pvalsSig.index:
        pvalsSig = pvalsSig.drop(index='Isolated')
    if 'Miscellaneous' in pvalsSig.index:
        pvalsSig = pvalsSig.drop(index='Miscellaneous')
    data=pvalsSig
    csize=0.9
    rsize=0.2
    nc = data.shape[1]
    nr = data.shape[0]
    fcl = 12 / nc
    sns.set(font_scale=fcl)
    colrs = sns.color_palette("dark:white_r", as_cmap=True)
    res = sns.clustermap(data, cmap=colrs, linecolor='black', xticklabels=True, yticklabels=True, linewidths=0.01, figsize=(csize * nc, rsize * nr), dendrogram_ratio=(0.2, 0.05), cbar_pos = (0.77, 0.90, 0.01, 0.05), rasterized = True, cbar_kws={'ticks':[int(data.min().min()), (data.max().max()-0.2)]})
    hm = res.ax_heatmap.get_position()
    col = res.ax_col_dendrogram.get_position()
    res.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*1.1, hm.height])
    res.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*1.1, col.height*0.95])
    res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_ymajorticklabels(), fontsize=230 / nr)
    res.ax_heatmap.set_xticklabels(res.ax_heatmap.get_xmajorticklabels(), fontsize=255*1.2 / nr)
    res.ax_cbar.tick_params(labelsize=8, pad= 3)
    res.ax_cbar.set_yticklabels([str(int(data.min().min())), str(int(data.max().max()))])
    #plt.show()
    res.savefig(pathSig, dpi=300)
    plt.close('all')




