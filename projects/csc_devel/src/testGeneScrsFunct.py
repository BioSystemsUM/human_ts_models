### Modules:
import os
import pandas as pd
import seaborn as sns
import numpy as np
import pickle
import itertools
import cobra
from functools import reduce
import matplotlib.pyplot as plt
from scipy.spatial import distance
from troppo.omics.core import OmicsContainer
from cobamp.utilities.parallel import batch_run
from sklearn.preprocessing import MinMaxScaler
from troppo.methods_wrappers import ReconstructionWrapper
from multiprocessing import cpu_count
from troppo.omics.core import IdentifierMapping, TypedOmicsMeasurementSet
from itertools import product, chain
from multiprocessing.dummy import Pool
from functools import partial
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from faker import Factory

### Functions:

def scaleGenExp(modelProcPath, modelMedAdapPath, BaseDir):
    '''
    - applies min-max normalization to individual studies' gene expression data (to all genes/ just genes in generic metabolic model)
      saves to .csv files and creates clustermaps
    - returns generic model and model adapted for medium composition
    :param modelProcPath: param of preprocessiHuman function
    :param modelMedAdapPath: param of preprocessiHuman function
    :param BaseDir: base directory
    :return model: generic metabolic model
    :return modelMedAdap: metabolic model adapted for medium composition
    '''
    # get list with all genes in medium adapted metabolic model, human1:
    model = cobra.io.read_sbml_model(modelProcPath)
    modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)
    geneLst = list({g.id for g in modelMedAdap.genes})
    # for each study pre-process gene expression data:
    Datapath = os.path.join(BaseDir, 'data/studies')
    Dirs = [files for files in os.listdir(Datapath) if os.path.isdir(os.path.join(Datapath, files))]
    #sys.setrecursionlimit(1000000) # needed to plot clusterplot with high number of rows
    for StudyNumber in Dirs:
        print(StudyNumber)
        StudyDir = os.path.join(BaseDir,'data/studies/' + str(StudyNumber)) # directory with data of study
        ResStudyDir = os.path.join(BaseDir,'results/studies/' + str(StudyNumber)) # directory where to save results of study
        SupStudyDir = os.path.join(BaseDir,'support/studies/' + str(StudyNumber)) # directory where to save intermediary/supporting files
        if not os.path.exists(ResStudyDir):
            os.makedirs(ResStudyDir)
        if not os.path.exists(SupStudyDir):
            os.makedirs(SupStudyDir)
        GeneExpressionPath = os.path.join(StudyDir, 'NormData', 'NormGeneExp_Groups.tab')
        studyExpVal = pd.read_csv(GeneExpressionPath, sep='\t')
        studyExpVal = studyExpVal.iloc[:,[0] + list(range(3, studyExpVal.shape[1]))]
        names = studyExpVal.columns.tolist()
        samples = [StudyNumber + '.' +  x for x in names[1:]]
        samples.insert(0,names[0])
        studyExpVal.columns = samples
        studyExpVal.set_index('Ensembl', inplace=True)
        studyExpVal = studyExpVal.dropna() # remove rows with NaN if they exist
        cLstudyExpVal = studyExpVal.copy()
        # drop genes with 0s for all samples:
        excBl = cLstudyExpVal.sum(axis=1) == 0
        exc = cLstudyExpVal[excBl].index
        cLstudyExpVal.drop(exc, axis=0, inplace=True)
        # select genes in metabolic model:
        sdExpMetb = cLstudyExpVal[cLstudyExpVal.index.map(lambda x: x in geneLst)]
        # apply min-max normalization just to genes in metabolic model:
        TExpMetb = sdExpMetb.T
        scaler = MinMaxScaler()
        scaledMtb = scaler.fit_transform(TExpMetb)
        fnTMtb = pd.DataFrame(scaledMtb, columns=TExpMetb.columns, index=TExpMetb.index)
        fnMtb = fnTMtb.T
        fnMtb.index.name = ' '
        # clustermap of genes of human model with min-max normalization:
        ncMtb = fnMtb.shape[1]
        nrMtb = fnMtb.shape[0]
        fclMtb = 15/ncMtb/2.3
        sns.set(font_scale=fclMtb)
        resMtb = sns.clustermap(fnMtb,  cmap='coolwarm', linecolor='black', xticklabels=True, yticklabels=False, linewidths=0.00001, figsize=(0.77/4*ncMtb, 0.0065/4*nrMtb), dendrogram_ratio=(0.2,0.05),  cbar_pos = (0.91 - 0.05, 0.90, 0.02, 0.05))
        hm = resMtb.ax_heatmap.get_position()
        resMtb.ax_heatmap.set_position([hm.x0, hm.y0 + 0.05, hm.width*0.90, hm.height*0.95])
        col = resMtb.ax_col_dendrogram.get_position()
        resMtb.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.90, col.height*0.95])
        row = resMtb.ax_row_dendrogram.get_position()
        resMtb.ax_row_dendrogram.set_position([row.x0, row.y0 + 0.05, row.width, row.height*0.95])
        resMtb.ax_cbar.tick_params(labelsize=14/2.3, pad= 5/3)
        #plt.show()
        resMtb.savefig(os.path.join(ResStudyDir, 'explCluster_MetbGenes.png'), dpi=300)
        resMtb.savefig(os.path.join(ResStudyDir, 'explCluster_MetbGenes.svg'), dpi=300)
        '''
        # without min-max norm:
        sdExpMetb.index.name = ' '
        ncMtb = sdExpMetb.shape[1]
        nrMtb = sdExpMetb.shape[0]
        fclMtb = 15/ncMtb/2.3
        sns.set(font_scale=fclMtb)
        resMtb = sns.clustermap(sdExpMetb,  cmap='coolwarm', linecolor='black', xticklabels=True, yticklabels=False, linewidths=0.00001, figsize=(0.77/4*ncMtb, 0.0065/4*nrMtb), dendrogram_ratio=(0.2,0.05),  cbar_pos = (0.84-0.06, 0.90, 0.02, 0.05), cbar_kws={'ticks':[sdExpMetb.min().min(), np.round(sdExpMetb.mean().mean(), 2), int(sdExpMetb.max().max())]})
        hm = resMtb.ax_heatmap.get_position()
        resMtb.ax_heatmap.set_position([hm.x0 - 0.05, hm.y0 + 0.05, hm.width*0.90, hm.height*0.95])
        col = resMtb.ax_col_dendrogram.get_position()
        resMtb.ax_col_dendrogram.set_position([col.x0 - 0.05, col.y0, col.width*0.90, col.height*0.95])
        row = resMtb.ax_row_dendrogram.get_position()
        resMtb.ax_row_dendrogram.set_position([row.x0, row.y0 + 0.05, row.width*0.8, row.height*0.95])
        resMtb.ax_cbar.tick_params(labelsize=14/2.3, pad= 5/3)
        resMtb.ax_cbar.set_yticklabels([str(sdExpMetb.min().min()), str(''), str(int(sdExpMetb.max().max()))])
        #plt.show()
        resMtb.savefig(os.path.join(ResStudyDir, 'explCluster_MetbGenes_noNorm.png'), dpi=300)
        resMtb.savefig(os.path.join(ResStudyDir, 'explCluster_MetbGenes_noNorm.svg'), dpi=300)
        '''
        # save gene expression normalized wih min-max but including 0s, for genes in generic metabolic model:
        ExpMetbZeros = studyExpVal[studyExpVal.index.map(lambda x: x in geneLst)]  # select genes that are in generic metabolic model
        TExpMetbZeros = ExpMetbZeros.T
        scaledMtbZeros = scaler.fit_transform(TExpMetbZeros)
        fnTMtbZeros = pd.DataFrame(scaledMtbZeros, columns=TExpMetbZeros.columns, index=TExpMetbZeros.index)
        fnTMtbZeros.to_csv(os.path.join(SupStudyDir, 'MinMaxNormExp_MetbGenes.csv'))
        # save gene expression normalized wih min-max but including 0s, for all genes (even those not in generic metabolic model):
        TExpAllZeros = studyExpVal.T
        scaledAllZeros = scaler.fit_transform(TExpAllZeros)
        fnTAllZeros = pd.DataFrame(scaledAllZeros, columns=TExpAllZeros.columns, index=TExpAllZeros.index)
        fnTAllZeros.to_csv(os.path.join(SupStudyDir, 'MinMaxNormExp_AllGenes.csv'))
    return model, modelMedAdap

def TestScores(expC, sample, GenExpName):
    '''
    - determines gene scores for different expression thresholds
    :param expC: gene expression values (pre-normalized or not) of one study. columns are genes and rows are samples
    :param sample: name of sample for which gene scores are being calculated
    :param GenExpName: string with name of file with expression data, from which we retrieve info on whether all genes of just metabolic genes are analysed
    :return: list of dictionaries where each dict corresponds to a threshold strategy (global, local1, local1B, local2, local2B)
    '''
    IdGenesType = GenExpName.split('_')[1].split('.')[0] # string indicates if using all genes or just metabolic genes to determine thresholds
    qvalues = [0.1, 0.25, 0.5, 0.75, 0.9]
    quantiles = expC.quantile(qvalues)  # get quantile values for each gene across all samples. rows are quantiles, columns are genes
    ### global strategy - one global threshold:
    global_thresholds = quantiles.T.apply(lambda x: x.mean())  # each global threshold is the mean of a type of quantile values across all genes
    sample_series = expC.loc[sample, :]  # subset expression values representing a sample
    global_dicts = dict()
    for i, g in enumerate(global_thresholds):
        maxExp = (expC/g).apply(np.log).max().max()
        globAct = (sample_series / g).apply(np.log).clip(lower=-maxExp, upper=maxExp)
        global_dicts[(IdGenesType, 'global', i, None, None,)] = globAct.to_dict()
        # divided expression values of sample by threshold and logaritmized.
        # 'clip()' avoids -Inf/Inf values. -maxExp is used to avoid very low negative values, which would otherwise penalize too much low expressed genes when doing 'tinit'
        # negative values: represent genes with scores bellow global threshold, and positive values: genes which scores are above threshold
        # global_dicts is dict where keys are id of a global threshold (a number) and values are dict.
        # lower level dict has 'geneSymbol (ENSGENE)' as key and as value a gene score
    ### local T1 strategy - one global threshold and one local threshold:
    param_combinations = list(product(*[range(len(qvalues))] * 2)) # [range(len(quantiles))]*2 -> it is [range(0, 5), range(0, 5)]
        # list of tupples where 1st element represents id of global and the 2nd element the id of local threshold
        # each number represents a quantile: 0 -> quantile 0.1; 1 -> quantile 0.25; 2 -> quantile 0.5; 3 -> quantile 0.75; 4 -> quantile 0.9
    ## giving equal weight to "off" from local and global thresholds:
    local1_dicts = {}
    for k, v in param_combinations:
        gt, lt = global_thresholds.iloc[k], quantiles.iloc[v, :]  # selects a combination of actual global and local thresholds
        maxExp = (expC / gt).apply(np.log).max().max()
        activity = (sample_series / gt).apply(np.log) # genes with expression bellow global threshold are deemed 'off'', so get negative scores
        gt_active = activity >= 0 # genes with expression equal or above global threshold are 'maybe on' and score for local threshold is determined
        activity[gt_active] = (sample_series[gt_active] / lt[gt_active]).apply(np.log) # if a 'maybe on' gene gets a positive score for local threshold, then gene is 'on'
        activity = activity.clip(lower=-maxExp, upper=maxExp) # 'clip()' avoids -Inf/Inf values. -maxExp is used to avoid very low negative values, which would otherwise penalize too much low expressed genes when doing 'tinit'
        local1_dicts[(IdGenesType, 'local1', k, None, v)] = activity.to_dict()
        # dict with names of local and threshold combination used as keys and as values a dict.
        # lower level dict has genes as keys and gene scores as values
    ## giving more weight to "off" from global than local thresholds:
    local1B_dicts = {}
    for k, v in param_combinations:
        gt, lt = global_thresholds.iloc[k], quantiles.iloc[v,:]
        maxExp = (expC / gt).apply(np.log).max().max()
        activity = (sample_series / gt).apply(np.log) - 1 # genes with expression bellow global threshold are deemed 'off'', so get scores < -1
        gt_active = activity >= -1  # genes with expression equal or above global threshold are 'maybe on' and score for local threshold is determined
        activity[gt_active] = (sample_series[gt_active] / lt[gt_active]).apply(np.log)  # if a 'maybe on' gene gets a positive score for local threshold, then gene is 'on'
        activity[gt_active] = activity[gt_active].clip(lower=-1, upper=maxExp) # genes that are "off" according with local threshold can't have sores more negative than genes "off" according with global threshold
        activity = activity.clip(lower= -maxExp - 1, upper=maxExp)
        local1B_dicts[(IdGenesType, 'local1B', k, None, v)] = activity.to_dict()
    ### local T2 strategy - two global thresholds and one local threshold:
    global_lt2_params = list(product(range(len(qvalues)), list(zip(*np.where(np.fromfunction(lambda i, j: i < j, [len(qvalues)] * 2))))))
        # list of tuples where 1st element represents a local threshold and 2nd element is a tuple representing the two global thresholds
        # note that for the two global thresholds names we don't consider order - (0,1) or (1,0) is same so we just have (0,1)
        # each number represents a quantile: 0 -> quantile 0.1; 1 -> quantile 0.25; 2 -> quantile 0.5; 3 -> quantile 0.75; 4 -> quantile 0.9
    ## giving equal weight to "on/off" from local and global thresholds:
    local2_dicts = {}
    for v, k in global_lt2_params:
       gtl, gtu, lt = global_thresholds.iloc[k[0]], global_thresholds.iloc[k[1]], quantiles.iloc[v,:]  # selects combination of two global and one local threshold
       maxExp = (expC / gtu).apply(np.log).max().max()
       upp_activity = (sample_series/gtu).apply(np.log) # genes with expression equal or above global upper threshold are "on"
       gtu_inactive = upp_activity < 0 # genes with expr bellow upper global threshold are considered potential 'maybe on'
       low_activity = (sample_series/gtl).apply(np.log) # genes with expression bellow global lower threshold are "off"
       gtl_maybes, gtl_lows = (low_activity >= 0) & gtu_inactive, low_activity < 0 # gtl_maybes are genes 'maybe on' - with expr between lower and upper global thresholds
       upp_activity[gtl_lows] = low_activity[gtl_lows]
       activity_maybe = (sample_series[gtl_maybes]/lt[gtl_maybes]).apply(np.log) # local threshold determines whether genes 'maybe on' are in fact 'on' (positive score) or 'off' (negative score)
       upp_activity[activity_maybe.index] = activity_maybe
       upp_activity = upp_activity.clip(lower=-maxExp, upper=maxExp) # to avoid -inf/Inf and at same time avoid very low negative values, which would otherwise penalize too much low expressed genes when doing 'tinit'
       local2_dicts[(IdGenesType,'local2',k[0],k[1],v)] = upp_activity.to_dict()
        # dict with names of local and global thresholds combination used as keys and as values a dict.
        # lower level dict has genes as keys and gene scores as values
    ## giving more weight to "on/off" from global than local thresholds:
    local2B_dicts = {}
    for v, k in global_lt2_params:
       gtl, gtu, lt = global_thresholds.iloc[k[0]], global_thresholds.iloc[k[1]], quantiles.iloc[v,:]  # selects combination of two global and one local threshold
       maxExp = (expC / gtu).apply(np.log).max().max()
       upp_activity = (sample_series/gtu).apply(np.log) + 1 # genes with expression equal or above global upper threshold (>1) are "on"
       gtu_inactive = upp_activity < 1 # genes with expr bellow upper global threshold are considered potential 'maybe on'
       low_activity = (sample_series/gtl).apply(np.log) - 1 # genes with expression bellow global lower threshold (<-1) are "off"
       gtl_maybes, gtl_lows = (low_activity >= -1) & gtu_inactive, low_activity < -1 # gtl_maybes are genes 'maybe on' - with expr between lower and upper global thresholds
       upp_activity[gtl_lows] = low_activity[gtl_lows]
       activity_maybe = (sample_series[gtl_maybes]/lt[gtl_maybes]).apply(np.log) # local threshold determines whether genes 'maybe on' are in fact 'on' (positive score) or 'off' (negative score)
       activity_maybe = activity_maybe.clip(lower=-1, upper=1) # genes 'on' by local threshold get scores between [0, 1], genes 'off' by local threshold get scores between [0, -1], genes 'on' by upper global threshold get scores > 1, genes 'off' by lower global threshold get scores < -1
       upp_activity[activity_maybe.index] = activity_maybe
       upp_activity = upp_activity.clip(lower=-maxExp - 1, upper=maxExp + 1) # to avoid -inf/Inf and at same time avoid very low negative values, which would otherwise penalize too much low expressed genes when doing 'tinit'
       local2B_dicts[(IdGenesType,'local2B',k[0],k[1],v)] = upp_activity.to_dict()
    thrStrLst = [global_dicts, local1_dicts, local1B_dicts, local2_dicts, local2B_dicts]
    return thrStrLst

def CombineThreshAND_OR_rules(exp_data):
    '''
    - Create dictionary with all combinations of threshold strategies + and/or rules
    :param exp_data: dataframe with gene scores. columns are genes and rows are thresholds (multiindex)
    :return: dict where key is id of combination of threshold strategy with and/or functions and value is a tuple. tupple has as first value another tuple ('and' function, 'or' function) and 2nd value is a dict where key is gene id and value is gene score
    '''
    # Create an identifier mapping object - only used if gene ids need to be converted to model gene id (ensembl)
    mapping = IdentifierMapping('human_transcriptomics', pd.read_csv('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt', index_col=0, sep = '\t'))
    # Create an omics measurement set object with the gene scores dataframe:
    omics_mset = TypedOmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values, mapping)
    data_dicts = {'_'.join(map(str,k)):v for k,v in omics_mset.data.T.to_dict().items()}
        # dict where key is a gene threshold strategy and values are other dict with: ensembl gene id vs gene score
    # Define a dict of functions to apply AND/OR rules when generating reaction scores:
    #funcs = {'minsum':((lambda x: min([max(y, 0) for y in x])), sum), 'minmax':((lambda x: min([max(y, 0) for y in x])), max)}
    funcs = {'minsum':(min, sum),'minmax':(min,max)}
      # dict where keys are 'minsum' or 'minmax' and values are corresponding functions.
      # for eg. 'minsum' is a tuple (function of minimum, function of sum)
    # Get dictionary with all combinations of threshold strategies + and/or rules:
    runs = dict(chain(*[[((k,n),(aof, v)) for k,v in data_dicts.items()] for n, aof in funcs.items()]))
        # dict where key is id of combination of threshold strategy with and/or functions and value is a tuple
        # tupple has as first value another tuple ('and' function, 'or' function) and 2nd value is a dict where key is gene id and value is gene score
    return runs

def saveResult(res, path):
    '''
    - saves results as pickle object
    '''
    file = open(path, 'wb')
    pickle.dump(res, file)
    file.close()

def ReacScoresFunc(score_tuple, params):
    '''
    - funct used in getRCScores function
    :param score_tuple: (('and' function, 'or' function), dict with gene ids as keys and gene score as values)
    :param params:
    :return:
    '''
    aofx, data_dict = score_tuple
        # aofx is tupple with and/or rule functions
        # data_dict is dict with gene id vs gene score of a combination to test
    # create omics container:
    oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x', data=data_dict, nomenclature='custom')
    rw = params['rw']  # load parameters
    mr = rw.model_reader
    reacScores = oc_sample.get_integrated_data_map(model_reader=mr, and_func=aofx[0], or_func=aofx[1])
    return reacScores

def getRCScores(runs, StdReacScores, sample, rw, NTHREADS):
    '''
    - determines reaction scores for a combination of thresholds + and/or rules
    :param runs: dict where key is id of combination of threshold strategy with and/or functions and value is a tuple.
                 tupple has as first value another tuple ('and' function, 'or' function) and 2nd value is a dict where key is gene id and value is gene score
    :param StdReacScores: dict that is going to be updated with gene scores (for all thresholds) of one sample
    :param sample: sample name
    :param rw: reconstruction wrapper of generic model
    :param NTHREADS: number of threads available
    :return: updated dict with sample ids and inner dict has threshold+and/or rule as keys and values are inner dicts with reaction ids as keys and reaction scores as values
    '''
    labs, iters = zip(*runs.items())
    # labs is labels of combination to test
    # iters is tuple of tupples.
    # lower tupple is ((function 'and', function 'or'), dict with gene id vs gene score)
    # get reaction scores:
    output = batch_run(ReacScoresFunc, iters, {'rw': rw}, threads=min(len(runs), NTHREADS)) # get a list where each element corresponds to a threshold+and/or rule combination
    Lstout = list()
    for out in output: # out corresponds to one threshold and/or rule combination
        scD = dict()
        for rc, sc in out.get_scores().items(): # for each reaction id and reaction score
            scD[rc] = sc
        Lstout.append(scD)
    # create dictionary with threshold combinations dictionaries,
    # inner dictionaries have reation ids as keys and reaction scores as values:
    batch_Rsc_res = dict(zip(labs, Lstout))
    StdReacScores[sample] = batch_Rsc_res
    return StdReacScores

def GetThresAndRcScores(BaseDir, model, GenExpAllName, GenExpMetbName):
    Datapath = os.path.join(BaseDir, 'support/studies')
    Dirs = [files for files in os.listdir(Datapath) if os.path.isdir(os.path.join(Datapath, files))]
    rw = ReconstructionWrapper(model, ttg_ratio=9999)
    NTHREADS = int(cpu_count() - (cpu_count()/4))
    StdReacScores_EachStudy = dict()
    for StudyNumber in Dirs:
        StudyDir = os.path.join(BaseDir,'support/studies/' + str(StudyNumber)) # directory with data of study
        GeneExpressionPathAll = os.path.join(StudyDir, GenExpAllName)
        GeneExpressionPathMetb = os.path.join(StudyDir, GenExpMetbName)
        studyExpValAll = pd.read_csv(GeneExpressionPathAll, index_col=0)
        studyExpValMetb = pd.read_csv(GeneExpressionPathMetb, index_col=0)
        samples = studyExpValAll.index.tolist()
        for sample in samples:
            thrStrLstAll = TestScores(studyExpValAll, sample, GenExpAllName) # determine gene scores for different expression thresholds  - for all genes
            thrStrLstMetb = TestScores(studyExpValMetb, sample, GenExpMetbName)  # determine gene scores for different expression thresholds  - for metabolic genes only (those in generic metabolic model)
            thrStrLst = thrStrLstAll + thrStrLstMetb # join list of dicts with gene scores for all genes and for metabolic genes only
            expData = pd.concat([pd.DataFrame.from_dict(d).T for d in thrStrLst], sort=True) # dataframe where each row is a threshold (eg:Allgenes global 0 NaN NaN, local2 3 4.0 0.0) and each column is a gene
            runs_EachStudy = CombineThreshAND_OR_rules(expData)# dict where key is id of combination of threshold strategy with and/or functions and value is a tuple. tupple has as first value another tuple ('and' function, 'or' function) and 2nd value is a dict where key is gene id and value is gene score
            StdReacScores_EachStudy = getRCScores(runs_EachStudy, StdReacScores_EachStudy, sample, rw, NTHREADS) # dict with sample ids and inner dict has threshold+and/or rule as keys and values are inner dicts with reaction ids as keys and reaction scores as values
    if not os.path.exists(os.path.join(BaseDir, 'support/TestGeneScores')):
        os.makedirs(os.path.join(BaseDir, 'support/TestGeneScores'))
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'StdReacScores.pkl')
    saveResult(res=StdReacScores_EachStudy, path=path)

def GetDictOfReacScoresDF(BaseDir):
    '''
    - save dictionary with combination of thresholds + and/or rules as key and dataframe as value,
      dataframe has react ids as rows and samples as columns, values are reaction scores
    - saves also same output but with median values replacing NAs
    :param BaseDir: base directory
    :return: above mentioned dict
    '''
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'StdReacScores.pkl')
    infile = open(path,'rb')
    StdReacScores_EachStudy = pickle.load(infile)
    infile.close()
    finalRcD = dict()
    finalRcDMed = dict()
    exampleSmp = list(StdReacScores_EachStudy.keys())[0]
    thrs = StdReacScores_EachStudy[exampleSmp].keys()
    for thr in thrs:
        odfLst = list()
        for sample,v in StdReacScores_EachStudy.items():
            ndf = pd.DataFrame.from_dict(v[thr], orient='index', columns=[sample])
            odfLst.append(ndf)
        odf = reduce(lambda x,y: x.join(y, how='outer'), odfLst)
        # exclude genes with NaNs across all samples - when GPR rule doesn't exist for that reaction:
        odf.dropna(how='all', inplace=True)
        # replace remaining NaN values by median of corresponding gene expression over all remaining samples,
        # otherwise we couldn't determine distance between combinations of samples 2 by 2 (as done below) nor do PCA when there is NaN values:
        odfMed = odf.T.fillna(odf.T.median(axis=0),axis=0).T
        finalRcDMed[thr] = odfMed # save with NAs replaced by median: to use for eucledian dist calculation and PCA
        finalRcD[thr] = odf # save with NAs: to use when actually building models
    Finalpath = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct_medianCorr.pkl')
    saveResult(res=finalRcDMed, path=Finalpath)
    Finalpath = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct.pkl')
    saveResult(res=finalRcD, path=Finalpath)

def FindEucDist(GrpDct, finalRcDMed):
    '''
    - funct to be used inside EucDstWithinSampleGroups and SimInf2ObservedDs functions
    :param GrpDct: dict with group ids as keys and corresponding model ids as values
    :param finalRcDMed: dict with thresholds as keys and reaction scores dataframes as values
    '''
    GrpCbDct = {id: [cb for cb in itertools.combinations(g, 2)] for id, g in GrpDct.items()} # 2 by 2 combinations of samples in each group of samples
    ThrDistDic = dict()
    for thr, thrDf in finalRcDMed.items(): # for each threshold and corresponding dataframe (with reaction scores for each sample)
        DistDct = dict()
        for id, g in GrpCbDct.items(): # for each sample group and corresponding 2 by 2 sample combinations
            Lst = list()
            for cb in g:
                Lst.append(distance.euclidean(thrDf[cb[0]], thrDf[cb[1]])) # determine eucledian distance for each 2 by 2 sample combination
            DistDct[id] = pd.Series(Lst).mean() # mean eucledian distance of sample combinations of same group of samples
        ThrDistDic[thr] = DistDct
    ThrDistDf = pd.DataFrame.from_dict(ThrDistDic)
    return ThrDistDf

def EucDstWithinSampleGroups(BaseDir):
    '''
    - create dataframe with average eucledian distance of reaction scores between dif. donors/cell lines of same cell type (CCs/CSCs) and study/tissue
    :param BaseDir: base directory
    :return: saves dataframe with average eucledian distance, where columns are thresholds and rows are groups
    '''
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct_medianCorr.pkl') # dict with thresholds as keys and dataframe with react scores as values, where NAs where replaced by median gene expression
    infile = open(path,'rb')
    finalRcDMed = pickle.load(infile)
    infile.close()
    # create a dataframe with thresholds as columns and groups as rows, where values are average of real euclidean distance between samples of same group:
    GrpDf = pd.read_excel(os.path.join(BaseDir, 'data/expectedSampleGroups.xlsx'), sheet_name='groups_by_study_cell_type')
    GrpDct = {el : list(GrpDf[el].dropna()) for el in GrpDf.columns}
    ThrDistDf = FindEucDist(GrpDct, finalRcDMed)
    ThrDistDf = ThrDistDf.loc[GrpDf.columns,list(finalRcDMed.keys())] # to order rows and columns
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'MeasuredDistanceDf.pkl')
    saveResult(res=ThrDistDf, path=path)
    return GrpDct

def SimInf2ObservedDst(AllGrpSmp, GrpLenDctR, GrpLenDctM, ThrDistDf, finalRcDMed):
    '''
    - to be used by function SimEucDistWithinSampleGroups
    - for one random simulation gets dataframe with 'True' values when average simulated distances are inferior to observed distances of corresponding group
    and gets average of distance between samples of same simulated groups
    :param AllGrpSmp: zip of one list of randomly ordered samples of rnaseq and of one list of randomly ordered samples of microarray
    :param GrpLenDctR: dict with group ids as keys and group length as values - for groups with rnaseq samples
    :param GrpLenDctM: dict with group ids as keys and group length as values - for groups with microarray samples
    :param ThrDistDf: dataframe with thresholds as columns and groups as rows with averaged observed distances (between reaction scores) for each group as values
    :param finalRcDMed: dict with thresholds as keys and dataframes with reaction scores as values
    :return InfDistDf: dataframe with 'True' values when average simulated distances are inferior to observed distances of corresponding group
    :return SimThrDistDf: dataframe with thresholds as columns and groups as rows, where values are average of distance between samples of same simulated group
    '''
    AllGrpSmpR = AllGrpSmp[0]
    AllGrpSmpM = AllGrpSmp[1]
    SimGrpDctR = dict()  # get dict with group number as key and as values we have a list of specific length of random sample names - for groups with rnaseq samples
    i = 0
    for id, v in GrpLenDctR.items():
        f = i + v
        SimGrpDctR[id] = AllGrpSmpR[i:f]
        i = i + v
    SimGrpDctM = dict()  # get dict with group number as key and as values we have a list of specific length of random sample names - for groups with microarray samples
    i = 0
    for id, v in GrpLenDctM.items():
        f = i + v
        SimGrpDctM[id] = AllGrpSmpM[i:f]
        i = i + v
    # create a dataframe with thresholds as columns and groups as rows, where values are average of distance between samples of same simulated group:
    SimThrDistDfM = FindEucDist(SimGrpDctM, finalRcDMed)
    SimThrDistDfR = FindEucDist(SimGrpDctR, finalRcDMed)
    SimThrDistDf = pd.concat([SimThrDistDfM, SimThrDistDfR])
    SimThrDistDf = SimThrDistDf.loc[ThrDistDf.index, ThrDistDf.columns]  # order rows/columns of dataframe of simulated distances according to dataframe of real distances - to latter on stack them with matching positions
    # create dataframe with 'True' values when average simulated distances are inferior to observed distances of corresponding group:
    InfDistDf = SimThrDistDf < ThrDistDf
    return InfDistDf, SimThrDistDf

def SimEucDistWithinSampleGroups(BaseDir, GrpDct, RGrps, MGrps, maxSim = 1000):
    '''
    - save list of dataframes where each dataframe represents a simulation and 'True' values exist when average simulated distances for a group are inferior to observed distances for that group
    - save list of dataframes where each dataframe represents a simulation and values are average of distance between samples of same simulated group
    :param BaseDir: base directory
    :param GrpDct: dict with group id as keys and corresponding donors/cell lines as values
    :param RGrps: list of ids for groups with rnaseq samples
    :param MGrps: list of ids for groups with microarray samples
    :param maxSim: number of simulations
    '''
    pathf1 = os.path.join(BaseDir, 'support/TestGeneScores', 'SimulatedDistanceDfList.pkl')
    pathf2 = os.path.join(BaseDir, 'support/TestGeneScores', 'SimInferior2RealDistanceDfList.pkl')
    if (not os.path.exists(pathf1)) or (not os.path.exists(pathf2)):
        path = os.path.join(BaseDir, 'support/TestGeneScores', 'MeasuredDistanceDf.pkl') # dataframe with averaged observed distances in each group and threshold
        infile = open(path,'rb')
        ThrDistDf = pickle.load(infile)
        infile.close()
        path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct_medianCorr.pkl') # dict with thresholds as keys and dataframe with react scores as values, where NAs where replaced by median gene expression
        infile = open(path,'rb')
        finalRcDMed = pickle.load(infile)
        infile.close()
        GrpLenDctR = dict() # dict where keys are group names and values are corresponding group sizes - for groups with rnaseq samples
        GrpLenDctM = dict() # dict where keys are group names and values are corresponding group sizes - for groups with microarray samples
        for id, v in GrpDct.items():
            if id in RGrps:
                GrpLenDctR[id] = len(v)
            elif id in MGrps:
                GrpLenDctM[id] = len(v)
        AllGrpSmp = [smp for v in GrpDct.values() for smp in v] # list of all samples
        AllGrpSmpR = [smp for smp in AllGrpSmp if smp.startswith('R')] # list of all samples which are to be distributed among groups of rnaseq samples
        AllGrpSmpM = [smp for smp in AllGrpSmp if smp.startswith('M')] # list of all samples which are to be distrbuted among groups of microarray samples
        SimDistDfLst = list()
        InfDistDfLst = list()
        # create list of lists in which each list is samples ordered in a different way - to get dif. samples in same position in each simulation:
        SmpRandomOrdLstR = list()
        SmpRandomOrdLstM = list()
        for n in range(0, maxSim):
            lstR = AllGrpSmpR[:]
            np.random.shuffle(lstR) # shuffles sample order within list
            SmpRandomOrdLstR.append(lstR)
        for n in range(0, maxSim):
            lstM = AllGrpSmpM[:]
            np.random.shuffle(lstM) # shuffles sample order within list
            SmpRandomOrdLstM.append(lstM)
        # use multithreading:
        pool = Pool(processes=int(cpu_count() - (cpu_count()/2)))
        fc = partial(SimInf2ObservedDst, GrpLenDctR = GrpLenDctR, GrpLenDctM = GrpLenDctM, ThrDistDf=ThrDistDf, finalRcDMed=finalRcDMed)
        results = pool.map(fc, zip(SmpRandomOrdLstR, SmpRandomOrdLstM))
        pool.close()
        pool.join()
        for tr, dst in results: # tr = InfDistDf (True when sim dist < observed dist for each sample group); dst = SimThrDistDf (sim dist for each sample group)
            InfDistDfLst.append(tr)
            SimDistDfLst.append(dst)
        file1 = open(pathf1, 'wb')
        pickle.dump(SimDistDfLst, file1)
        file1.close()
        file2 = open(pathf2, 'wb')
        pickle.dump(InfDistDfLst, file2)
        file2.close()

def GetPvalueDfs(BaseDir):
    '''
    - get pvalues -> number of simulations where simulated distance was smaller than observed distance divided by total number of simulations
    :param BaseDir: basic directory
    :return PvalSimDstBellowObs: dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group and threshold combination
    '''
    # load list 1000 simulated dataframes where each dataframe has 'Trues' when distance of samples in a simulated group is lower than observed distance for same group:
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'SimInferior2RealDistanceDfList.pkl')
    infile = open(path,'rb')
    InfDistDfLst = pickle.load(infile)
    infile.close()
    # load dataframe with average observed distance between samples of same observed sample groups:
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'MeasuredDistanceDf.pkl')
    file = open(path, 'rb')
    ThrDistDf = pickle.load(file)
    file.close()
    # get dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group and threshold combination.
    # a.k.a. dataframe of p-values:
    SimStck = np.stack(InfDistDfLst)
    maxSim = SimStck.shape[0]
    PvalSimDstBellowObs = np.sum(SimStck, axis = 0)/(maxSim)
    PvalSimDstBellowObs = pd.DataFrame(PvalSimDstBellowObs, index=ThrDistDf.index, columns=ThrDistDf.columns)
    return PvalSimDstBellowObs

def getThrsLwstPval(PvalSimDstBellowObs, BaseDir):
    '''
    - save the thresholds with lowest until and including 3rd lowest p-value (percentage of simulations where average distance was lower than average observed distance)
      of sample groups for which we have cell lines with lethality score data (G1 and G17)
    :param PvalSimDstBellowObs: dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group and threshold combination
    :param BaseDir: basic directory
    '''
    IncOrd = PvalSimDstBellowObs.apply(lambda x: np.sort(x.unique()), axis=1)
    pos = 4 # excluded
    lowest = IncOrd.apply(lambda x: x[:pos])
    df=PvalSimDstBellowObs.T # to iterate over rows
    D = {x: list(df[x][df[x].isin(y)].index) for x,y in zip(df,lowest)} # dict with best thresholds (with lowest pvalues) for each group of samples
    # select just groups with samples for which we have cell lines with lethality score data (G1 and G17)
    # lethality scores are only available for CCs cell lines and some groups are composed of samples from donnors (not cell lines).
    # as we have lethality scores for a cell line of G1 (a microarray) and a cell line of G17 (rnaseq study) we test the correlation
    # with lethality scores for those studies and use the best thresholds for all remaining microarray and rnaseq studies, respectively.
    groups= ['G1', 'G17']
    Ds = {k:v for k,v in D.items() if k in groups}
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'thres2testGroups.pkl')
    saveResult(Ds, path)

def CreateBoxplot(data, hue, path):
    '''
    - creates formatted boxplot graph and saves
    :param data: dataframe with data
    :param hue: variable representing different colors in graph
    :param path: where to save
    :return:
    '''
    fig, ax = plt.subplots(figsize=(20/2,12/2))
    axM = sns.boxplot(x='Group', y='value', hue=hue, data=data)
    handles, labels = axM.get_legend_handles_labels()
    plt.legend(handles, labels, bbox_to_anchor=(1.015, 1), loc=2, borderaxespad=0., fontsize=15)
    axM.set_xlabel('', fontsize=20)
    axM.set_ylabel('', fontsize=20, labelpad=5)
    axM.tick_params(labelsize=15)
    fig.tight_layout()
    #plt.show()
    plt.savefig(path, dpi=300)
    plt.close('all')

def BoxplotTestThres(BaseDir, PvalSimDstBellowObs, MGrps, RGrps):
    '''
    - creates/saves boxplots of reaction scores to compare dif. threshold strategies
    :param BaseDir: basic directory
    :param PvalSimDstBellowObs: dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group and threshold combination
    :param MGrps: list of ids for groups with microarray samples
    :param RGrps: list of ids for groups with rnaseq samples
    '''
    BoxDf = PvalSimDstBellowObs.copy()
    BoxDf.columns = ['_'.join(t) for t in BoxDf.columns.to_list()]
    BoxDf['Group'] = BoxDf.index
    BoxDf = pd.melt(BoxDf, id_vars=['Group'])
    BoxDf.index = BoxDf['Group']
    BoxDf['Genes'] = BoxDf['variable'].apply(lambda x: x.split('_')[0])
    BoxDf['Threshold'] = BoxDf['variable'].apply(lambda x: '_'.join(x.split('_')[1:-1]))
    BoxDf['And/Or Rule'] = BoxDf['variable'].apply(lambda x: x.split('_')[-1])
    BoxDf.drop(columns='variable', inplace=True)
    BoxDf['Threshold Strategy'] = BoxDf['Threshold'].apply(lambda x: x.split('_')[0])
    BoxDfM = BoxDf.loc[MGrps]
    BoxDfR = BoxDf.loc[RGrps]
    # we later decided to exclude from the analysis samples of normal cells and normals stem cells,
    # as those were only available for one tissue:
    BoxDfR.drop(index=['G13', 'G14'], inplace=True)
    rd = {'G15':'G13','G16':'G14','G17':'G15','G18':'G16','G19':'G17','G20':'G18','G21':'G19','G22':'G20'}
    BoxDfR.rename(index=rd, inplace=True) # to not show an interval in numbers in graph
    BoxDfR.replace({'Group': rd}, inplace=True)
    ### Boxplot of pvalues of each group of samples for all genes vs only metabolic genes:
    # microarray:
    pathM = os.path.join(BaseDir, 'results/TestGeneScores', 'AllvsMetbGenesM.svg')
    CreateBoxplot(data=BoxDfM, hue='Genes', path=pathM)
    # rnaseq:
    pathR = os.path.join(BaseDir, 'results/TestGeneScores', 'AllvsMetbGenesR.svg')
    CreateBoxplot(data=BoxDfR, hue='Genes', path=pathR)
    ### Boxplot of pvalues of each group of samples for minsum vs minmax:
    # microarray:
    pathM = os.path.join(BaseDir, 'results/TestGeneScores', 'AndOrStrategiesM.svg')
    CreateBoxplot(data=BoxDfM, hue='And/Or Rule', path=pathM)
    # rnaseq:
    pathR = os.path.join(BaseDir, 'results/TestGeneScores', 'AndOrStrategiesR.svg')
    CreateBoxplot(data=BoxDfR, hue='And/Or Rule', path=pathR)
    ### Boxplot of pvalues of each group of samples for diff. threshold strategies:
    # microarray:
    pathM = os.path.join(BaseDir, 'results/TestGeneScores', 'ThrsM.svg')
    CreateBoxplot(data=BoxDfM, hue='Threshold Strategy', path=pathM)
    #rnaseq:
    pathR = os.path.join(BaseDir, 'results/TestGeneScores', 'ThrsR.svg')
    CreateBoxplot(data=BoxDfR, hue='Threshold Strategy', path=pathR)

def CreateAllSmpThrCmbDf(BaseDir):
    '''
    - creates dataframe with all samples + threshold combinations on rows and reactions on column,
    with reaction scores as values - where NAs are replaced by median
    :param Basedir: basic diretory
    '''
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct_medianCorr.pkl')
    infile = open(path, 'rb')
    finalRcDMed = pickle.load(infile)
    infile.close()
    thrDfLst = list()
    for k,v in finalRcDMed.items():
        colnms = [smp + '.' + '.'.join(k) for smp in v.columns] # list of samplenames+threshold comb
        v.columns = colnms # rename columns so that we can join dataframes of dif thresholds below without having duplicate column names
        # add column to dataframe with threshold combination of that dataframe:
        thrR = pd.DataFrame([k[0]] * v.shape[1], columns=['Threshold'], index = v.columns).T
        # add column to dataframe with And/Or rules combination of that dataframe:
        aofR = pd.DataFrame([k[1]] * v.shape[1], columns=['And_OR_rule'], index = v.columns).T
        v = v.append(thrR)
        v = v.append(aofR)
        # join dataframe of one threshold+and/or rule to list of dataframes:
        thrDfLst.append(v)
    # join list of dataframes into a mega dataframe with all samples of all threshold+And/Or rule comb:
    AlThrSmp = reduce(lambda x,y: x.join(y,how='outer'),thrDfLst).T
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresRcScoreDf_median.pkl')
    file = open(path, 'wb')
    pickle.dump(AlThrSmp, file)
    file.close()

def createPCAgraph(valuesDf, completeDf, classes, classesName, where2Save, legend=False, Ncomp=2, comp2plot=['PC1', 'PC2']):
    '''
    - creates formatted PCA plots and saves
    :param valuesDf:
    :param completeDf:
    :param classes: list of names of classes to split data (used in legend of scatter plot), i.e. local1, local2, etc
    :param classesName: name of column in 'completeDf' containing 'classes' for all observations
    :param where2Save: path to save plot
    :param legend: boolean indicating whether to add a legend or not
    :param Ncomp: number of components o test
    :param comp2plot: list with components to present in PCA axis, default is ['PC1', 'PC2']
    '''
    # scale features/columns/reactions to normally distributed values with a mean of 0 and std of 1 (centered to have zero mean):
    x = StandardScaler().fit_transform(valuesDf)
    # get dataframe with pc (principal component) values for all samples:
    pca = PCA(n_components=Ncomp) # pca with 2 PCs
    pc_best = pca.fit_transform(x)
    compCol = list(map(lambda x: 'PC' + str(x), range(1,Ncomp+1)))
    pcDf = pd.DataFrame(data = pc_best, columns = compCol, index = completeDf.index)
    # get amount of variance each pc holds:
    pcVar = pca.explained_variance_ratio_
    # graph:
    plt.figure(figsize=(10/2,10/2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('PC1({}%)'.format(round(pcVar[0]*100,2)),fontsize=15)
    plt.ylabel('PC2({}%)'.format(round(pcVar[1]*100,2)),fontsize=15)
    fake = Factory.create()
    clLst = [fake.hex_color() for n in range(0,len(classes))]
    for g,c in zip(classes,clLst):
        keepInd = completeDf[classesName] == g
        plt.scatter(pcDf.loc[keepInd, comp2plot[0]], pcDf.loc[keepInd, comp2plot[1]], color=c, s=50)
    if legend:
        plt.legend(classes, prop={'size': 10})
    plt.tight_layout()
    #plt.show()
    plt.savefig(where2Save, dpi=300)
    plt.close('all')

def PCATestThres(BaseDir):
    '''
    :param BaseDir: basic directory
    - creates/saves PCA plots of reaction scores to compare dif. threshold strategies
    '''
    # create dataframe with all samples + threshold combinations on rows and reactions on column, with reaction scores as values - where NAs are replaced by median:
    if not os.path.exists(os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresRcScoreDf_median.pkl')):
        CreateAllSmpThrCmbDf(BaseDir)
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresRcScoreDf_median.pkl')
    infile = open(path,'rb')
    AlThrSmp = pickle.load(infile)
    infile.close()
    # prepare dataframe for PCA graphs:
    FAlThrSmp = AlThrSmp.copy()
    FAlThrSmp['Threshold Strategy'] = FAlThrSmp['Threshold'].apply(lambda y: y.split('_')[1]) # add column with threshold strategy
    FAlThrSmp['Genes'] = FAlThrSmp['Threshold'].apply(lambda y: y.split('_')[0])
    thresholds = list(FAlThrSmp['Threshold Strategy'].unique())
    aor = list(FAlThrSmp['And_OR_rule'].unique())
    genes = list(FAlThrSmp['Genes'].unique())
    ## PCA of threshold strategies (global, local1 and local2):
    # rnaseq:
    FAlThrSmpR = FAlThrSmp[list(map(lambda x: x.startswith('R'), FAlThrSmp.index))]
    FAlThrSmpRValues = FAlThrSmpR.iloc[:,:- 4].values # from previous dataframe get only reaction scores as numpy array
    path=os.path.join(BaseDir, 'results/TestGeneScores','PCAthresR_PC1PC2.svg')
    createPCAgraph(valuesDf=FAlThrSmpRValues, completeDf=FAlThrSmpR, classes= thresholds, classesName='Threshold Strategy', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
    # microarray:
    FAlThrSmpM = FAlThrSmp[list(map(lambda x: x.startswith('M'), FAlThrSmp.index))]
    FAlThrSmpMValues = FAlThrSmpM.iloc[:,:- 4].values
    path=os.path.join(BaseDir, 'results/TestGeneScores','PCAthresM_PC1PC2.svg')
    createPCAgraph(valuesDf=FAlThrSmpMValues, completeDf=FAlThrSmpM, classes= thresholds, classesName='Threshold Strategy', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
    ## PCA of AND/OR rules:
    # rnaseq:
    path = os.path.join(BaseDir, 'results/TestGeneScores','PCAandorR_PC1PC2.svg')
    createPCAgraph(valuesDf=FAlThrSmpRValues, completeDf=FAlThrSmpR, classes= aor, classesName='And_OR_rule', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
    # microarray:
    path = os.path.join(BaseDir, 'results/TestGeneScores', 'PCAandorM_PC1PC2.svg')
    createPCAgraph(valuesDf=FAlThrSmpMValues, completeDf=FAlThrSmpM, classes= aor, classesName='And_OR_rule', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
    ## PCA of AllGenes/MetbGenes:
    # rnaseq:
    path = os.path.join(BaseDir, 'results/TestGeneScores','PCAgenesR_PC1PC2.svg')
    createPCAgraph(valuesDf=FAlThrSmpRValues, completeDf=FAlThrSmpR, classes= genes, classesName='Genes', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])
    # microarray:
    path = os.path.join(BaseDir, 'results/TestGeneScores','PCAgenesM_PC1PC2.svg')
    createPCAgraph(valuesDf=FAlThrSmpMValues, completeDf=FAlThrSmpM, classes= genes, classesName='Genes', where2Save=path, legend=True, Ncomp=2, comp2plot=['PC1', 'PC2'])