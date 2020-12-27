### Modules:
import os
import pandas as pd
import seaborn as sns
import numpy as np
import pickle
import itertools
from functools import reduce
from scipy.spatial import distance
from troppo.omics.core import OmicsContainer
from cobamp.utilities.parallel import batch_run
from ModelInfo import preprocessiHuman
from sklearn.preprocessing import MinMaxScaler
from troppo.methods_wrappers import ReconstructionWrapper
from multiprocessing import cpu_count
from troppo.omics.core import IdentifierMapping, TypedOmicsMeasurementSet
from itertools import product, chain
from multiprocessing.dummy import Pool
from functools import partial

### Functions:

def scaleGenExp(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, HamsPath, BaseMediumPath, MediumPath, BiomassID, BaseDir):
    '''
    - applies min-max normalization to individual studies' gene expression data (to all genes/ just genes in generic metabolic model)
      saves to .csv files and creates clustermaps
    - returns generic model and model adapted for medium composition
    :param modelPathSBML: param of preprocessiHuman function
    :param modelProcPath: param of preprocessiHuman function
    :param modelProcInfoPath: param of preprocessiHuman function
    :param modelMedAdapPath: param of preprocessiHuman function
    :param modelInfoPath: param of preprocessiHuman function
    :param HamsPath: param of preprocessiHuman function
    :param BaseMediumPath: param of preprocessiHuman function
    :param MediumPath: param of preprocessiHuman function
    :param BiomassID: param of preprocessiHuman function
    :param BaseDir: base directory
    :return model: generic metabolic model
    :return modelMedAdap: metabolic model adapted for medium composition
    '''
    # get list with all genes in medium adapted metabolic model, human1:
    model, modelMedAdap = preprocessiHuman(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, HamsPath, BaseMediumPath, MediumPath, BiomassID)
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
        studyExpVal = studyExpVal.iloc[:,[0] + list(range(3,studyExpVal.shape[1]))]
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
        # clustermap of genes of human model with min-max normalization:
        ncMtb = fnMtb.shape[1]
        nrMtb = fnMtb.shape[0]
        fclMtb = 10/ncMtb
        sns.set(font_scale=fclMtb)
        resMtb = sns.clustermap(fnMtb,  cmap='coolwarm', linecolor='black', xticklabels=True, yticklabels=False, linewidths=0.00001, figsize=(0.77*ncMtb, 0.0065*nrMtb), dendrogram_ratio=(0.2,0.05))
        resMtb.fig.subplots_adjust(bottom=0.27 * fclMtb)
        resMtb.savefig(os.path.join(ResStudyDir, 'explCluster_MetbGenes.png'), dpi=300)
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
    GrpDf = pd.read_excel(os.path.join(BaseDir, 'data/expectedSampleGroups.xlsx'), sheet_name='groups_by_study_cell_type', index=False)
    GrpDct = {el : list(GrpDf[el].dropna()) for el in GrpDf.columns}
    ThrDistDf = FindEucDist(GrpDct, finalRcDMed)
    ThrDistDf = ThrDistDf.loc[GrpDf.columns,list(finalRcDMed.keys())] # to order rows and columns
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'MeasuredDistanceDf.pkl')
    saveResult(res=ThrDistDf, path=path)
    return GrpDct

def SimInf2ObservedDst(AllGrpSmp, GrpLenDct, ThrDistDf, finalRcDMed):
    '''
    - to be used by function SimEucDistWithinSampleGroups
    - for one random simulation gets dataframe with 'True' values when average simulated distances are inferior to observed distances of corresponding group
    and gets average of distance between samples of same simulated groups
    :param AllGrpSmp: one list of randomly ordered samples
    :param GrpLenDct: dict with group ids as keys and group length as values
    :param ThrDistDf: dataframe with thresholds as columns and groups as rows with averaged observed distances (between reaction scores) for each group as values
    :param finalRcDMed: dict with thresholds as keys and dataframes with reaction scores as values
    :return InfDistDf: dataframe with 'True' values when average simulated distances are inferior to observed distances of corresponding group
    :return SimThrDistDf: dataframe with thresholds as columns and groups as rows, where values are average of distance between samples of same simulated group
    '''
    SimGrpDct = dict()  # get dict with group number as key and as values we have a list of specific length of random sample names
    i = 0
    for id, v in GrpLenDct.items():
        f = i + v
        SimGrpDct[id] = AllGrpSmp[i:f]
        i = i + v
    # create a dataframe with thresholds as columns and groups as rows, where values are average of distance between samples of same simulated group:
    SimThrDistDf = FindEucDist(SimGrpDct, finalRcDMed)
    SimThrDistDf = SimThrDistDf.loc[ThrDistDf.index, ThrDistDf.columns]  # order rows/columns of dataframe of simulated distances according to dataframe of real distances - to latter on stack them with matching positions
    # create dataframe with 'True' values when average simulated distances are inferior to observed distances of corresponding group:
    InfDistDf = SimThrDistDf < ThrDistDf
    return InfDistDf, SimThrDistDf

def SimEucDistWithinSampleGroups(BaseDir, GrpDct, maxSim = 1000):
    '''
    - save list of dataframes where each dataframe represents a simulation and 'True' values exist when average simulated distances for a group are inferior to observed distances for that group
    - save list of dataframes where each dataframe represents a simulation and values are average of distance between samples of same simulated group
    :param BaseDir: base directory
    :param GrpDct: dict with group id as keys and corresponding donors/cell lines as values
    :param maxSim: number of simulations
    '''
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'MeasuredDistanceDf.pkl') # dataframe with averaged observed distances in each group and threshold
    infile = open(path,'rb')
    ThrDistDf = pickle.load(infile)
    infile.close()
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct_medianCorr.pkl') # dict with thresholds as keys and dataframe with react scores as values, where NAs where replaced by median gene expression
    infile = open(path,'rb')
    finalRcDMed = pickle.load(infile)
    infile.close()
    GrpLenDct = {id : len(v) for id,v in GrpDct.items()} # dict where keys are group names and values are corresponding group sizes
    AllGrpSmp = [smp for v in GrpDct.values() for smp in v] # list of all samples which are to be distributed among groups
    SimDistDfLst = list()
    InfDistDfLst = list()
    # create list of lists in which each list is samples ordered in a different way - to get dif. samples in same position in each simulation:
    SmpRandomOrdLst = list()
    for n in range(0, maxSim):
        np.random.shuffle(AllGrpSmp) # shuffles sample order within list
        SmpRandomOrdLst.append(AllGrpSmp)
    # use multithreading:
    pool = Pool(processes=int(cpu_count() - (cpu_count()/4)))
    fc = partial(SimInf2ObservedDst, GrpLenDct = GrpLenDct, ThrDistDf=ThrDistDf, finalRcDMed=finalRcDMed)
    results = pool.map(fc, SmpRandomOrdLst)
    pool.close()
    pool.join()
    for tr, dst in results: # tr = InfDistDf (True when sim dist < observed dist for each sample group); dst = SimThrDistDf (sim dist for each sample group)
        InfDistDfLst.append(tr)
        SimDistDfLst.append(dst)
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'SimulatedDistanceDfList.pkl')
    file = open(path, 'wb')
    pickle.dump(SimDistDfLst, file)
    file.close()
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'SimInferior2RealDistanceDfList.pkl')
    file = open(path, 'wb')
    pickle.dump(InfDistDfLst, file)
    file.close()

def GetPvalueDfs(BaseDir, usrVal=1E-3):
    '''
    - get pvalues -> number of simulations where simulated distance was smaller than observed distance divided by total number of simulations
    - get dataframe with 'Trues' when a group pvalue is below a user defined threshold
    :param BaseDir: basic directory
    :param usrVal: user defined threshold to select groups where pvalue is bellow it
    :return PvalSimDstBellowObs: dataframe with number of simulations where simulated distance was smaller than observed distance divided by total number of simulations, for each sample group and threshold combination
    :return IsBellow: dataframe with 'Trues' when a group pvalue is below a user defined threshold. rows are groups and columns are react scores' thresholds
    :return usrVal: userval chosen above
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
    # get dataframe with 'Trues' when a group pvalue is below a user defined threshold.
    # rows are groups and columns are thresholds + And/Or rule combinations:
    IsBellow = pd.DataFrame(PvalSimDstBellowObs < usrVal, index=ThrDistDf.index, columns=ThrDistDf.columns)
    return PvalSimDstBellowObs, IsBellow, usrVal

def BestHighestNbGroups(IsBellow, BaseDir):
    '''
    - saves best performing thresholds+and/or rules for all groups of samples (the ones with 3 highest number of groups)
    :param IsBellow: dataframe with 'Trues' when a group pvalue is below a user defined threshold. rows are groups and columns are gene scores' thresholds
    :param BaseDir: basic directory
    '''
    #larg = IsBellow.sum().nlargest(2, keep='all')
    #larg = IsBellow.sum()[IsBellow.sum() == IsBellow.sum().max()]
    ###
    Tv = -np.sort(-IsBellow.sum().unique())
    rank = 3
    t = Tv[:rank].tolist()
    Blarg = IsBellow.sum()[IsBellow.sum().isin(t)]
    #Blarg = pd.DataFrame([i for i in Blarg.index if i[1] == 'minsum']).set_index(0)[1]  # select only minsum
    Blarg.to_csv(os.path.join(BaseDir, 'support/TestGeneScores', 'threshighestNbGroups'), sep='\t', index_label=[0,1], header=False)

'''
Mv = -np.sort(-IsBellow.loc[MGrps].sum().unique()) # sort sums
Rv = -np.sort(-IsBellow.loc[RGrps].sum().unique())
Tv = -np.sort(-IsBellow.sum().unique())
rank = 1
tM, tR ,t = Mv[:rank].tolist(), Rv[:rank].tolist(), Tv[:rank].tolist()
largM = IsBellow.loc[MGrps].sum()[IsBellow.loc[MGrps].sum().isin(tM)]
BlargM = [i for i in largM.index if i[1] == 'minsum']
BlargM
len(BlargM)
largR = IsBellow.loc[RGrps].sum()[IsBellow.loc[RGrps].sum().isin(tR)]
BlargR = [i for i in largR.index if i[1] == 'minsum']
BlargR
len(BlargR)
larg = IsBellow.sum()[IsBellow.sum().isin(t)]
Blarg = [i for i in larg.index if i[1] == 'minsum']
Blarg
len(Blarg)

int = set(BlargM) & set(BlargR)
len(int)

'''

