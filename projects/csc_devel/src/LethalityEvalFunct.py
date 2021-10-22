### Modules:
import os
import pandas as pd
import pickle
import ast
import re
import numpy as np
from sklearn.metrics import matthews_corrcoef
from functools import reduce
from troppo.methods_wrappers import ReconstructionWrapper, GapfillWrapper
from multiprocessing import cpu_count
from cobamp.utilities.parallel import batch_run
from troppo.tasks.core import TaskEvaluator
from tasksFunct import AddTaskCobrapy, TaskGapfill, EvalAllTasksAtSameTime
from cobamp.core.optimization import BatchOptimizer
from cobamp.wrappers.external_wrappers import get_model_reader
from troppo.omics.core import IdentifierMapping, TypedOmicsMeasurementSet
from cobamp.utilities.parallel import MP_THREADS
from retrying import retry

### Functions:
def saveResult(res, path):  # save results as pickle
    file = open(path, 'wb')
    pickle.dump(res, file)
    file.close()

def ThrsRcScores4aSample(BaseDir, model):
    '''
    - save a list of threshold + and/or rule names and list of dictionaries with corresponding reaction:reaction-scores,
      for best thresholds to test in G1 and G17
    :param BaseDir: base directory
    :param model: generic model
    '''
    crspD = {'G1': 'M2.2.13.CCs_Huh7', 'G17': 'R5.8.2.TN.CCs_HCC1937'}
    # get dict with threshold + and/or rule combinations vs dataframe with samples in columns and reactions in rows, where values are react scores:
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'AllStudiesAllthresholdsDct.pkl')
    infile = open(path,'rb')
    StdReacScores = pickle.load(infile)
    infile.close()
    # load all thresholds to test for each study/sample group
    path = os.path.join(BaseDir, 'support/TestGeneScores', 'thres2testGroups.pkl')
    infile = open(path, 'rb')
    thrsD = pickle.load(infile)
    infile.close()
    for k,v in thrsD.items(): # for each group of samples/study
        # get a list of thresholds to test and corresponding dicts with reaction scores of sample for which we want to reconstruct models:
        rcLst = list()
        thrLst = list()
        for thr in v:
            rcD = StdReacScores[thr][crspD[k]].fillna(0).to_dict() # replace 'NaN' (corresponding to genes present in one microarray/rnaseq platform that are absent in another) by 0 (it's the score that is given to 'None' values bellow)
            rcLst.append(rcD)
            thrLst.append(thr)
        # add scores 'None' to reactions that do not have gene-reaction-rule (cause those had been excluded from 'AllStudiesAllthresholdsDct.pkl' in function GetDictOfReacScoresDF of testGeneScrs.py):
        mdRcIds = {r.id for r in model.reactions}
        WithGeneReacRule = set(rcLst[0].keys())
        NoGeneReacRule = mdRcIds - WithGeneReacRule
        NoGeneReacRuleDct = {rid: None for rid in NoGeneReacRule}
        for d in rcLst:
            d.update(NoGeneReacRuleDct)
        thrdir = os.path.join(BaseDir, 'support/TestGeneScores/models2Test')
        rcdir = os.path.join(BaseDir, 'support/TestGeneScores/rcScores2Test')
        if not os.path.exists(thrdir):
            os.makedirs(thrdir)
        if not os.path.exists(rcdir):
            os.makedirs(rcdir)
        saveResult(res=thrLst, path=os.path.join(thrdir, crspD[k] + '.pkl')) # save list of thresholds to test
        saveResult(res=rcLst, path=os.path.join(rcdir, crspD[k] + '.pkl')) # save list of corresponding reaction scores to test
        pd.DataFrame([len(thrLst)]).to_csv(os.path.join(BaseDir, 'support/TestGeneScores/models2Test/NumbModels' + crspD[k] + '.tab'), sep='\t', header=False, index=False)


def createLstRct2Test(model, BaseDir):
    '''
    - split model reaction ids (except boundaries/exchanges) into different files to be used for analysis with cluster
      (cluster uses different jobs, ech job processes one group of reactions. 1 file)
    :param model: generic metabolic model
    :param BaseDir: basic directory
    '''
    # close model exchange/boundary reactions (in "human1" all boundaries are exchange reactions):
    taskModel = model.copy()
    for r in taskModel.boundary:
        r.knock_out()
    # get all reactions of generic model except boundaries/exchange reactions:
    rids = {r.id for r in taskModel.reactions}
    bdrids = {r.id for r in taskModel.boundary}
    frids = list(rids - bdrids)
    # split reactions in groups (each group is going to be run as a job in the cluster):
    grps = 30
    integ = len(frids)//grps # how many elements each group is going to have
    remain = len(frids)%grps # how many elements the last group has
    intv = grps*[integ] + [remain]
    i = 0
    dir = os.path.join(BaseDir, 'support/LethalityEval/ModelRc2Test')
    if not os.path.exists(dir):
        os.makedirs(dir)
    for v in intv: # v is size of each group
        print(i)
        f = i + v
        pd.Series(frids[i:f]).to_csv(os.path.join(BaseDir, 'support/LethalityEval/ModelRc2Test/EssLst' + str(i)), index=False, header=False) # saves each group in a file
        i = i + v


def fastcoreReacfunc_ReacScores(Rscore, params):
    '''
    - returns series with booleans. 'True' when reaction is to be kept in reconstructed model according to FASTCORE algorithm, 'False' when to be excluded
      reactions associated with biomass production are kept (protected) and 'None' values' reactions are excluded.
    :param Rscore: dictionary for one model with reaction ids as keys and reaction scores as values
    :param params: dictionary with key 'rw' and a reconstruction wrapper of generic model as value
    '''
    rw = params['rw']  # load parameters
    t = 0 # reaction scores above threshold, are reaction scores above 0
    protected = ['biomass_human', 'HMR_10023', 'HMR_10024']
    # for each reaction if score is not None or above threshold or in 'Protected' (like biomass id) we get the score:
    try:
        scores = [[k for k, v in Rscore.items() if (v is not None and v > t) or k in protected]]
        return rw.run_from_omics(omics_data=scores, integration_strategy='', algorithm='fastcore', solver='CPLEX')
    except Exception as e:
        print(e)
        return {r: False for r in rw.model_reader.r_ids}

def initReacfunc_ReacScores(Rscore, params):
    '''
    - returns series with booleans. 'True' when reaction is to be kept in reconstructed model according to INIT algorithm, 'False' when to be excluded
      reactions associated with biomass production are kept (protected - are given highest score) and 'None' reactions are given score of 0.
    :param Rscore: dictionary for one model with reaction ids as keys and reaction scores as values
    :param params: dictionary with key 'rw' and a reconstruction wrapper of generic model as value
    '''
    rw = params['rw']  # load parameters
    protected = ['biomass_human', 'HMR_10023', 'HMR_10024']
    try:
        scores = {k:(v/2 if v < 0 else v) if v is not None else 0 for k, v in Rscore.items()}
        scores.update({x: max(scores.values()) for x in protected})
        return rw.run_from_omics(omics_data=scores, integration_strategy='', algorithm='tinit', solver='CPLEX')
    except Exception as e:
        print(e)
        return {r: False for r in rw.model_reader.r_ids}

def ReconstructionReac_ReacScores(model, labs, iters, path, alg, testThres):
    '''
    - get/save dataframe where columns are reactions and rows are tested threshold+and/or rule combinations. 'True' value means reaction is to be kept in reconstructed model.
    :param model: generic metabolic model unconstrained (with exchange reactions open) - cause if your core set contains reactions that are incompatible with the media, the algorithm won't be able to reconstruct the model.
    :param labs: list of strings with id of threshold+and/or rule combination/ model name to test
    :param iters: list of dictionaries with reaction scores, each element in this list corresponds to one element in labs list
    :param path: path to file that will save results
    :param alg: string can be 'both' (test both init and fastcore), 'init' or 'fastcore'
    :param testThres: use the function to test different thresholds (True) or to reconstruct different models with same threshold (False)
    '''
    # Create a reconstruction wrapper instance that loads pre-processed metabolic model:
    # when finding reactions needed to reconstruct models, you should use it unconstrained with open exchanges - cause if your core set contains reactions that are incompatible with the media, the algorithm won't be able to reconstruct the model.
    rw = ReconstructionWrapper(model, ttg_ratio=9999)
    NTHREADS = int(cpu_count() - (cpu_count()/4))
    result_dicts = {}
    if testThres:
        flabs = [tuple(['fastcore'] + list(l)) for l in labs]  # labels of models and and/or rules for fastcore
        tlabs = [tuple(['init'] + list(l)) for l in labs]  # labels of models and and/or rules for INIT
    else:
        flabs = [tuple(['fastcore', l]) for l in labs]  # labels of models and and/or rules for fastcore
        tlabs = [tuple(['init', l]) for l in labs]  # labels of models and and/or rules for INIT
    ## test both algorithmns by default, otherwise test just the user defined algorithm:
    if alg == 'both':
        # Get reactions to exclude when reconstruting with FASTCORE:
        output = batch_run(fastcoreReacfunc_ReacScores, iters, {'rw': rw}, threads=min(len(iters), NTHREADS))
        batch_fastcore_res = dict(zip(flabs, output)) # create dictionary with fastcore combinations' labels and dictionaries mapping reactionId to bolean (True - means reaction kept after reconstruction)
        result_dicts.update(batch_fastcore_res) # to merge above dictionary with a similar one below for tinit
        # Get reactions to exclude when reconstructing with INIT:
        output = batch_run(initReacfunc_ReacScores, iters, {'rw': rw}, threads=min(len(iters), 2))
        batch_tinit_res = dict(zip(tlabs, output))
        result_dicts.update(batch_tinit_res)
    elif alg == 'fastcore':
        # Get reactions to exclude when reconstructing with FASTCORE:
        output = batch_run(fastcoreReacfunc_ReacScores, iters, {'rw': rw}, threads=min(len(iters), NTHREADS))
        batch_fastcore_res = dict(zip(flabs,output))  # create dictionary with fastcore combinations' labels and dictionaries mapping reactionId to bolean (True - means reaction kept after reconstruction)
        result_dicts.update(batch_fastcore_res)  # to merge above dictionary with a similar one below for tinit
    elif alg == 'init':
        # Get reactions to exclude when reconstructing with INIT:
        output = batch_run(initReacfunc_ReacScores, iters, {'rw': rw}, threads=min(len(iters), 2))
        batch_tinit_res = dict(zip(tlabs, output))
        result_dicts.update(batch_tinit_res)
    # Create dataframe with results from all models:
    reacRes = pd.DataFrame.from_dict(result_dicts, orient='index')
    reacRes.to_csv(path)

def retry_if_RuntimeError(exception):
    """Return True if we should retry (in this case when it's an RuntimeError), False otherwise"""
    return isinstance(exception, RuntimeError)

def Reconst_taskEval(modelMedAdap, task_list, BaseDir, BiomassMetbExtID, protected, Finalpath, BoolDict, StudyNumb, fname):
    '''
    - gapfill reconstructed models with EFM methods and then gapfill for tasks with cobrapy gapfill,
    - saves: * info on models that did not grow or did not performed tasks and info on gafill solution reactions
             * dataframe where rows are reactionIDs, columns are models' names, values are bounds of reconstructed and gapfilled models on specific medium
    :param modelMedAdap: model adapted for specific medium composition (exchange reaction of medium metabolites are open and remaining have bounds LB=0, UB=1000)
    :param task_list: list of tasks to test
    :param BaseDir: basic directory
    :param BiomassMetbExtID: id of biomass metabolite in external compartment
    :param protected: reactions related with biomass production and exchange/transport to in/out of cell
    :param Finalpath: directory where to save results
    :param BoolDict: dict where key is ('algorithm', 'threshold', 'And/Or function') and value another dict with (reactID: True/False - to include or to exclude)
    :param StudyNumb: when testing different thresholds for same study is a string with study number, otherwise use ''
    :param fname: name of file with reconst. alg. output
    '''
    exg_rx_ids = {r.id for r in modelMedAdap.exchanges}  # model exchange reactions ids
    tasksDic = {task.name: task for task in task_list}  # dict with task name : task
    gpfl_wrapper = GapfillWrapper(modelMedAdap)
    NotReconst = list()
    NotProdBmss = list()
    FailedTasks = dict()
    RecModelGapReac = dict()
    GapFillTasksReac = dict()
    RecModelReacBds = dict()
    EssRcTasksPath = os.path.join(BaseDir, 'support/MetbTasks/essentialTasks/EssentialReactions4Tasks.tab')
    EssRcTasks = pd.read_csv(EssRcTasksPath, sep= '\t').iloc[:,0]
    for modelID, d in BoolDict.items():
        print(modelID)
        # if algorithm return 'False' for all reactions in model (it didn't find a solution), then continue to next model:
        if sum(d.values()) == 0:  # if algorithm return 'False' for all reactions in model (it didn't find a solution), then continue to next model
            print('algorithm did not find solution for model ' + str(modelID))
            NotReconst.append('_'.join(modelID))
            continue
        # list of reactions to be removed according to reconstruction algorithm,
        # excludes exchange reactions (to preserve adaptation to medium composition),
        # excludes reactions essential for essential tasks in generic model and reactions related with biomass production/transport:
        Id2Remove = list({id for id, bol in d.items() if not bol} - exg_rx_ids - set(protected) - set(EssRcTasks))
        # reconstruct model:
        taskModel = modelMedAdap.copy()  # model used as reconstructed model and for task gapfill - adapted for medium composition
        taskModel.remove_reactions(Id2Remove, remove_orphans=False) # do not remove orphan metabolites (that have no association with any reaction) for now, to be able to test tasks that use these metabolites
        ## evaluate/gapfill each essential task on reconstructed model (sequential gapfill):
        GapFillTasksReacIds = list()
        UnivTaskModel = modelMedAdap.copy()  # generic model to use for model/task gapfill - adapted for medium composition
        for k in taskModel.boundary:  # close reconstructed model boundaries to test tasks
            k.knock_out()
        for k in UnivTaskModel.boundary: # close generic model boundaries to test tasks
            k.knock_out()
        for tsk in task_list:
            print(tsk)
            tasksEval = TaskEvaluator(model=taskModel, tasks=task_list, solver='CPLEX')  # create a task evaluator instance. task gapfill is sequential so this need to be included in task for loop, to update task model evaluator
            tn = tsk.name # get task name
            tasksEval.current_task = tn # set task to evaluate
            taskres = tasksEval.evaluate()[0] # evaluate task
            tasksEval.current_task = None  # remove current task from evaluation - back to normal
            if not taskres: # if task fails, gapfill
                # add a task to generic model copy (UnivTaskModel) and reconstructed model copy (taskModel):
                UnivTaskModel, taskModel, taskReacIds = AddTaskCobrapy(tasksDic=tasksDic, UnivModel=UnivTaskModel, IncompModel=taskModel, Tname=tn)
                # find gapfill solution:
                # if task gapfill fails due to 'RuntimeError: failed to validate gapfilled model, try lowering the integer_threshold',
                # try again untill it works and up to max of 3 times
                @retry(retry_on_exception=retry_if_RuntimeError, stop_max_attempt_number=3)
                def retry_3times():
                    print('b')
                    gapfillTasksolIds = TaskGapfill(UnivTaskModel, taskModel)  # list of reaction ids needed in reconstructed model for task completion
                    return gapfillTasksolIds
                try:
                    gapfillTasksolIds = retry_3times()
                except:
                    # remove task from generic model copy (UnivTaskModel) and taskModel:
                    taskModel.remove_reactions(taskReacIds, remove_orphans=False)  # do not remove orphan metabolites (that have no association with any reaction) for now, to be able to test tasks that use these metabolites
                    UnivTaskModel.remove_reactions(taskReacIds, remove_orphans=False)  # do not remove orphan metabolites (that have no association with any reaction) for now, to be able to test tasks that use these metabolites
                    continue
                else:
                    # save task gapfill solution to list with reactions' ids of all tasks' gapfill solutions:
                    GapFillTasksReacIds.append(gapfillTasksolIds)
                    # add task gapfill solution to taskModel:
                    reacTaskCopy = [UnivTaskModel.reactions.get_by_id(id).copy() for id in gapfillTasksolIds]
                    taskModel.add_reactions(reacTaskCopy)
                    # remove task from generic model copy (UnivTaskModel) and taskModel:
                    taskModel.remove_reactions(taskReacIds, remove_orphans=False)  # do not remove orphan metabolites (that have no association with any reaction) for now, to be able to test tasks that use these metabolites
                    UnivTaskModel.remove_reactions(taskReacIds, remove_orphans=False)  # do not remove orphan metabolites (that have no association with any reaction) for now, to be able to test tasks that use these metabolites
        # get gapfill solution of all tasks:
        if len(GapFillTasksReacIds) != 0:
            GapFillTasksReacIds = reduce(lambda x,y: x+y, GapFillTasksReacIds) # flats list of lists into a list
        GapFillTasksReac['_'.join(modelID)] = GapFillTasksReacIds
        # try to gapfill model with EFMs when model doesn't grow:
        nonConsReac = {exc.id for exc in modelMedAdap.exchanges if exc.bounds[0] >= 0 and exc.bounds[1] > 0}  # exchange reactions where metabolites are non-consumed (lb>=0, ub>0)
        nonConsMetb = [list(modelMedAdap.reactions.get_by_id(r).metabolites)[0].id for r in nonConsReac]  # metabolites of exchange reactions where metabolites are non-consumed
        nonConsMetb = set(nonConsMetb) - {BiomassMetbExtID}  # biomass has to be excluded from non-consumed (although it is not consumed) cause it is excreted (lb>0,ub>0, so lb is dif. from 0)
        ls_override = {'produced': {BiomassMetbExtID}, 'non_consumed': nonConsMetb}  # produced is biomass boundary. it's same as optimize for biomass production reaction, as there is only one reaction that produces biomass in cytosol. No need to turn biomass boundary irreversible (biomass[s] -> ) as it is already done when adapting medium
        # list of reactions to be removed according to reconstruction algorithm,
        # excludes exchange reactions (to preserve adaptation to medium composition),
        # excludes reactions essential for essential tasks in generic model, reactions related with biomass production/transport and reactions of task gapfill solution:
        Id2Remove = list({id for id, bol in d.items() if not bol} - exg_rx_ids - set(protected) - set(EssRcTasks) - set(GapFillTasksReacIds))
        # find model gapfill solution with EFM method if it can/if needed:
        gapfillsol = gpfl_wrapper.run(avbl_fluxes=Id2Remove, ls_override=ls_override, algorithm='efm')  # model gapfill solution
        if len(gapfillsol) != 0:  # if gapfillsol is not empty (if models needs gapfill to grow)
            gapfillsol = gapfillsol[0]
            RecModelGapReac['_'.join(modelID)] = gapfillsol
        # finally, reconstruct model with all info:
        reconstModel = modelMedAdap.copy()  # model to use for reconstruction
        finalId2Remove = list(set(Id2Remove) - set(gapfillsol)) # final reactions 2 remove: already excludes boundaries, reac. for biomass, essential 4 tasks in generic model, model and task gapfill solution
        for rid in finalId2Remove:  # remove final list of reactions to reconstruct model
            reconstModel.reactions.get_by_id(rid).knock_out()
        # test if reconstructed model produces biomass after adding task gapfill reactions:
        opt = reconstModel.optimize()
        status = opt.status
        val = opt.objective_value
        if status != 'optimal' or val < 1E-9:
            print('model ' + '_'.join(modelID) + ' not producing biomass after tasks gapfill')
            NotProdBmss.append(', '.join(['_'.join(modelID), 'feasibility: ' + status, 'objective_value: ' + str(val)]))
            continue  # move to next model
        # test tasks again:
        for k in reconstModel.boundary:  # close reconstructed model boundaries to test tasks
            k.knock_out()
        tasksName, batch_res_tasks = EvalAllTasksAtSameTime(reconstModel, task_list)
        TasksFailedAfterGapfill = list()
        for Tname, res in zip(tasksName, batch_res_tasks):
            print(Tname, res)
            if not res[0]:
                print('model ' + '_'.join(modelID) + ' failed task ' + Tname + ' after gapfill')
                TasksFailedAfterGapfill.append(Tname)
        if len(TasksFailedAfterGapfill) != 0: # if any of the tasks failed:
            TasksFailedAfterGapfill = ','.join(TasksFailedAfterGapfill)
            FailedTasks['_'.join(modelID)] = [TasksFailedAfterGapfill]  # save failed tasks
            continue # move to next model if any essential task failed (may fail when task gapfill takes longer than 20 min - cause moves on without offering solution)
        # put boundaries back accordingly to media composition:
        for k in modelMedAdap.boundary:
            reconstModel.reactions.get_by_id(k.id).bounds = k.bounds
        # dictionary with samples as keys and as values recontructed models' final reaction bounds:
            RecModelReacBds['_'.join(modelID)] = {r.id: reconstModel.reactions.get_by_id(r.id).bounds for r in reconstModel.reactions}
    ftdir=os.path.join(Finalpath, StudyNumb, 'FailedTasks')
    tgpdir=os.path.join(Finalpath, StudyNumb, 'TasksGapFillSol')
    bdir=os.path.join(Finalpath, StudyNumb, 'RecModelBounds')
    nosoldir=os.path.join(Finalpath, StudyNumb, 'ModelNoAlgmSol')
    nogrowdir=os.path.join(Finalpath, StudyNumb, 'ModelNoGrowNoFeasibility')
    mdpdir=os.path.join(Finalpath, StudyNumb, 'ModelGapFillSol')
    if not os.path.exists(ftdir): # create dir if does not exist
        os.makedirs(ftdir)
    if not os.path.exists(tgpdir):
        os.makedirs(tgpdir)
    if not os.path.exists(bdir):
        os.makedirs(bdir)
    if not os.path.exists(nosoldir):
        os.makedirs(nosoldir)
    if not os.path.exists(nogrowdir):
        os.makedirs(nogrowdir)
    if not os.path.exists(mdpdir):
        os.makedirs(mdpdir)
    # dataframe with models that failed tasks and info on which tasks were failed for those models:
    pd.DataFrame.from_dict(FailedTasks).to_csv(os.path.join(ftdir, fname), sep='\t', index=False)
    # dataframe with reactions added during tasks' gapfill. when model did not require tasks gapfill, value is empty:
    pd.DataFrame({k: pd.Series(v) for k, v in GapFillTasksReac.items()}).to_csv(os.path.join(tgpdir, fname), sep='\t', index=False)
    # dataframe where rows are reactionIDs, columns are models' names, values are bounds of reconstructed and gapfilled models on specific medium:
    RecModelBounds = pd.DataFrame(RecModelReacBds)
    RecModelBounds.to_csv(os.path.join(bdir, fname), sep='\t')
    # dataframe with model names for which algorithm could not find a solution:
    pd.DataFrame(NotReconst).to_csv(os.path.join(nosoldir, fname), sep='\t', index=False)
    # dataframe with feasibility status and objective value of models that could not produce biomass/or not be feasible after gapfill:
    pd.DataFrame(NotProdBmss).to_csv(os.path.join(nogrowdir, fname), sep='\t', index=False)
    # dataframe with reactions added during model gapfill, for models that needed it:
    pd.DataFrame({k: pd.Series(v) for k, v in RecModelGapReac.items()}).to_csv(os.path.join(mdpdir, fname), sep='\t', index=False)

def RcnstThrsEval (BaseDir, StudyNumberR, StudyNumberM, model, i, path2MdNames, path2MdScores, pathDir, testThres, alg):
    '''
    - Reconstruct models for thresholds+and/or rules combinations to test and for studies
    :param BaseDir: basic directory
    :param model: generic model Not adapted for medium composition
    :param i: index to select one threshold/model to test from list of all those to test
    :param path2MdNames: path to model names list
    :param path2MdScores: path to model reaction scores list
    :param pathDir: path to directory where to save results
    :param testThres: boolean indicating whether to use function 'ReconstructionReac_ReacScores' to test thresholds or to build models for studies
    :param alg: string can be 'both' (test both init and fastcore), 'init' or 'fastcore', or '' (when we do Not want to test thresholds)
    '''
    # load a list of thresholds to test and corresponding dicts with (reaction:reaction scores) of the sample for which we want to reconstruct models:
    infile = open(path2MdNames, 'rb')
    thrLst = pickle.load(infile)
    infile.close()
    infile2 = open(path2MdScores, 'rb')
    rcLst = pickle.load(infile2)
    infile2.close()
    # select one threshold/model to test at a time (to be able to run one model in one job of the cluster):
    thr = [thrLst[i]]
    thrD = [rcLst[i]]
    if not testThres: # if we want to build models for studies (no threshold testing)
        if thr[0].startswith('R'):
            thrsDf = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumberR, 'CorrScores.tab'), sep='\t', index_col=0)
            besthres = thrsDf.loc[thrsDf['value'].idxmax()].drop(['value', 'biomass_cutoff', 'essential_cutoff'])
        else:
            thrsDf = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumberM, 'CorrScores.tab'), sep='\t', index_col=0)
            besthres = thrsDf.loc[thrsDf['value'].idxmax()].drop(['value', 'biomass_cutoff', 'essential_cutoff'])
        alg=besthres.algorithm # best algorithm for rnaseq studies
        pathRc = os.path.join(pathDir, thr[0])
    else:
        pathRc = os.path.join(pathDir, '_'.join(thr[0]))
    # get reactions to include/exclude from reconstructed models according to FASTCORE/INIT:
    if not os.path.exists(pathDir):
        os.makedirs(pathDir)
    ReconstructionReac_ReacScores(model, labs=thr, iters=thrD, path=pathRc, alg=alg, testThres=testThres)

def RcnstGapFill4ThrsEval (BaseDir, BiomassMetbExtID, StudyNumber, modelMedAdap, task_list, protected, Intpath, Finalpath, testThres, fname):
    '''
    - Reconstruct models for and gapfill for essential tasks and model growth
    :param BaseDir: basic directory
    :param BiomassMetbExtID: id of biomass metabolite in external compartment
    :param StudyNumber: study number
    :param modelMedAdap: model adapted for specific medium composition (exchange reaction of medium metabolites are open and remaining have bounds LB=0, UB=1000)
    :param task_list: list of tasks to test
    :param protected: reactions related with biomass production and exchange/transport to in/out of cell
    :param Intpath: path to directory constaining files. each file represents one sample/threshold and has "True/False" indicating if reactions are to be kept or not
    :param Finalpath: path to directory where to save results
    :param testThres: boolean indicating whether models to reconstruct are for threshold testing or not
    :param fname: name of file with reconst. alg. output
    '''
    # reconstruct models, gapfill reconstructed models with "EFM" method, and gapfill for tasks using cobrapy gapfill:
    if testThres: # if reconstructed models are for threshold testing
        BoolDict = pd.read_csv(Intpath, header=0, index_col=[0, 1, 2]).T.to_dict() # dict will have models/thresholds as keys and inner dictionary has reaction ids as keys and "true/false" as values
    else: # if reconstructed models are not for threshold testing
        BoolDict = pd.read_csv(Intpath, header=0, index_col=[0, 1]).T.to_dict()
    Reconst_taskEval(modelMedAdap, task_list, BaseDir, BiomassMetbExtID, protected, Finalpath, BoolDict, StudyNumb=StudyNumber, fname=fname)

def ProcessExpLethalityScoresDf(LethalCellLinePath, modelMedAdap):
    '''
    - preprocess dataframe with experimental lethality scores for genes, in dif. cell lines
    :param LethalCellLinePath: path to file with experimental lethality scores for different cell lines
    :param modelMedAdap: model adapted for specific medium composition (exchange reaction of medium metabolites are open and remaining have bounds LB=0, UB=1000)
    :return ach_mset: processed omics measurement set object of experimental lethality scores
    :return cobamp_model: modelMedAdap model in cobamp format
    :return model_genes: all genes in cobamp_model
    '''
    # get cobamp model from cobrapy model:
    cobamp_model = get_model_reader(modelMedAdap, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
    # all genes in cobamp model:
    model_genes = cobamp_model.gpr.get_genes()
    # get dataframe where rows are cell lines, columns are 'genes (symbol - entrez id)' and values are score centered on zero,
    # where negative values are lethal genes:
    ach_df = pd.read_csv(LethalCellLinePath, index_col=0)
    # convert entrez_id into string:
    hgnc_df = pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt', index_col=0, sep='\t') # gene id conversion dataframe
    hgnc_df['entrez_id'] = [str(int(x)) if not np.isnan(x) else np.nan for x in hgnc_df['entrez_id']]
    mapping = IdentifierMapping('human_transcriptomics', hgnc_df)
    # Create an omics measurement set object with cell lines' dataset components:
    ach_mset = TypedOmicsMeasurementSet(ach_df.index, ach_df.columns, ach_df.values, mapping)
    # keep entrez gene id only in column names:
    ensembl_patt = re.compile('\([0-9]*\)')
    ach_mset.column_names = [str(ensembl_patt.findall(k)[0].replace('(', '').replace(')', '')) for k in ach_mset.column_names]
    # convert entrez to ensembl:
    conv_dict = mapping.get_id_table(ach_mset.column_names, 'entrez_id').set_index('entrez_id')['ensembl_gene_id'].to_dict() # dict where achiles entrez gene ids are keys and corresponding ensembl gene ids are values
    ach_mset.column_names = [conv_dict[k] if k in conv_dict else np.nan for k in ach_mset.column_names] # keep only ensembl gene ids, when convertion from entrez didn't work it was put as nan - conversion
    # drops columns corresponding to genes not in metabolic model:
    ach_mset.drop(columns=ach_mset.data.columns[~ach_mset.data.columns.isin(model_genes)]) # '~' symbol turns Trues to Falses
    ach_mset.transform(lambda x: x.fillna(0)) # each nan value in achiles dataframe is replaced by 0
    return ach_mset, cobamp_model, model_genes

def InactReactFromGeneKo(cobamp_model, model_genes):
    '''
    - get reactions that become inactive when a gene is knockout
    :param cobamp_model: model adapted for specific medium composition in cobamp format
    :param model_genes: all genes of cobamp_model
    :return kos_to_test: list of combination of reactions that are inactive, depending on which gene is knockout
    :return reaction_states_inactive: dict where key is gene to knockout and value is id of reaction(s) that are inactive in universal model when that knockout happens
    '''
    # func bellow creates dict where all model genes are True except gene we want to knockout which is False:
    def get_state_dict_from_ko(ko):
        d = {k: True for k in cobamp_model.gpr.get_genes()}
        for k in ko:
            d[k] = False
        return d
    # create dict where key is gene of model and value is a dict
    # inner dict has gene as key and True as value, except for the knockout:
    state_dicts = {gko: get_state_dict_from_ko([gko]) for gko in model_genes}
    # create dict where key is a gene to knockout and value is a dict
    # inner dict has as key a reaction that is inactive or no GPR associated and value is 'False' or 'None'
    reaction_states = {gko: {i: j for i, j in {k: cobamp_model.gpr.eval_gpr(ind, state_dicts[gko]) for ind, k in enumerate(cobamp_model.reaction_names)}.items() if not j} for gko in state_dicts.keys()}
        # {k: cobamp_model.gpr.eval_gpr(ind, state_dicts[gko]) for ind, k in enumerate(cobamp_model.reaction_names)}
        # is a dict for a certain gene knockout where key is reaction and value is True/False/None (active/inactive/no GPR rule associated)
        # {i: j for i, j in {k: cobamp_model.gpr.eval_gpr(ind, state_dicts[gko]) for ind, k in enumerate(cobamp_model.reaction_names)}.items() if not j}
        # is a dict for a certain gene knockout for reactions that are false (inactive) or None (there's no GPR rule associated) under that gene knockout
        # key is reaction and value is False/None
    # create dict where key is gene to knockout and value is id of reaction(s) that are inactive in universal model when that knockout happens:
    reaction_states_inactive = {k: tuple(frozenset({i for i, j in v.items() if j is not None})) for k, v in reaction_states.items()}
    # create dict where key is a inactive or group of inactive reactions and value is a gene or group of inactive genes associated
    # can be read as, reaction(s) is/are inactive when any of the genes in corresponding set is knockout:
    rko_sets = {}
    for k, v in reaction_states_inactive.items():
        if v not in rko_sets.keys():
            rko_sets[v] = {k}
        else:
            rko_sets[v] |= {k} # |= in dict means update: if inactive reaction is already in rko_sets, joins the new gene to the set of other genes
    # exclude empty reactions, for when a gene knockout doesn't have inactive associated reaction:
    kos_to_test = [l for l in list(rko_sets.keys()) if len(l) > 0]
    # kos_to_test is a list of combination of reactions that are inactive, depending on which gene is knockout
    return kos_to_test, reaction_states_inactive

def batch_ko_optimize(model, kos, context, objective, objective_sense, threads, EssMtb=False):
    '''
    - list where each element represent a comb of reaction kos corresponding to one gene ko, containing simulated lethality scores
      or
      list where each element represent a comb of reaction which bounds change when one metabolite is ko, containing simulated lethality scores
    :param model: generic model adapted for medium composition in cobamp format
    :param kos: list of combination of reactions that are ko when a gene is ko
                or
                list of combination of (reaction, reaction lb, reaction up) which bounds change when a metabolite is knockout
    :param context: dict of (LB,UB) reactions of a reconstructed model
    :param objective: coefficient of model objective
    :param objective_sense: if False, means to maximize the objective
    :param threads: maximum threads available
    :param EssMtb: boolean indicating whether to work with lethal metabolites or not (to work with lethal genes)
    :return: - list where each element represent a comb of reaction kos or of reactions which bounds change.
               * if ko/bounds change is essential (gives infeasible model, or model that doesn't grow) value is 0
               * if ko/bounds change gives positive objective, value is objective value after ko/ bounds change divided by objective value before ko/bound change:
                 in that case if value is below 1 means biomass decreased with ko and ko affects biomass, although not necessaily essential
    '''
    with model as m:
        for crx, Bd in context.items():
            Bd = ast.literal_eval(Bd) # values are as str so eval is needed
            m.set_reaction_bounds(crx, lb=Bd[0], ub=Bd[1]) # change reaction bounds
        sol = m.optimize(objective, objective_sense)
        ofv = sol.objective_value()
        if sol.status() == 'optimal' and ofv > 0:
            if EssMtb:
                # list of dicts where each dict corresponds to set of reactions which bounds have to be change when a metabolite is ko
                # dic is reactionIndex of reaction vs bounds to change to:
                gko_to_rko = [{model.decode_index(rko[0], 'reaction'): (rko[1], rko[2]) for rko in ko_set} for ko_set in kos]  # ko_set is comb of reactions which bounds are to changed when a metabolite is ko
            else:
                # list of dicts where each dict corresponds to set of reactions ko when a gene is ko
                # dic is reactionIndex of reaction to ko as keys and (0,0) as values:
                gko_to_rko = [{model.decode_index(rko, 'reaction'): (0, 0) for rko in ko_set} for ko_set in kos] # ko_set is comb of reactions that are ko when a gene is ko
            # create model object:
            bopt = BatchOptimizer(m.model, threads=int(min(len(gko_to_rko), threads, MP_THREADS)))
            # simulates several models for the same reconstructed model (1 threshold comb) where one gene/metabolite is ko at each time:
            opt_result = bopt.batch_optimize(gko_to_rko, [{model.decode_index(rko, 'reaction'): v for rko, v in objective.items()}] * len(gko_to_rko), [objective_sense] * len(gko_to_rko))
                # objective = {'biomass_human': 1}, '1' is coef of objective
                # objective_sense is False, which means to optimize
                # opt_result - list of 'optimize' solutions - for each of reaction/reaction combinations
            return [k.objective_value() / ofv if (k.status() == 'optimal') and ofv > 0 else 0 for k in opt_result] # if infeasible, model can't get a steady-state distribution, so ko is lethal
        else:
            print('\tReconstructed model objective is 0... skipping')
            return []

def CorrValue(BaseDir, StudyNumber, cobamp_model, kos_to_test, reaction_states_inactive, ach_mset, LethalityCellLine):
    '''
    - get dict where key is threshold combination and value is a dataframe with experimental and simulation scores for each model
    :param BaseDir: basic directory
    :param StudyNumber: study number
    :param cobamp_model: model adpated for medium composition in cobamp format
    :param kos_to_test: list of combination of reactions that are inactive, depending on which gene is knockout
    :param reaction_states_inactive: dict where key is gene to knockout and value is id of reaction(s) that are inactive in universal model when that knockout happens
    :param ach_mset: processed omics measurement set object of experimental lethality scores
    :param LethalityCellLine: id of cell line in lethality experimental data
    :return corr_coef_dict: dict where key is threshold combination and value is a dataframe with experimental and simulation scores for each model
    '''
    dirpath = os.path.join(BaseDir, 'support/LethalityEval/GapFillRes', StudyNumber, 'RecModelBounds')
    # create dict where key is threshold combination and value is dict
    # inner dict has reaction as key and (LB,UB) as value:
    model_dicts= dict()
    for file in os.listdir(dirpath):
        df = pd.read_csv(os.path.join(dirpath, file), sep='\t')
        if not df.empty:
            df = df.T
            df.columns = df.iloc[0]
            df = df.drop(df.index[0]).T
            model_dicts.update(df.to_dict())
    # rebuild reconsructed models and ko combinations of inactive reactions:
    corr_coef_dict = {}
    for mkey, mdict in model_dicts.items():  # for each threshold combination/model
        print(mkey)
        context_dict = {k: v for k, v in mdict.items()} # make a copy of dict with reaction as key and (LB,UB) as value
        NTHREADS = cpu_count() - cpu_count()/4
        bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context=context_dict, objective={'biomass_human': 1}, objective_sense=False, threads=NTHREADS)
        # create dict where key is comb of reactions ko and value is score representing if that combination was essential (0) or affected biomass (<1) or not (>=1):
        opt_res = dict(zip(kos_to_test, bkopt))
        # create dict where key is gene ko and value is score representing if was essential (0) or affected biomass (<1) or not (>=1):
        of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in reaction_states_inactive.items()}
        # if else condition  guarantees that genes
        # which don't have associated inactive reaction when ko get a score of 1 - meaning they don't affect biomass
        # create dataframe where:
        # 1st col is experimental score centered on zero, where negative values are lethal genes
        # 2nd col is score calculated from reconstructed model where essential/lethal gene has score 0, and gene that affected biomass score from 0 to <1, and gene not decreasing biomass as core of 1 or above
        ddf = ach_mset.data.reindex(index=[LethalityCellLine], columns=of_res).append(pd.Series(of_res, name='biomass')).T
        # create dict where key is threshold comb and value is a dataframe with experimental and simulation scores:
        corr_coef_dict[mkey] = ddf
    path = os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumber, 'LethalitySimvsExpScoreDct.pkl')
    if not os.path.exists(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumber)): # create dir if does not exist
        os.makedirs(os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumber))
    saveResult(corr_coef_dict, path)
    return corr_coef_dict

def CorrExpSimulatedEssentialGenes(LethalCellLinePath, modelMedAdap, BaseDir, StudyNumber, LethalityCellLine):
    '''
    - determines correlation between experimental and simulated lethal genes
    :param LethalCellLinePath: path to file with experimental lethality scores for different cell lines
    :param modelMedAdap: model adapted for specific medium composition (exchange reaction of medium metabolites are open and remaining have bounds LB=0, UB=1000)
    :param BaseDir: basic directory
    :param StudyNumber: study number
    :param LethalityCellLine: id of cell line in lethality experimental data
    :return finaldf: dataframe with threshold comb and correlation scores; best threshold combination is printed
    '''
    # preprocess dataframe with experimental lethality scores for genes in dif. cell lines:
    ach_mset, cobamp_model, model_genes = ProcessExpLethalityScoresDf(LethalCellLinePath, modelMedAdap)
    # get inactive reactions when a gene is knockout:
    kos_to_test, reaction_states_inactive = InactReactFromGeneKo(cobamp_model, model_genes)
    # get dataframe with genes in rows, simulation vs experimental tested on columns and values are scores representing lethality, for each model:
    corr_coef_dict = CorrValue(BaseDir, StudyNumber, cobamp_model, kos_to_test, reaction_states_inactive, ach_mset, LethalityCellLine)
    # determine correlation values between essential(lethal) genes in simulation vs experimental data, for dif. lethality thresholds:
    corr_coef_params = {}
    # get dict where key is threshold comb/model and value a dict
    # lower level dict has as key: numbers representing score thresholds to consider gene as lethal; value: correlation coef
    for k, ddf in corr_coef_dict.items(): # for each threshold comb and dataframe of lethality scores
        corr_coef_params[k] = {} # add threshold comb id as key and empty dict as value
        for i,lv in enumerate([-0.5, -0.6, -0.7, -0.8, -0.9]):
            print(i, lv)
            for j,bv in enumerate([0.001, 0.999]):
                print(j, bv)
                corr_coef_params[k][(i,j)] = matthews_corrcoef((ddf['biomass'] < bv), (ddf[LethalityCellLine] < lv))
    # get dataframe of correlation coef where row is threshold comb and column is comb of score cutoffs to consider gene lethal:
    score_dfs = pd.DataFrame.from_dict(corr_coef_params).T.reset_index()
    # split threshold combinations' column into individual columns, each one with one threshold:
    score_dfs[['algorithm', 'Genes', 'thres', 'globalmin', 'globalmax', 'local', 'int_function']] = pd.DataFrame(score_dfs['index'].apply(lambda x: x.split('_') if isinstance(x, str) else [np.nan]*4).to_list())
    score_dfs.drop(columns='index', inplace=True)
    # unpivots the table and renames columns:
    score_df_melt = pd.melt(score_dfs, id_vars=['Genes','thres','globalmin', 'globalmax', 'local']+['algorithm','int_function']).rename(columns={'variable_0': 'essential_cutoff', 'variable_1': 'biomass_cutoff'})
    # find best comb of cutoffs:
    bestCutOffsInd = score_df_melt.groupby(['essential_cutoff', 'biomass_cutoff']).mean().idxmax()[0]
    # apply best comb of cutoffs:
    finaldf = score_df_melt[score_df_melt['essential_cutoff'] == bestCutOffsInd[0]]
    finaldf = finaldf[finaldf['biomass_cutoff'] == bestCutOffsInd[1]]
    fnd = os.path.join(BaseDir, 'support/LethalityEval/SimvsExpScore', StudyNumber)
    if not os.path.exists(fnd):
        os.makedirs(fnd)
    finaldf.to_csv(os.path.join(fnd, 'CorrScores.tab'), sep='\t', header=True)
    besthres = finaldf.loc[finaldf['value'].idxmax()]
    quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    quantiles = list(map(str, quantiles))
    if np.isnan(float(besthres.globalmin)):
        gbmin = str(besthres.globalmin)
    else:
        gbmin = quantiles[int(float(besthres.globalmin))]
    if np.isnan(float(besthres.globalmax)):
        gbmax = str(besthres.globalmax)
    else:
        gbmax = quantiles[int(float(besthres.globalmax))]
    if np.isnan(float(besthres.local)):
        lc = str(besthres.local)
    else:
        lc = quantiles[int(float(besthres.local))]
    print('study:' + StudyNumber, 'strategy:' + besthres.thres, 'globalmin:' + gbmin, 'globalmax:' + gbmax, 'local:' + lc, 'GPR_rules:'+ besthres.int_function, 'algorithm:'+ besthres.algorithm, 'genes:'+besthres.Genes)







