### Modules:
import os
import pandas as pd
import numpy as np
import json5
import copy
from troppo.tasks.task_io import JSONTaskIO
from cobamp.utilities.parallel import batch_run
from cobra import Reaction, flux_analysis
from troppo.tasks.core import Task, TaskEvaluator

### Functions:
def ConvertExcelTasks2Json(TasksExcelPath):
    '''
    - convert tasks in excell format to appropriate troppo tasks object format in Json file
    :param TasksExcelPath: path to tasks' file in .xlsx format
    :return: saves .json file with tasks
    '''
    TasksLst = list()
    ## For each task (row in tasks table)
    metbTasks = pd.read_excel(TasksExcelPath).dropna(subset=['task'])
    metbTasks['task'] = metbTasks['task'].apply(int)
    for k in range(metbTasks.shape[0]):
        print(k)
        TaskInfo = metbTasks.iloc[k, :]  # extract info on the task
        TaskName = str(TaskInfo.loc['task'])
        inflowMet = str(TaskInfo.loc['inputIds']).split(', ') # inflow metabolites
        LBin = str(TaskInfo.loc['LBin']).split(', ')
        UBin = str(TaskInfo.loc['UBin']).split(', ')
        outflowMet = str(TaskInfo.loc['outputIds']).split(', ') # outflow metabolites
        LBout = TaskInfo.loc['LBout'].split(', ')
        UBout = TaskInfo.loc['UBout'].split(', ')
        ReacName = str(TaskInfo.loc['equRctId'])
        ReacSubs = str(TaskInfo.loc['equationSubs']).split(', ')
        ReacProd = str(TaskInfo.loc['equationProd']).split(', ')
        ReacCoefSubs = str(TaskInfo.loc['equationCoeficientSubs']).split(', ')
        ReacCoefProd = str(TaskInfo.loc['equationCoeficientProd']).split(', ')
        LBreac = TaskInfo.loc['LBequ']
        UBreac = TaskInfo.loc['UBequ']
        inflowDic = dict()
        outflowDic = dict()
        for i in range(len(inflowMet)):
            inbounds = list()
            metName = inflowMet[i]
            inbounds.append(float(LBin[i]))
            inbounds.append(float(UBin[i]))
            inflowDic[metName] = inbounds
        for j in range(len(outflowMet)):
            outbounds = list()
            metName = outflowMet[j]
            outbounds.append(float(LBout[j]))
            outbounds.append(float(UBout[j]))
            outflowDic[metName] = outbounds
        if ReacSubs != ['nan']:
            metSubsCoefDic = {ReacSubs[l]: (-1)*int(ReacCoefSubs[l]) for l in range(len(ReacSubs))}
            metProdCoefDic = {ReacProd[l]: int(ReacCoefProd[l]) for l in range(len(ReacProd))}
            metCoefDic = metSubsCoefDic.copy()
            metCoefDic.update(metProdCoefDic)
            reactionDic = {ReacName: tuple((metCoefDic, (float(LBreac), float(UBreac))))}
            task = Task(
            should_fail = False,
            inflow_dict = inflowDic,
            outflow_dict = outflowDic,
            reaction_dict = reactionDic,
            name = TaskName,
            flux_constraints = {},
            mandatory_activity = [],
            annotations = {'subsystem': str(TaskInfo.loc['subsystem']), 'system': str(TaskInfo.loc['system']),
            'id': str(TaskInfo.loc['task']), 'description': str(TaskInfo.loc['description'])})
        else:
            task = Task(
            should_fail = False,
            inflow_dict = inflowDic,
            outflow_dict = outflowDic,
            reaction_dict= {},
            name = TaskName,
            flux_constraints= {},
            mandatory_activity = [],
            annotations = {'subsystem': str(TaskInfo.loc['subsystem']), 'system': str(TaskInfo.loc['system']),
            'id': str(TaskInfo.loc['task']), 'description': str(TaskInfo.loc['description'])})
        TasksLst.append(task)
    return TasksLst

def preprocessTasks(tasks, model):
    '''
    - process tasks: exclude tasks with metabolites not present in generic model adapted for medium composition;
                     and deal with task metabolites that enter and go out of system at same time
    :param tasks: list with tasks
    :param model: generic model adapted for medium composition
    :return task_list: processed tasks
    '''
    # if a task's metabolites are not in the generic model with medium adapted, exclude task:
    task_list = [t for t in tasks if len((set(t.inflow_dict) | set(t.outflow_dict)) - set([m.id for m in model.metabolites])) == 0]
    # if a task's metabolite is both entering and going out of model, then bounds are set to (-1000, 1000) and it is just kept in the list of 'inflow' metabolites:
    for task in task_list:
        task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in task.inflow_dict.items()}  # k is metabolite id, v is bounds
        task.outflow_dict = {k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}
        #task.mandatory_activity = [] # and mandatory_activity (list of reactions that should have flux when doing that task) is set to [].cause we are not interested in using this
    return task_list

def RemoveFailedTasks(tasksName, batch_res_tasks, task_list, tasksDir):
    '''
    - remove tasks failed in generic model from list of tasks to test in reconstructed models
    :param tasksName: tasks' names
    :param batch_res_tasks: info on whether task fails (False) or not
    :param task_list: list of tasks
    :param tasksDir: path to directory where to save outup files
    :return task_listNew: list of tasks without those that failed on generic model although should not fail (should_fail=False)
    '''
    # get tasksName and batch_res_tasks tasks ordered as task_list (from task 1 to end):
    ind = np.argsort(list(map(int, tasksName)))  # get indexes to get tasks in order
    tasksName = [tasksName[i] for i in ind]  # get tasks names in order
    batch_res_tasks = [batch_res_tasks[i] for i in ind]  # as batch_res_tasks elements have same order as tasksName (initially), we can order batch_res_tasks in this way
    tasksFailed = list()
    for taskN, task, t in zip(tasksName, batch_res_tasks, task_list): # checks which tasks fail
        if (not task[0]) and (not t.should_fail): # tasks that fail in universal model but that should_fail=True (they should fail), may not fail in reconstructed models, so should be kept to bet tested in reconstructed models
            print('Task ' + taskN + ' failed on generic model. It will not be tested on reconstructed models')
            tasksFailed.append(taskN)
    task_listNew = [el for el in task_list if el.name not in tasksFailed]
    pd.Series(tasksFailed).to_csv(tasksDir + '/GenericModelFailedTasks.tab', sep='\t', index=False, header=False)  # saves failed tasks into a file
    return task_listNew

def EvalAllTasksAtSameTime(taskModel, task_list):
    '''
    - evaluate all tasks in generic model
    :param taskModel: model in which tasks are tested
    :param task_list: list of tasks to test
    :return tasksName: tasks' names
    return batch_res_tasks: info on whether each task failed (False) or not (True)
    '''
    tasksEval = TaskEvaluator(model=taskModel, tasks=task_list, solver='CPLEX')
    tasksName = tasksEval.tasks
    batch_res_tasks = batch_run(TaskEvaluator.batch_function, tasksName, {'tev': tasksEval}, threads=10) # does evaluation of all tasks at same time
    return tasksName, batch_res_tasks

def TestTasksUniversalModel(TasksExcelPath, TasksJsonPath, model):
    '''
    - convert tasks from excell to json, when json file format is not available
    - evaluate tasks on generic model
    - save info on tasks that were failed by generic model and remove those from the list of tasks to further test in reconstructed models
    :param TasksExcelPath: path to tasks' file in .xlsx format
    :param TasksJsonPath: path to .json file where to save/get tasks
    :param model: generic model adapted for medium composition
    :return task_list: group of tasks that generic model was able to perform
    '''
    if not os.path.exists(TasksJsonPath):
        ## convert tasks in excell format to appropriate format in Json file
        tasks = ConvertExcelTasks2Json(TasksExcelPath)
        ## preprocess tasks:
        taskModel = model.copy()
        # close opened boundary reactions of the model to evaluate tasks:
        for k in taskModel.boundary:
            k.knock_out()
        task_list = preprocessTasks(tasks, taskModel)
        # evaluate tasks on generic model:
        tasksName, batch_res_tasks = EvalAllTasksAtSameTime(taskModel, task_list)
        # failed tasks on universal model will be removed from list of tasks to test in reconstructed model:
        tasksDir = '/'.join(TasksJsonPath.split('/')[:len(TasksJsonPath.split('/'))-1])
        task_list = RemoveFailedTasks(tasksName, batch_res_tasks, task_list, tasksDir)
        JSONTaskIO().write_task(TasksJsonPath, task_list)

def AddTaskCobrapy(tasksDic, UnivModel, IncompModel, Tname):
    '''
    - adds a task to generic metabolic model and reconstructed model that have had their boundary reactions previously closed
    :param tasksDic: dictionary with name of task and corresponding task. each task is in troppo tasks object format
    :param UnivModel: generic model which has their boundary reactions closed
    :param IncompModel: reconstructed model which has their boundary reactions closed
    :param Tname: name of task to add
    :return UnivModel: above mentioned UnivModel but with task added
    :return IncompModel: above mentioned IncompModel but with task added
    :return taskReacIds: ids of reactions of the task that was added
    '''
    taskInfo = tasksDic[Tname]
    taskReacIds = list() # list will have all ids of reactions added in one task
    ## Add inflow reactions:
    for m, v in taskInfo.inflow_dict.items():
        rctId = m + '_IN'
        taskReacIds.append(rctId)
        met = UnivModel.metabolites.get_by_id(m).copy()
        rc = Reaction(rctId) # create a new reaction with reaction id as 'rctId'
        LB = float(v[1]) * (-1) # bounds of 'in' metabolites in taskDict are as for eg. (2,1000), as these are 'in' reactions in model corresponds to (-1000, -2)
        UB = float(v[0]) * (-1)
        rc.add_metabolites({met: -1.0})  # metabolite on left side of equation
        rc.bounds = (LB, UB)
        IncompModel.add_reactions([rc])
        rc2 = rc.copy() # needed cause if we apply same reaction to both incomplete and universal models, when we try to remove the reactions in both models gives error
        UnivModel.add_reactions([rc2])
    for m, v in taskInfo.outflow_dict.items():
        rctId = m + '_OUT'
        taskReacIds.append(rctId)
        met = UnivModel.metabolites.get_by_id(m).copy()
        rc = Reaction(rctId)
        LB = float(v[0])
        UB = float(v[1])
        rc.add_metabolites({met: -1.0})  # metabolite on left side of equation
        rc.bounds = (LB, UB)
        IncompModel.add_reactions([rc])
        rc2 = rc.copy()
        UnivModel.add_reactions([rc2])
    if bool(taskInfo.reaction_dict):  # if task has a additional reaction
        for k, v in taskInfo.reaction_dict.items():
            rctId = k
            metdic = v[0]
            taskReacIds.append(rctId)
            rc = Reaction(rctId)
            rc.add_metabolites({UnivModel.metabolites.get_by_id(m).copy(): ind for m, ind in metdic.items()})
            rc.bounds = tuple(v[1])
            IncompModel.add_reactions([rc])
            rc2 = rc.copy()
            UnivModel.add_reactions([rc2])
    return UnivModel, IncompModel, taskReacIds

def TaskGapfill(UnivModel, IncompModel):
    '''
    - gapfill for tasks using cobrapy gapfill
    :param UnivModel: generic model with task to gapfill added and boundary reactions closed
    :param IncompModel: reconstructed model with task to gapfill added and boundary reactions closed
    :return gapfillTasksolIds: ids of reactions that are needed to gapfill the task
    '''
    gapfiller = flux_analysis.gapfilling.GapFiller(IncompModel, universal=UnivModel, exchange_reactions=False, demand_reactions=False, integer_threshold=1e-9, lower_bound=0)  # lower_bound = 0: for a task to be passed model has to be feasible, but no need for any objective reaction flux to be positive. there are no objectives in tasks
    gapfiller.model.solver.problem.parameters.timelimit.set(1200) # gapfill can take a max of 20 min (1200 sec) otherwise it will continue without solution
    gapfiller.model.solver.configuration.tolerances.feasibility = 1e-9
    gapfiller.model.solver.configuration.tolerances.integrality = 1e-9
    gapfiller.model.solver.configuration.tolerances.optimality = 1e-9
    gapfillTasksol = gapfiller.fill()[0]
    gapfillTasksolIds = [reac.id for reac in gapfillTasksol]
    return gapfillTasksolIds

def EssTaskRcFct(path, BaseDir):
    '''
    - saves all reactions in generic model that if excluded cause one or more of the essential tasks to fail
    Note: prior to run this function, 'prepare4cluster.py' script has to be run and 'submitJobsEssRcKo.sh' has to be run on a cluster
    :param BaseDir: basic directory
    :param path: path to file with results
    '''
    TasksReacKoLst = list()
    fileLst = os.listdir(os.path.join(BaseDir, 'support/LethalityEval/EssRcKo'))
    for file in fileLst:
        rs = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/EssRcKo', file)).iloc[:,0].tolist()
        TasksReacKoLst = TasksReacKoLst + rs
    pd.Series(TasksReacKoLst).to_csv(path, sep='\t', index=False, header= False)

### Functions:
def EvalTasksReacKo(rid, md, task_list):
    '''
    - knock out (ko) to a given reaction and if any task fails returns True, otherwise returns False
    :param rid: id of reaction to knock out
    :param md: model provided
    :param task_list: list of tasks to test
    '''
    r = md.reactions.get_by_id(rid)
    bd = r.bounds
    r.knock_out()
    tasksEval = TaskEvaluator(model=md, tasks=task_list, solver='CPLEX')  # create a task evaluator instance
    for tsk in task_list:
        tn = tsk.name # get task name
        tasksEval.current_task = tn # set task to evaluate
        taskres = tasksEval.evaluate()[0] # evaluate task
        tasksEval.current_task = None  # remove current task from evaluation - back to normal
        if not taskres: # if task fails, append rcid to list and move to next reaction
            r.bounds = bd
            return rid, True
    else:
        r.bounds = bd
        return rid, False

def ConvertTasksRecon2Human1(TJsonPath, MetAssPath, finalpath, model):
    '''
    - convert tasks with metabolite ids of recon3d to tasks with metabolite ids of human1,
      also excludes tasks not working on generic model from list of tasks to test
    :param TJsonPath: path to.json file with tasks from recon3d (consensus) list
    :param MetAssPath: path to.json file with association between recon3d and human1 metabolites' ids
    :param finalpath: path to json file that contains the result, a list of tasks that can be tested on human1 model
    :param model: generic metabolic model
    '''
    ## get tasks with metabolite ids in recon3d format:
    tasks = JSONTaskIO().read_task(TJsonPath)
    ## get dataframe with association between metabolites ids in recon3d and ihuman and compartments:
    with open(MetAssPath) as json_file:
        metbass = json5.load(json_file)
    df = pd.concat([pd.Series(metbass['mets']), pd.Series(metbass['metRecon3DID'])], axis=1)
    df.columns = ['Human1', 'Recon3D']
    df['Compartment'] = df['Human1'].apply(lambda x: x[-1])
    ## exclude tasks with metabolites that do not exist in human1:
    t2Remove = list()
    for t in tasks:
        INB = [True for im, v in t.inflow_dict.items() if im[:-3] not in list(df['Recon3D'])]  # get True if not all IN metabolites can be converted to iHuman ID
        OUTB = [True for im, v in t.outflow_dict.items() if im[:-3] not in list(df['Recon3D'])]  # get True if not all OUT metabolites can be converted to iHuman ID
        if len(INB) != 0 or len(OUTB) != 0:
            t2Remove.append(t.name)
    tasks = [t for t in tasks if t.name not in t2Remove]  # remove tasks
    ## replace ids of metabolites from recon3d to human1:
    # note: all tasks to test have reaction_dict empty: [t.name for t in tasks if len(t.reaction_dict) != 0] is []
    Mt2Remove = list()
    tasksInD = dict()
    tasksOutD = dict()
    for t in tasks:
        inDct = {''.join([im[:-2], 's]']) if (im[-2] == 'x' or im[-2] == 'e') else im: v for im, v in t.inflow_dict.items()}  # replace x and e compartments by s (both correspond to extracellular comp in human1)
        outDct = {''.join([im[:-2], 's]']) if (im[-2] == 'x' or im[-2] == 'e') else im: v for im, v in t.outflow_dict.items()}  # replace x and e compartments by s (both correspond to extracellular comp in human1)
        inD = dict()
        outD = dict()
        for im, v in inDct.items():  # for inflow metabolites
            if len(df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]) == 0:  # if combination of metabolite id + compartment do not exist in human1
                Mt2Remove.append(t.name)  # append task, to latter on remove it
            else:  # else, replace recon3d metabolite id by id in human1
                inD[df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]['Human1'].iloc[0]] = v
        for im, v in outDct.items():  # for outflow metabolites
            if len(df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]) == 0:  # if combination of metabolite id + compartment do not exist in human1
                Mt2Remove.append(t.name)  # append task, to latter on remove it
            else:  # else, replace recon3d metabolite id by id in human1
                outD[df[(df['Recon3D'] == im[:-3]) & (df['Compartment'] == im[-2])]['Human1'].iloc[0]] = v
        tasksInD[t.name] = inD
        tasksOutD[t.name] = outD
    tasks = [t for t in tasks if t.name not in set(Mt2Remove)]  # remove reactions with combination of metabolite id + compartment that do not exist in human1
    tasksN = copy.deepcopy(tasks)
    for t in tasksN:
        t.inflow_dict = tasksInD[t.name]
        t.outflow_dict = tasksOutD[t.name]
        t.mandatory_activity = [] # mandatory activity is specific reactions that are required to have flux. But it is not necessary to test the task
    ## preprocess tasks:
    taskModel = model.copy()
    # close opened boundary reactions of the model to evaluate tasks:
    for k in taskModel.boundary:
        k.knock_out()
    task_list = preprocessTasks(tasksN, taskModel)
    # evaluate tasks on generic model:
    tasksName, batch_res_tasks = EvalAllTasksAtSameTime(taskModel, task_list)
    # failed tasks on universal model will be removed from list of tasks to test in reconstructed model:
    tasksDir = '/'.join(finalpath.split('/')[:len(finalpath.split('/')) -1])
    task_list = RemoveFailedTasks(tasksName, batch_res_tasks, task_list, tasksDir)
    JSONTaskIO().write_task(os.path.join(tasksDir, 'processedTasks.json'), task_list)







