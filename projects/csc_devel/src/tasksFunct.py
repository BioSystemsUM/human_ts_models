### Modules:
import os
import pandas as pd
import numpy as np
from troppo.tasks.task_io import JSONTaskIO
from troppo.tasks.core import Task, TaskEvaluator
from cobamp.utilities.parallel import batch_run
from multiprocessing.dummy import Pool
from multiprocessing import cpu_count
from functools import partial
from scaleGenExpFunct import saveResult
from cobra import Reaction, flux_analysis

### Functions:
def ConvertExcelTasks2Json(TasksExcelPath, TasksJsonPath):
    '''
    - convert tasks in excell format to appropriate troppo tasks object format in Json file
    :param TasksExcelPath: path to tasks' file in .xlsx format
    :param TasksJsonPath: path to .json file where to save/get tasks
    :return: saves .json file with tasks
    '''
    TasksLst = list()
    ## For each task (row in tasks table)
    metbTasks = pd.read_excel(TasksExcelPath)
    for k in range(metbTasks.shape[0]):
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
    JSONTaskIO().write_task(TasksJsonPath, TasksLst)

def preprocessTasks(TasksJsonPath, model):
    '''
    - process tasks: exclude tasks with metabolites not present in generic model adapted for medium composition;
                     and deal with task metabolites that enter and go out of system at same time
    :param TasksJsonPath: path to .json file with tasks
    :param model: generic model adapted for medium composition
    :return task_list: processed tasks
    '''
    tasks = JSONTaskIO().read_task(TasksJsonPath)
    # if a task's metabolites are not in the generic model with medium adapted, exclude task:
    task_list = [t for t in tasks if len((set(t.inflow_dict) | set(t.outflow_dict)) - set([m.id for m in model.metabolites])) == 0]
    # if a task's metabolite is both entering and going out of model, then bounds are set to (-1000, 1000) and it is just kept in the list of 'inflow' metabolites:
    for task in task_list:
        task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in task.inflow_dict.items()}  # k is metabolite id, v is bounds
        task.outflow_dict = {k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}
        #task.mandatory_activity = [] # and mandatory_activity (list of reactions that should have flux when doing that task) is set to [].cause we are not interested in using this
    return task_list

def RemoveFailedTasks(tasksName, batch_res_tasks, task_list, TasksJsonPath):
    '''
    - remove tasks failed in generic model from list of tasks to test in reconstructed models
    :param tasksName: tasks' names
    :param batch_res_tasks: info on whether task fails (False) or not
    :param task_list: list of tasks
    :param TasksJsonPath: path to .json file with tasks (used to retrieve tasks directory)
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
    tasksDir = '/'.join(TasksJsonPath.split('/')[:len(TasksJsonPath.split('/'))-1])
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
        ConvertExcelTasks2Json(TasksExcelPath, TasksJsonPath)
    ## preprocess tasks:
    taskModel = model.copy()
    # close opened boundary reactions of the model to evaluate tasks:
    for k in taskModel.boundary:
        k.knock_out()
    task_list = preprocessTasks(TasksJsonPath, taskModel)
    # evaluate tasks on generic model:
    tasksName, batch_res_tasks = EvalAllTasksAtSameTime(taskModel, task_list)
    # failed tasks on universal model will be removed from list of tasks to test in reconstructed model:
    task_list = RemoveFailedTasks(tasksName, batch_res_tasks, task_list, TasksJsonPath)
    return task_list

'''
def EssTaskRcFct(model, task_list, path):
    
    - saves a .pkl object with list of reactions that are essential for tasks performed in a given model
    :param model: where tasks are performed
    :param task_list: list of tasks to perform
    :param path: path to pickle file where to save results
   
    taskModel = modelMedAdap.copy()
    EssTaskRcLst = list()
    for r in taskModel.boundary:
        r.knock_out()
    r2nk = set(taskModel.reactions) - set(taskModel.boundary)
    for r in r2nk:
        print(r.id)
        rb = r.bounds
        print(taskModel.reactions.get_by_id(r.id).bounds)
        r.knock_out()
        print(taskModel.reactions.get_by_id(r.id).bounds)
        tasksName, batch_res_tasks = EvalAllTasksAtSameTime(taskModel, task_list)
        for taskN, task in zip(tasksName, batch_res_tasks):  # checks which tasks fail
            if not task[0]:
                EssTaskRcLst.append(r.id)
                break
        r.bounds = rb
        print(taskModel.reactions.get_by_id(r.id).bounds)
'''
def EssTaskRcFct(model, task_list, path, BaseDir):
    '''
    - Identify/saves reactions in generic model that if excluded cause one or more of the essential tasks to fail
    :param model: generic model may be adapted or not for medium composition)
    :param task_list: list of tasks to test
    :param path: path to file with results
    '''
    # function that given a rid does knockout (ko) to corresponding reaction and checks whether that ko is essential for at least one task:
    def EvalTasksReacKo(rid, model, task_list):
        print(rid)
        md = model.copy() # to be able to to ko to one reaction at a time while using multiprocessing
        md.reactions.get_by_id(rid).knock_out()
        tasksName, batch_res_tasks = EvalAllTasksAtSameTime(taskModel, task_list)
        for taskN, task in zip(tasksName, batch_res_tasks):  # checks which tasks fail
            if not task[0]:
                print('end', rid)
                return rid, True
        else:
            print('end', rid)
            return rid, False
    # close model exchange/boundary reactions (in "human1" all boundaries are exchange reactions):
    taskModel = model.copy()
    for r in taskModel.boundary:
        r.knock_out()
    # get all reactions of generic model except boundaries/exchange reactions:
    rids = {r.id for r in taskModel.reactions}
    bdrids = {r.id for r in taskModel.boundary}
    frids = list(rids - bdrids)
    grps = 30
    integ = len(frids)//grps
    remain = len(frids)%grps
    intv = grps*[integ] + [remain]
    i = 0
    for v in intv:
        print(i)
        f = i + v
        pd.Series(frids[i:f]).to_csv(os.path.join(BaseDir, 'support/LethalityEval/EssLst' + str(i)), index=False)
        i = i + v
    TasksReacKoLst = list()
    fileLst = os.listdir(os.path.join(BaseDir, 'support/LethalityEval/EssRcKo'))
    for file in fileLst:
        rs = pd.read_csv(os.path.join(BaseDir, 'support/LethalityEval/EssRcKo', file)).iloc[:,0].tolist()
        TasksReacKoLst = TasksReacKoLst + rs
    len(TasksReacKoLst)
    pd.Series(TasksReacKoLst).to_csv(path, sep='\t', index=False, header= False)
    # use multiprocessing to make faster the evaluation of each reaction ko:
    pool = Pool(processes=int(cpu_count() - (cpu_count()/4)))
    fc = partial(EvalTasksReacKo, model=taskModel, task_list=task_list)
    results = pool.map(fc, frids)
    pool.close()
    pool.join()
    # save reactions essential to at least one task into a list:
    TasksReacKoLst = list()
    for rid, bool in results:
        if bool:
            TasksReacKoLst.append(rid)
    pd.Series(TasksReacKoLst).to_csv(path, sep='\t', index=False)
'''

11:42 - 59
23seg * 10400=239200
239200 / 60 /60 =66 horas... 3 dias...
    def EvalTasksReacKo(rid, model, task_list):
        print(rid)
        md = model.copy() # to be able to to ko to one reaction at a time while using multiprocessing
        md.reactions.get_by_id(rid).knock_out()
        tasksEval = TaskEvaluator(model=md, tasks=task_list, solver='CPLEX')  # create a task evaluator instance
        for tsk in task_list:
            tn = tsk.name # get task name
            tasksEval.current_task = tn # set task to evaluate
            taskres = tasksEval.evaluate()[0] # evaluate task
            tasksEval.current_task = None  # remove current task from evaluation - back to normal
            if not taskres: # if task fails, append rcid to list and move to next reaction
                return rid, True
        else:
            return rid, False
        '''

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








