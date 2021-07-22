import sys, os

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])


from cobamp.wrappers.external_wrappers import get_model_reader


import os
import pandas as pd
import numpy as np

from troppo.methods_wrappers import GapfillWrapper

from cobamp.utilities.file_io import pickle_object, read_pickle
from cobamp.utilities.parallel import MP_THREADS
from troppo.tasks.core import TaskEvaluator
from projects.breast_mcf7.scripts.human1.models import get_human1_model, get_human1_essential_tasks

AVAILABLE_THREADS = 60

ROOT_FOLDER = 'projects/breast_mcf7'
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/mcf7_comparison')
MODEL_PATH = os.path.join(ROOT_FOLDER, 'support/Human-GEM.xml')
TASK_RESULTS_FOLDER = os.path.join(CS_MODEL_DF_FOLDER, 'task_evaluation')
if not os.path.exists(TASK_RESULTS_FOLDER): os.makedirs(TASK_RESULTS_FOLDER)

CS_MODEL_SET_NAME = 'cs_models_all_combinations_biomass'

TASK_RESULTS_PATH = os.path.join(TASK_RESULTS_FOLDER, CS_MODEL_SET_NAME+'_taskeval.pkl')
TASK_RESULTS_POST = os.path.join(TASK_RESULTS_FOLDER, CS_MODEL_SET_NAME+'_taskeval_postgapfill.pkl')
TASK_GAPFILL_PATH = os.path.join(TASK_RESULTS_FOLDER, CS_MODEL_SET_NAME+'_taskgapfill.pkl')
MODEL_DF_PATH = os.path.join(CS_MODEL_DF_FOLDER, CS_MODEL_SET_NAME+'.csv')
MODEL_DF_PATH_AFTER_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER,CS_MODEL_SET_NAME+'_postgapfill.csv')

ESS_RESULTS_FOLDER = os.path.join(CS_MODEL_DF_FOLDER, 'essentiality')
if not os.path.exists(ESS_RESULTS_FOLDER): os.makedirs(ESS_RESULTS_FOLDER)

TASK_ESS_FOLDER = os.path.join(ESS_RESULTS_FOLDER, CS_MODEL_SET_NAME+'_task_essentials/')
if not os.path.exists(TASK_ESS_FOLDER): os.makedirs(TASK_ESS_FOLDER)

model = get_human1_model()
task_list = get_human1_essential_tasks(model)
model_df = pd.read_csv(MODEL_DF_PATH, index_col=[0,1,2])

task_model = model.copy()
for rx in task_model.boundary: rx.bounds = (0, 0)
task_model_no_boundary = task_model.copy()
task_model_no_boundary.remove_reactions(task_model_no_boundary.boundary)
task_model_no_boundary.remove_reactions([k for k in task_model_no_boundary.reactions if len(k.metabolites) == 0])


## Evaluate tasks
if not os.path.exists(TASK_RESULTS_PATH+'.'+'0'):
    tev = TaskEvaluator(model=task_model_no_boundary, tasks=task_list)

    task_eval_results = {}
    for i, sub_model_df in enumerate(np.array_split(model_df[[r.id for r in task_model_no_boundary.reactions]], AVAILABLE_THREADS)):
        print('Model chunk #'+str(i+1))
        bound_changes = [{k:(0, 0) for k,v in row.to_dict().items() if not v}
                         for i, row in sub_model_df.iterrows()]
        task_eval_results.update({sub_model_df.index[k]: {i:j[0] for i,j in v.items()}
                                  for k,v in tev.batch_evaluate(bound_changes, AVAILABLE_THREADS).items()})
        pickle_object(task_eval_results,TASK_RESULTS_PATH+'.'+str(i))
else:
    from cobamp.utilities.file_io import read_pickle
    task_eval_results = {}
    for k in range(12):
        cur_len = len(task_eval_results)
        task_eval_results.update(read_pickle(TASK_RESULTS_PATH+'.'+str(k)))

tevr = pd.DataFrame(task_eval_results).T
tevr.sum(axis=1)

tev = TaskEvaluator(model=task_model_no_boundary, tasks=task_list)
cob_cbm = get_model_reader(task_model_no_boundary, ttg_ratio=99999).to_cobamp_cbm('CPLEX')
all_genes = cob_cbm.gpr.get_genes()
state_dicts = (dict([(k,True) for k in all_genes]+[(g, False)]) for g in all_genes)
evals = {i:[cob_cbm.gpr.eval_gpr(j, state) for j in range(len(cob_cbm.gpr))] for i, state in enumerate(state_dicts)}
kos = {all_genes[k]:frozenset(cob_cbm.reaction_names[k] for k,v in enumerate(d) if v == False) for k,d in evals.items()}
unique_kos = list(set(kos.values()))

rx_set_to_gene_map = {k:[] for k in kos.values()}
for gene_ko, rx_ko_set in kos.items():
    rx_set_to_gene_map[rx_ko_set].append(gene_ko)

model_dicts_eval = {k:{x for x,y in v.items() if not y}
                    for k,v in model_df[[r.id for r in task_model_no_boundary.reactions]].T.to_dict().items()}

model_dict_keys = list(model_dicts_eval.keys())
model_dicts_eval_list = [model_dicts_eval[k] for k in model_dict_keys]
unique_ko_dict_list = [{i:(0,0) for i in k} for k in unique_kos]

model_dfs = {}
for model_name, row in zip(model_dict_keys, model_dicts_eval_list):
    print(model_name)
    with task_model_no_boundary as m:
        for rxid in row:
            m.reactions.get_by_id(rxid).bounds = (0,0)
        tev = TaskEvaluator(model=m, tasks=task_list)
    dfs = []
    ko_task_evaluation = tev.batch_evaluate(
        bound_changes=unique_ko_dict_list,threads=AVAILABLE_THREADS, output_sol=False, mp_batch_size=50000)
    dfs.append(pd.DataFrame([{k:v[0] for k,v in ko_task_evaluation[i].items()} for i in range(len(unique_ko_dict_list))]))
    model_dfs[model_name] = pd.concat(dfs)
    model_dfs[model_name].index = unique_kos
    model_dfs[model_name].to_csv(os.path.join(TASK_ESS_FOLDER, '-'.join(model_name)+'.csv'))

# def task_flows_to_overrides(task):
#     ls_override = {k:[] for k in ['non_consumed', 'consumed', 'non_produced', 'produced']}
#
#     for k,v in task.inflow_dict.items():
#         if v[0] == 0: ls_override['non_produced'].append(k)
#         else: ls_override['consumed'].append(k)
#
#     for k, v in task.outflow_dict.items():
#         if v[0] == 0: ls_override['non_consumed'].append(k)
#         else: ls_override['produced'].append(k)
#
#     return ls_override
#
#
# def gapfill_local(missing, params):
#     print('Gapfilling from', len(missing))
#     steady_state_override = params['sso']
#     cfg = params['cfg']
#     gw = params['gw']
#     return gw.run(avbl_fluxes=missing, algorithm='efm',
#                   ls_override=steady_state_override, kshproperties=cfg)
#
# # Gapfill for tasks
# if not os.path.exists(TASK_GAPFILL_PATH):
#     from troppo.methods.gapfill.efm import DEFAULT_CONFIG
#
#     DEFAULT_CONFIG['BIGMVALUE'] = 1e6
#     DEFAULT_CONFIG['N_THREADS'] = 1
#     gw_generic = GapfillWrapper(task_model_no_boundary)
#     task_gapfill_results = {}
#     for task in task_list:
#         print('Gapfilling for',task)
#         task_gapfill_results[task.name] = {}
#         models_to_gapfill = tevr[tevr[task.name] == False].index
#         if len(models_to_gapfill) > 0:
#             task_gapfill_model_df = model_df[set(model_df.columns) & set([r.id for r in task_model_no_boundary.reactions])]
#             if len(task.reaction_dict) > 0:
#                 tev = TaskEvaluator(model=task_model_no_boundary, tasks=[task])
#                 tev.model.initialize_optimizer()
#                 tev.model.remove_reactions({k for k in task.get_task_bounds() if 'flow' in k.split('_')[-1]})
#                 gw = GapfillWrapper(tev.model)
#             else:
#                 gw = gw_generic
#
#             missing_reactions = [{k for k,v in task_gapfill_model_df.loc[cs_model_index,:].to_dict().items() if not v}
#                                  for cs_model_index in models_to_gapfill]
#
#             from cobamp.utilities.parallel import batch_run
#
#             mp_gapfill_result = batch_run(gapfill_local, missing_reactions,
#                                           paramargs={'sso':task_flows_to_overrides(task),
#                                                      'gw': gw,
#                                                      'cfg': DEFAULT_CONFIG},
#                                           threads=min(len(missing_reactions), AVAILABLE_THREADS//3))
#
#             task_gapfill_results[task.name] = dict(zip(models_to_gapfill, mp_gapfill_result))
#     pickle_object(task_gapfill_results, TASK_GAPFILL_PATH)
#
# else:
#     task_gapfill_results = read_pickle(TASK_GAPFILL_PATH)
#
#
#
# gapfill_reactions_df = pd.DataFrame(task_gapfill_results).T.apply(
#     lambda y: [x[0] if isinstance(x, list) and len(x) > 0  else [] for x in y])
#
# from itertools import chain
#
#
# model_df_gapfill = model_df.copy()
# change_dict = {k: set(chain(*gapfill_reactions_df[k].to_list())) for k in gapfill_reactions_df}
# for k,v in change_dict.items():
#     model_df_gapfill.loc[k,v] = True
#
# model_df_gapfill.to_csv(MODEL_DF_PATH_AFTER_GAPFILL)

# Evaluate tasks again
# if not os.path.exists(TASK_RESULTS_POST+'.'+'0'):
#     tev = TaskEvaluator(model=task_model_no_boundary, tasks=task_list)
#
#     task_eval_results = {}
#     for i, sub_model_df in enumerate(np.array_split(model_df_gapfill[[r.id for r in task_model_no_boundary.reactions]], 6)):
#         print('Model chunk #'+str(i+1))
#         bound_changes = [{k:(0, 0) for k,v in row.to_dict().items() if not v}
#                          for i, row in sub_model_df.iterrows()]
#         task_eval_results.update({sub_model_df.index[k]: {i:j[0] for i,j in v.items()}
#                                   for k,v in tev.batch_evaluate(bound_changes, AVAILABLE_THREADS).items()})
#         pickle_object(task_eval_results,TASK_RESULTS_POST+'.'+str(i))
# else:
#     from cobamp.utilities.file_io import read_pickle
#     task_eval_results = {}
#     for k in range(6):
#         cur_len = len(task_eval_results)
#         task_eval_results.update(read_pickle(TASK_RESULTS_POST+'.'+str(k)))
#

# pd.DataFrame(task_eval_results).sum(axis=1)