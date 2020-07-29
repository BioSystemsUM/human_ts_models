import sys

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

import os
import re
import pandas as pd
import numpy as np

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model
from sklearn.metrics import matthews_corrcoef

from cobamp.core.optimization import BatchOptimizer
from cobamp.utilities.file_io import pickle_object, read_pickle
from cobamp.utilities.parallel import MP_THREADS
from cobamp.wrappers.external_wrappers import get_model_reader
from troppo.omics.core import IdentifierMapping, TypedOmicsMeasurementSet
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO

ROOT_FOLDER = 'projects/breast_mcf7'
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
TASKS_PATH = 'shared/task_sets/nl2019_tasks_h1_compact.json'
TASK_RESULTS_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/task_evaluation/')
TASK_RESULTS_PATH = os.path.join(TASK_RESULTS_FOLDER, 'cs_models_lt2_taskeval.json')
ACHILLES_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/Achilles_gene_effect.csv'
BOF_ESS_PATH = 'projects/breast_mcf7/results/human1/reconstructions/cs_models_all_combinations_essentiality.pkl'
BOF_PAR_ESS_PATH = 'projects/breast_mcf7/results/human1/reconstructions/cs_models_all_combinations_essentiality_par.pkl'
ESS_SCORE_PATH = 'projects/breast_mcf7/results/human1/reconstructions/cs_models_all_combinations_essentiality_par.csv'
TASK_ESS_PATH = 'projects/breast_mcf7/results/human1/reconstructions/cs_models_all_combinations_task_essentiality.pkl'

if not os.path.exists(MODEL_PATH):
    path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/modelFiles/xml/HumanGEM.xml')
    model = read_sbml_model(path)
    write_sbml_model(model, MODEL_PATH)
else:
    model = read_sbml_model(MODEL_PATH)
model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
blocked = find_blocked_reactions(model)
model.remove_reactions(blocked)

cobamp_model = get_model_reader(model, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
model_genes = cobamp_model.gpr.get_genes()

model_df = pd.read_csv(os.path.join(CS_MODEL_DF_FOLDER, 'cs_models_all_combinations_efm_gapfill_open_media.csv'), index_col=[0, 1, 2])

from cobamp.utilities.file_io import read_pickle
baseline = read_pickle('projects/breast_mcf7/data/human1_gems/a.pkl')

model_df = model_df.append(pd.Series(model_df.columns.isin(baseline), index=model_df.columns,
                        name=('tinit_baseline', None, None, None)))

ach_df = pd.read_csv(ACHILLES_PATH, index_col=0)

hgnc_df = pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt',
                      index_col=0, sep='\t')
hgnc_df['entrez_id'] = [str(int(x)) if not np.isnan(x) else np.nan for x in hgnc_df['entrez_id']]

mapping = IdentifierMapping('human_transcriptomics', hgnc_df)
ach_mset = TypedOmicsMeasurementSet(ach_df.index, ach_df.columns, ach_df.values, mapping)
ensembl_patt = re.compile('\([0-9]*\)')

ach_mset.column_names = [str(ensembl_patt.findall(k)[0].replace('(', '').replace(')', '')) for k in
                         ach_mset.column_names]
conv_dict = mapping.get_id_table(ach_mset.column_names, 'entrez_id').set_index('entrez_id')['ensembl_gene_id'].to_dict()
ach_mset.column_names = [conv_dict[k] if k in conv_dict else np.nan for k in ach_mset.column_names]
ach_mset.drop(columns=ach_mset.data.columns[~ach_mset.data.columns.isin(model_genes)])
ach_mset.transform(lambda x: x.fillna(0))


## ESSENTIAL KOs for growth

def get_state_dict_from_ko(ko):
    d = {k: True for k in cobamp_model.gpr.get_genes()}
    for k in ko:
        d[k] = False
    return d


state_dicts = {gko: get_state_dict_from_ko([gko]) for gko in model_genes}
reaction_states = {gko: {i: j for i, j in {k: cobamp_model.gpr.eval_gpr(ind, state_dicts[gko])
                                           for ind, k in enumerate(cobamp_model.reaction_names)}.items() if not j}
                   for gko in state_dicts.keys()}

reaction_states_inactive = {k: tuple(frozenset({i for i, j in v.items() if j is not None})) for k, v in
                            reaction_states.items()}

rko_sets = {}
for k, v in reaction_states_inactive.items():
    if v not in rko_sets.keys():
        rko_sets[v] = {k}
    else:
        rko_sets[v] |= {k}

kos_to_test = [l for l in list(rko_sets.keys()) if len(l) > 0]



def batch_ko_optimize(model, kos, context, objective, objective_sense, threads, override_bounds=None, ppreturn=True):
    with model as m:
        if override_bounds is not None:
            print('\toverriding bounds...')
            for k, v in override_bounds.items():
                if k in model.reaction_names:
                    m.set_reaction_bounds(k, lb=v[0], ub=v[1])
        else:
            for crx, truth in context.items():
                if not truth: m.set_reaction_bounds(crx, lb=0, ub=0)

        sol = m.optimize(objective, objective_sense)
        ofv = sol.objective_value()
        if sol.status() == 'optimal' and ofv > 0:
            gko_to_rko = [{model.decode_index(rko, 'reaction'): (0, 0) for rko in ko_set} for ko_set in kos]

            bopt = BatchOptimizer(model.model, threads=min(len(gko_to_rko), threads, MP_THREADS))

            opt_result = bopt.batch_optimize(gko_to_rko,
                                             [{model.decode_index(rko, 'reaction'): v for rko, v in
                                               objective.items()}] * len(gko_to_rko),
                                             [objective_sense] * len(gko_to_rko))

            return [k.objective_value() / ofv if (k.status() == 'optimal') and ofv > 0 else 0 for k in
                    opt_result] if ppreturn else opt_result
        else:
            print('\tModel objective is 0... skipping')
            return []


# gapfill_pkl = read_pickle(os.path.join(CS_MODEL_DF_FOLDER, 'cs_models_lt2_prot_fba_gapfill.csv'))

model_dicts = model_df.T.to_dict()


corr_coef_dict = {}
for mkey, mdict in model_dicts.items():
    print(mkey)
    context_dict = {k: v for k, v in mdict.items()}
    #context_dict.update({r.id: True for r in model.boundary})
    bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context=context_dict,
                              objective={'biomass_human': 1}, objective_sense=False, threads=12)

    opt_res = dict(zip(kos_to_test, bkopt))
    of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in reaction_states_inactive.items()}
    ddf = ach_mset.data.reindex(index=['ACH-000019'], columns=of_res).append(pd.Series(of_res, name='biomass')).T
    corr_coef_dict[mkey] = (matthews_corrcoef((ddf['biomass'] < 0.01), (ddf['ACH-000019'] < -0.6)), ddf)
    print('\t', corr_coef_dict[mkey][0])

corr_coef_params = {}
for k, tup in corr_coef_dict.items():
    ddf = tup[1]
    corr_coef_params[k] = {}
    for i,lv in enumerate([-0.5, -0.6, -0.7, -0.8, -0.9]):
        for j,bv in enumerate([0.001, 0.999]):
            corr_coef_params[k][(i,j)] = matthews_corrcoef((ddf['biomass'] < bv), (ddf['ACH-000019'] < lv))

corr_coef_params = read_pickle(BOF_PAR_ESS_PATH)
score_dfs = pd.DataFrame.from_dict(corr_coef_params).T.reset_index()
score_dfs[['thres','globalmin', 'globalmax', 'local']] = \
    pd.DataFrame(score_dfs['level_1'].apply(lambda x: x.split('_') if isinstance(x, str) else [np.nan]*4).to_list())
score_dfs = score_dfs.drop(columns='level_1').rename(columns={'level_0': 'algorithm', 'level_2': 'int_function'})

score_df_melt = pd.melt(score_dfs, id_vars=['thres','globalmin', 'globalmax', 'local']+['algorithm','int_function'])\
    .rename(columns={'variable_0': 'essential_cutoff', 'variable_1': 'biomass_cutoff'})



import matplotlib.pyplot as plt
import seaborn as sns
pickle_object(corr_coef_dict, BOF_ESS_PATH)
pickle_object(corr_coef_params, BOF_PAR_ESS_PATH)

task_list = JSONTaskIO().read_task(TASKS_PATH)
for t in task_list: t.mandatory_activity = []

### RUN THIS TO EVALUATE TASKS

task_model = model.copy()
# task_results = evaluate_tasks(task_list, task_model, model_dicts, 12, 'CPLEX')
mkeys, mdicts = zip(*model_dicts.items())
mcontexts = [[k for k, v in m.items() if not v] for m in mdicts]

threads = 12
solver = 'CPLEX'
nthreads = min(len(task_list), threads if MP_THREADS >= threads > 0 else MP_THREADS)

all_reactions = set([r.id for r in task_model.reactions])

task_eval_results = {}

best_score_model = {best_score_id: model_dicts[best_score_id]}
for mkey, result in best_score_model.items():
    # using with statements to change the COBRA model temporarily
    # this is done to knock-out reaction not appearing the FASTCORE result
    with task_model as context_specific_model:
        protected = set([k for k, v in result.items() if v])  # get reactions included in the sample-specific model
        to_remove = all_reactions - protected  # get reactions except the protected ones
        for rid in to_remove:
            context_specific_model.reactions.get_by_id(rid).knock_out()  # knock-out reactions not in the model
        for boundary in task_model.boundary:
            boundary.knock_out()

        task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver=solver)
        task_names = task_eval.tasks

        task_ko_sols = {}
        for task_name in task_names:
            print('\t', task_name, mkey, sum(bkopt_sols.values()))
            task_eval.current_task = task_name
            bkopt = batch_ko_optimize(model=task_eval.model, kos=kos_to_test,
                                      context={}, objective={'biomass_human': 1},
                                      objective_sense=False, threads=12, ppreturn=False)
            bkopt_sols = dict(zip(kos_to_test, [sol.status() == 'optimal' for sol in bkopt]))
            task_ko_sols[task_name] = bkopt_sols
        task_eval_results[mkey] = task_ko_sols

pickle_object(task_eval_results,
              'projects/breast_mcf7/results/human1/reconstructions/cs_models_lt2_task_essentiality.pkl')

## RUN HERE TO PREDICT

corr_coef_dict = read_pickle(BOF_ESS_PATH)
corr_coef_params = pd.DataFrame(read_pickle(BOF_PAR_ESS_PATH)).T
task_eval_results = read_pickle(TASK_ESS_PATH)

task_dfs = {}
task_correlations = {}
for model_key, model_dict in task_eval_results.items():
    to_df = []
    for task, kod in model_dict.items():
        ind, val = zip(*kod.items())
        to_df.append(pd.Series(val, pd.Index(ind, dtype=tuple, tupleize_cols=False)))
    task_dfs[model_key] = pd.DataFrame(to_df).T

    test_mdl = task_dfs[model_key]
    mdl_genes_taskess = list(reaction_states_inactive.keys())
    gko_matrix = test_mdl.reindex(index=[reaction_states_inactive[k] for k in mdl_genes_taskess]).fillna(test_mdl.all())
    gko_matrix.index = mdl_genes_taskess

    from sklearn.linear_model import LogisticRegression

    X = (~gko_matrix).astype(int).loc[ach_mset.column_names, :]
    y = (ach_mset.data.loc['ACH-000019', :] > 0.5).astype(int)
    lm = LogisticRegression(C=100, max_iter=10000)
    lm.fit(X, y)
    task_correlations[model_key] = lm.score(X, y)
    print('Accuracy:', task_correlations[model_key])
    print('MCC:', matthews_corrcoef(y, lm.predict(X)))
    d1 = dict(enumerate(lm.coef_[0, :]))
    for k in sorted(d1, key=d1.get)[:10]:
        print(task_list[k], '=', d1[k])
    for k in sorted(d1, key=d1.get)[-10:]:
        print(task_list[k], '=', d1[k])
pd.Series(task_correlations)
