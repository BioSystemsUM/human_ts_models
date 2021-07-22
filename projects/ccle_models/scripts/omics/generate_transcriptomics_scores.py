import os
import pandas as pd
import numpy as np

from models import get_human1_model
from shared.src.thresholding import THRESHOLDING_FUNCTIONS

EXPFILE_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_full.csv'
EXPBEST_PATH = 'projects/ccle_models/data/ccle/DepMap Public 20Q1/CCLE_expression_best_scores.csv'
EXPRUN_PATH = 'projects/ccle_models/data/ccle/DepMap Public 20Q1/CCLE_expression_best_scores_run_map.csv'

ACHILLES_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/Achilles_gene_effect.csv'

for path in [EXPBEST_PATH, EXPRUN_PATH]:
    if not os.path.dirname(path):
        os.makedirs(os.path.dirname(path))

ach_df = pd.read_csv(ACHILLES_PATH, usecols=[0])
exp_df_ccle = pd.read_csv(EXPFILE_PATH, index_col=0)
exp_df_ccle = exp_df_ccle.loc[exp_df_ccle.index.isin(ach_df.iloc[:,0]),:]

qvalues = [0.1, 0.25, 0.5, 0.75, 0.9]
quantiles = exp_df_ccle.quantile(qvalues)
global_thresholds = quantiles.T.apply(lambda x: x.mean())

x = pd.read_csv(os.path.join('projects/breast_mcf7/results/human1/reconstructions/mcf7_comparison',
                             'cs_models_all_combinations'+'_best_params.csv'), index_col=0)
def coerce_to_number_if_possible(x):
    try:
        if x == 'nan':
            return np.nan
        else:
            return int(x)
    except:
        return x

combinations = pd.DataFrame(x.iloc[:,1].apply(lambda x: x.split('_')).tolist()).applymap(coerce_to_number_if_possible)
maxexp = np.log(exp_df_ccle.max().max())

nstuples = {}
for sample_id in exp_df_ccle.index:
    print(sample_id)
    sample = exp_df_ccle.loc[sample_id,:]
    for i, row in combinations.iterrows():
        thr_func = row[0]
        gtli, gtui, lti = map(lambda x: int(float(x)), row[1:].fillna(0))
        name = '_'.join(map(str, [thr_func,gtli,gtui,lti]))
        gtl, gtu, lti = global_thresholds.iloc[gtli], global_thresholds.iloc[gtui], quantiles.iloc[lti,]
        nstuples[sample.name+'_'+name] = THRESHOLDING_FUNCTIONS[thr_func](sample, gtl, gtu, lti, maxexp)



import re
patt = re.compile('ENSG[0-9]*')
model = get_human1_model()
dataset_genes = [g.id for g in model.genes]
patt_find = [(x, patt.findall(str(x))) for x in exp_df_ccle.columns]
new_ind = [x for x,k in patt_find if (len(k) > 0) and (k[0] in dataset_genes)]

nstuples_human1_genes = {k:{i:d[i] for i in new_ind} for k,d in nstuples.items()}

del nstuples

scores = pd.DataFrame(nstuples_human1_genes)
scores.T.to_csv(EXPBEST_PATH)

scores = pd.read_csv(EXPBEST_PATH, index_col=0)

run_map_combs = combinations.copy()
run_map_combs.iloc[:,1:] = run_map_combs.iloc[:,1:].fillna(0).astype(float).astype(int)
run_map_combs.T.apply(lambda x:'_'.join(map(str, x))).to_list()

run_map = x.copy()
run_map.iloc[:,1] = run_map_combs.T.apply(lambda x:'_'.join(map(str, x))).to_list()
dfs = []
for k in exp_df_ccle.index:
    subdf = run_map.copy()
    subdf.iloc[:,1] = k + '_' + subdf.iloc[:,1]
    dfs.append(subdf)

run_map_df = pd.concat(dfs)
run_map_df.to_csv(EXPRUN_PATH)

## split combinations
import numpy as np
for i, df in enumerate(np.array_split(run_map_df, 50)):
    df.to_csv('projects/ccle_models/configs/run_maps/CCLE_expression_best_scores_run_map'+
                      str(i)+'.csv')

import pandas as pd
run_map_df = pd.read_csv('projects/ccle_models/configs/run_maps/CCLE_expression_best_scores_run_map0.csv',index_col=0)

import pandas as pd

run_map_df = pd.concat([pd.read_csv('projects/ccle_models/configs/run_maps/CCLE_expression_best_scores_run_map' +str(i)
                                    + '.csv',index_col=0) for i in range(50)])

## cluster runs
tinit_runs = run_map_df[run_map_df.iloc[:,0] == 'tinit']
for i, df in enumerate(np.array_split(tinit_runs, 50)):
    df.to_csv('projects/ccle_models/configs/run_maps/CCLE_expression_best_scores_run_map_tinit'+str(i)+'.csv')

fc_runs = run_map_df[run_map_df.iloc[:,0] == 'fastcore']
fc_runs.to_csv('projects/ccle_models/configs/run_maps/CCLE_expression_best_scores_run_map_fastcore.csv')
