import os
import sys

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

import pandas as pd
import numpy as np

from cobamp.utilities.file_io import read_pickle

import matplotlib.pyplot as plt
import seaborn as sns

ROOT_FOLDER = 'projects/breast_mcf7'
MODEL_DF_NAME = 'cs_models_all_combinations'
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/mcf7_comparison')
MODEL_DF_PATH = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'.csv')

BOF_ESS_FOLDER = os.path.join(CS_MODEL_DF_FOLDER,'essentiality')
BOF_ESS_PATH = os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'_bof_essentiality.pkl')
BOF_PAR_ESS_PATH = os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'_bof_essentiality_params.pkl')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
NTHREADS = 12

corr_coef_params = read_pickle(BOF_PAR_ESS_PATH)
score_dfs = pd.DataFrame.from_dict(corr_coef_params).T.reset_index()
score_dfs[['thres','globalmin', 'globalmax', 'local']] = \
    pd.DataFrame(score_dfs['level_1'].apply(lambda x: x.split('_') if isinstance(x, str) else [np.nan]*4).to_list())
score_dfs = score_dfs.drop(columns='level_1').rename(columns={'level_0': 'algorithm', 'level_2': 'int_function','level_3': 'media_gapfill'})
idvs = ['thres','globalmin', 'globalmax', 'local']+['algorithm','int_function','media_gapfill']

score_df_melt = pd.melt(score_dfs, id_vars=idvs).rename(
    columns={'variable_0': 'essential_cutoff', 'variable_1': 'biomass_cutoff'})

score_dfs_no_id = score_dfs.drop(columns=idvs)
score_dfs_id = score_dfs[idvs].replace({'nan': np.nan})


blk_col = sns.colors.xkcd_rgb['black']
rcols = []
palettes = [
    list(map(lambda x: sns.colors.xkcd_rgb[x], ['seafoam green', 'light blue', 'blue', 'black'])), # thres
    sns.color_palette("Reds")[:5] + [blk_col], # global minimum
    [blk_col] + sns.color_palette("Blues")[:5], # global maximum
    [blk_col] + sns.color_palette("BuGn")[:5], # local threshold
    sns.color_palette("Set2")[:2] + [blk_col], # algorithm
    sns.color_palette("Set2")[2:4] + [blk_col],# integration_function
    sns.color_palette("Set2")[4:6] # media gapfill
]


for col, pal in zip(score_dfs_id.columns, palettes):
    series = score_dfs_id[col]
    unq = series.unique()
    print(col, unq)
    lut = dict(zip(unq, pal))
    rcols.append(series.map(lut))

row_colors = pd.concat(rcols,axis=1)
row_colors.iloc[score_dfs.shape[0]-1,:] = ['yellow']
row_colors.columns = ['Thresholding strategy', 'Unexpressed threshold', 'Expression threshold', 'Local threshold',
                      'Algorithm', 'Integration function', 'Growth media gapfill']

lethals_df = score_dfs_no_id.iloc[:,np.arange(0, 10, 2)]
lethals_df.columns = ['0.5','0.6','0.7','0.8','0.9']
g = sns.clustermap(lethals_df, row_colors=row_colors, yticklabels=False, col_cluster=False)
plt.gcf().set_size_inches(6, 8)
plt.tight_layout()
plt.savefig(os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'essentiality_par_lethals.pdf'))

decreased_df = score_dfs_no_id.iloc[:,np.arange(0, 10, 2)+1]
decreased_df.columns = ['0.5','0.6','0.7','0.8','0.9']
h = sns.clustermap(decreased_df, row_colors=row_colors, yticklabels=False, col_cluster=False)
plt.gcf().set_size_inches(6, 8)
plt.tight_layout()
plt.savefig(os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'essentiality_par_dec.pdf'))

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder

lrdf = score_dfs_id.fillna(-1).iloc[:score_dfs_id.shape[0]-1,:].replace({'nan': -1})
lrdf[['globalmin', 'globalmax', 'local']] = lrdf[['globalmin', 'globalmax', 'local']].astype(float).astype(int)

lr = LinearRegression()

ohe = OneHotEncoder()
ohe_mat = ohe.fit_transform(lrdf[['thres', 'algorithm', 'int_function','media_gapfill']]).toarray()
lrdf_fin = lrdf.drop(columns=['thres', 'algorithm', 'int_function','media_gapfill']).replace({-1: 0, 0: 0.1, 1: 0.25, 2:0.5, 3:0.75, 4:0.9})
lrdf_fin.columns = ['globalmin', 'globalmax', 'local']
lrdf_full = pd.concat([lrdf_fin,pd.DataFrame(ohe_mat, columns=ohe.get_feature_names())], axis=1)
lr.fit(lrdf_full, score_dfs_no_id.iloc[:score_dfs_id.shape[0]-1,np.arange(0, 10, 2)])


coef_srs = pd.DataFrame(lr.coef_, columns=lrdf_full.columns, index=lethals_df.columns)
coef_srs.columns = ['+ Unexpressed threshold', '+ Expression threshold', '+ Local threshold', 'Global thresholding approach',
                  'Local thresholding approach (one state)', 'Local thresholding approach (two state)',
                  'FASTCORE algorithm', 'tINIT algorithm', 'Integrate most expressed isoform',
                  'Integrate sum of isoform expression', 'Gapfill missing media metabolites']
coef_srs.T.plot(kind='barh')
plt.axvline(x=0, color='black')
plt.title = 'Feature importance'
plt.ylabel = 'Linear regression coefficient'
plt.xlabel = 'Features'
plt.tight_layout()
plt.gcf().set_size_inches(8, 4)
plt.savefig(os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'essentiality_par_importance.pdf'))
