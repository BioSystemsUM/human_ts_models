import sys
from itertools import chain

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

import os
import re
import pandas as pd
import numpy as np

from sklearn.metrics import matthews_corrcoef

from cobamp.core.optimization import BatchOptimizer
from cobamp.utilities.file_io import pickle_object
from cobamp.utilities.parallel import MP_THREADS
from cobamp.wrappers.external_wrappers import get_model_reader
from troppo.omics.core import IdentifierMapping, TypedOmicsMeasurementSet
from projects.breast_mcf7.scripts.human1.models import get_human1_model

ROOT_FOLDER = 'projects/breast_mcf7'
MODEL_DF_NAME = 'cs_models_all_combinations'
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/mcf7_comparison')
MODEL_DF_PATH = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'.csv')
MODEL_DF_PATH_BIOMASS_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'_biomass.csv')
MODEL_DF_PATH_FULL_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'biomass_media.csv')
MODEL_DF_FULL_GAPFILL_MEDIA = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'media_reactions.json')
ACHILLES_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/Achilles_gene_effect.csv'
ORIGINAL_MEDIA_PATH = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'original_media.txt')

BOF_ESS_FOLDER = os.path.join(CS_MODEL_DF_FOLDER,'essentiality')
BOF_ESS_PATH = os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'_bof_essentiality.pkl')
BOF_PAR_ESS_PATH = os.path.join(BOF_ESS_FOLDER,MODEL_DF_NAME+'_bof_essentiality_params.pkl')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
NTHREADS = 12


model = get_human1_model()
cobamp_model = get_model_reader(model, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
model_genes = cobamp_model.gpr.get_genes()

model_df = pd.read_csv(MODEL_DF_PATH_BIOMASS_GAPFILL, index_col=[0, 1, 2])

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
        for crx, truth in context.items():
            if not truth: m.set_reaction_bounds(crx, lb=0, ub=0)
        if override_bounds is not None:
            print('\toverriding bounds...')
            for k, v in override_bounds.items():
                if k in m.reaction_names:
                    m.set_reaction_bounds(k, lb=v[0], ub=v[1])


        sol = m.optimize(objective, objective_sense)
        ofv = sol.objective_value()
        if sol.status() == 'optimal' and ofv > 0:
            gko_to_rko = [{m.decode_index(rko, 'reaction'): (0, 0) for rko in ko_set} for ko_set in kos]

            bopt = BatchOptimizer(m.model, threads=min(len(gko_to_rko), threads, MP_THREADS))

            opt_result = bopt.batch_optimize(gko_to_rko,
                                             [{m.decode_index(rko, 'reaction'): v for rko, v in
                                               objective.items()}] * len(gko_to_rko),
                                             [objective_sense] * len(gko_to_rko))

            return [k.objective_value() / ofv if (k.status() == 'optimal') and ofv > 0 else 0 for k in
                    opt_result] if ppreturn else opt_result
        else:
            print('\tModel objective is 0... skipping')
            return []


# gapfill_pkl = read_pickle(os.path.join(CS_MODEL_DF_FOLDER, 'cs_models_lt2_prot_fba_gapfill.csv'))

model_dicts = model_df.T.to_dict()
from json import JSONDecoder
with open(MODEL_DF_FULL_GAPFILL_MEDIA, 'r') as f:
    model_dict_media = dict([(tuple(i),k) for i,k in JSONDecoder().decode(f.read())])

corr_coef_dict = {}
for mkey, mdict in model_dicts.items():
    def process_results(bkopt):
        opt_res = dict(zip(kos_to_test, bkopt))
        of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in reaction_states_inactive.items()}
        ddf = ach_mset.data.reindex(index=['ACH-000019'], columns=of_res).append(pd.Series(of_res, name='biomass')).T
        return matthews_corrcoef((ddf['biomass'] < 0.01), (ddf['ACH-000019'] < -0.6)), ddf

    print(mkey)
    context_dict = {k: v for k, v in mdict.items()}
    #context_dict.update({r.id: True for r in model.boundary})
    bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context=context_dict,
                              objective={'biomass_human': 1}, objective_sense=False, threads=12)
    res = process_results(bkopt)
    corr_coef_dict[tuple(list(mkey)+['no_media'])] = res
    print('\t without media:', res[0])

    if mkey in model_dict_media:
        exchanges = set(chain(*cobamp_model.get_boundary_reactions().values()))
        media_list = set([k.split('_flux_backwards')[0] for k in model_dict_media[mkey]])
        bound_override = {k:(0, 1000) for k in exchanges
                          - media_list - {'HMR_10024', 'HMR_10023', 'biomass_human'}}
        bound_override.update({k:(-1000, 1000) for k in media_list})

        bkopt_media = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context=context_dict,
                                        override_bounds=bound_override, objective={'biomass_human': 1},
                                        objective_sense=False, threads=12)
        res_med = process_results(bkopt_media)
        corr_coef_dict[tuple(list(mkey)+['media_gapfill'])] = res_med
        print('\t with media:', res_med[0])

pickle_object(corr_coef_dict, BOF_ESS_PATH)

corr_coef_params = {}
for k, tup in corr_coef_dict.items():
    ddf = tup[1]
    corr_coef_params[k] = {}
    for i,lv in enumerate([-0.5, -0.6, -0.7, -0.8, -0.9]):
        for j,bv in enumerate([0.001, 0.999]):
            corr_coef_params[k][(i,j)] = matthews_corrcoef((ddf['biomass'] < bv), (ddf['ACH-000019'] < lv))

pickle_object(corr_coef_params, BOF_PAR_ESS_PATH)