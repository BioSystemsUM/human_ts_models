import sys, os

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])


from cobamp.wrappers.external_wrappers import get_model_reader

import pandas as pd
import json

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.utilities.file_io import pickle_object, read_pickle

AVAILABLE_THREADS = 12

ROOT_FOLDER = 'projects/breast_mcf7'
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
FLX_PATH = 'projects/breast_mcf7/data/katzir2019/fluxes.csv'
RX_ASSOC_URL = 'https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/master/data/annotation/humanGEMRxnAssoc.JSON'
MODEL_DF_NAME = 'cs_models_all_combinations'
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/mcf7_comparison')
MODEL_DF_PATH = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'.csv')
MODEL_DF_PATH_BIOMASS_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'_biomass.csv')
MODEL_DF_FULL_GAPFILL_MEDIA = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'media_reactions.json')

NTHREADS = 12
if not os.path.exists(MODEL_PATH):
    path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/modelFiles/xml/HumanGEM.xml')
    model = read_sbml_model(path)
    write_sbml_model(model, MODEL_PATH)
else:
    model = read_sbml_model(MODEL_PATH)
model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
blocked = find_blocked_reactions(model)
model.remove_reactions(blocked)
model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)

model_df = pd.read_csv(MODEL_DF_PATH_BIOMASS_GAPFILL, index_col=[0, 1, 2])
cobamp_model = get_model_reader(model).to_cobamp_cbm('CPLEX')

flux_data = pd.read_csv(FLX_PATH)
flux_data.iloc[:,0] = list(map(
    lambda x: x.replace('_DASH_','_').replace('_LPAREN(','(').replace(')RPAREN_',')').replace('(','[').replace(')',']'),
    flux_data.iloc[:,0]))
flux_data = flux_data.replace('EX_glu_L[e]','EX_gln_L[e]')
rx_assoc_df = pd.read_json(RX_ASSOC_URL).set_index('rxns')

flux_data.iloc[:,0] = flux_data.iloc[:,0].replace(
    {v:k for k,v in
     rx_assoc_df[rx_assoc_df['rxnRecon3DID'].isin(flux_data.iloc[:,0])]['rxnRecon3DID'].to_dict().items()})
flux_data = flux_data.set_index(flux_data.columns[0]).T
flux_data_std = (flux_data - flux_data.mean())/flux_data.std()

with open(MODEL_DF_FULL_GAPFILL_MEDIA, 'r') as f:
    dct = dict([(tuple(k[0]),k[1]) for k in json.JSONDecoder().decode(f.read())])


scenarios = {
    'A': {},
    'B': {'HMR_9063': (0, 1000)},
    'C': {'HMR_6916': (0, 0)}
}

def compare_simulation(expected, predicted):
    predicted_std = (predicted - predicted.mean())/predicted.std()
    predicted_sign = (predicted_std / abs(predicted_std)).fillna(0)
    expected_sign = (expected/abs(expected)).fillna(0)
    return (expected_sign == predicted_sign).all()


from cobamp.utilities.parallel import batch_run
from cobra.flux_analysis.parsimonious import pfba


def simulate_model(name_row_tuple, params):
    ind, model_row = name_row_tuple
    scenarios, model, dct, flux_data = (params[k] for k in ['scenarios', 'model', 'dct', 'flux_data'])
    print('Model =',ind)
    sim_results = {}
    sim_results_bounded = {}
    boundary = {r.id for r in model.boundary}
    inactives = {k for k,v in model_row.to_dict().items() if not v}
    media_reactions = set([k.split('_flux_backwards')[0] for k in dct[ind]])
    with model as m:
        for rx in inactives - (boundary | {'HMR_10024', 'HMR_10023'}):
            m.reactions.get_by_id(rx).bounds = (0, 0)
        for sc_name, sc_bounds in scenarios.items():
            with m as scen_model:
                for rx, b in sc_bounds.items(): scen_model.reactions.get_by_id(rx).bounds = b
                try:
                    sol = pfba(scen_model).fluxes
                except:
                    sol = None
                sim_results[sc_name] = sol
                for rx in boundary - {'HMR_10024', 'HMR_10023'}: scen_model.reactions.get_by_id(rx).bounds = (0, 1000)
                for rx in media_reactions: scen_model.reactions.get_by_id(rx).bounds = (-1000, 1000)
                for rx, b in sc_bounds.items(): scen_model.reactions.get_by_id(rx).bounds = b
                try:
                    sol_b = pfba(scen_model).fluxes
                except:
                    print('Model',ind,'with media failed at',sc_name)
                    sol_b = None
                sim_results_bounded[sc_name] = sol_b
                #print(ind,sc_name, sol.fluxes['biomass_human'], sol_b.fluxes['biomass_human'])


        sim_results = pd.DataFrame(sim_results).loc[flux_data.columns, :].T
        sim_results_bounded = pd.DataFrame(sim_results_bounded).loc[flux_data.columns, :].T
        return sim_results, sim_results_bounded

inds, rows = list(zip(*model_df.iterrows()))
params = {'dct': dct, 'scenarios': scenarios, 'model': model, 'flux_data': flux_data}
results = dict(zip(inds,batch_run(simulate_model, list(zip(inds, rows)), params, threads=12)))


pickle_object(results, os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'_simulation_dicts.pkl'))

#sim_results, sim_results_bounded = read_pickle(os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'_simulation_dicts.pkl'))

results = read_pickle(os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'_simulation_dicts.pkl'))

def get_sign(df): return (df/df.abs()).fillna(0)

flux_data_sign = get_sign(flux_data)

sign_agreement = {}
for tested_model, result_df in results.items():
    sim, sim_bound = result_df
    sign_agreement[tested_model] = {'unbounded': (get_sign(sim) == flux_data_sign).sum(axis=0).to_dict(),
                                    'bounded': (get_sign(sim_bound.fillna(0)) == flux_data_sign).sum(axis=0).to_dict()}

adf = pd.DataFrame.from_dict(sign_agreement, orient='index').stack().to_frame()
flux_comparison_df = (pd.DataFrame(adf[0].values.tolist(), index=adf.index) > 1).sum(axis=1).sort_values()