import sys;


print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])


import os
import pandas as pd
import numpy as np
import re
from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.wrappers.external_wrappers import get_model_reader
from projects.breast_mcf7.scripts.omics import growth_media
from projects.breast_mcf7.scripts.human1.models import get_human1_model

from troppo.methods_wrappers import GapfillWrapper

import json

from itertools import chain

ROOT_FOLDER = 'projects/breast_mcf7'
MODEL_DF_NAME = 'cs_models_all_combinations'
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/mcf7_comparison')
MODEL_DF_PATH = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'.csv')
MODEL_DF_PATH_BIOMASS_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'_biomass.csv')
MODEL_DF_PATH_FULL_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'biomass_media.csv')
MODEL_DF_FULL_GAPFILL_MEDIA = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'media_reactions.json')
ORIGINAL_MEDIA_PATH = os.path.join(CS_MODEL_DF_FOLDER,MODEL_DF_NAME+'original_media.txt')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
NTHREADS = 12

## read the metabolic model and simplify it
model = get_human1_model()
cobamp_model = get_model_reader(model).to_cobamp_cbm('CPLEX')


model_df = pd.read_csv(MODEL_DF_PATH, index_col=[0,1,2])
# model_df[['HMR_10023', 'HMR_10024', 'biomass_human']] = True

## gapfill for growth under open medium bounds
gapfill_results = {}
gw = GapfillWrapper(model)

exchanges = set(chain(*cobamp_model.get_boundary_reactions().values()))


if not os.path.exists(MODEL_DF_PATH_BIOMASS_GAPFILL):
    for i, row in model_df.T.items():
        if i not in gapfill_results.keys():
            activity = row.to_dict()
            for k in exchanges:
                activity[k] = True
            inactive_rxs = set([k for k,v in activity.items() if not v])
            gapfill_results[i] = activity
            with cobamp_model as m:
                for rx in inactive_rxs:
                    m.set_reaction_bounds(rx, lb=0, ub=0)
                sol = m.optimize({'biomass_human': 1}, False)
            if sol.status() != 'optimal' or sol.objective_value() < 1e-3:
                print('Model',i,'did not pass:',sol)
                gapfill_sol = gw.run(avbl_fluxes=list(inactive_rxs), algorithm='efm', ls_override={'produced':['temp001s']})
                for k in gapfill_sol[0]:
                    activity[k] = True
                with cobamp_model as m:
                    for rx in [k for k,v in activity.items() if not v]:
                        m.set_reaction_bounds(rx, lb=0, ub=0)
                    gsol = m.optimize({'biomass_human': 1}, False)
                    print('Solution after gapfilling:',gsol)
            else:
                print('Model',i,'passed:',sol)
    biomass_gapfill_models = pd.DataFrame(gapfill_results).T
    biomass_gapfill_models.to_csv(MODEL_DF_PATH_BIOMASS_GAPFILL)
else:
    biomass_gapfill_models = pd.read_csv(MODEL_DF_PATH_BIOMASS_GAPFILL, index_col = [0, 1, 2])

#
media_metabolites = growth_media.get_medium_metabolites('MEM', 'iHuman_id')


extra_metabs = {'m10005s', 'm02560s'}
media_to_remove = {'m01628s'}
media_metabolite_ids = ((set(media_metabolites['iHuman_id'].dropna().unique())) - media_to_remove) & \
                       set([m.id for m in model.metabolites]) | extra_metabs

cobamp_model = get_model_reader(model).to_cobamp_cbm('CPLEX')
cobamp_model_boundrx = cobamp_model.get_boundary_reactions()
media_boundrx = {k: cobamp_model_boundrx[k][0] for k in media_metabolite_ids if k in cobamp_model_boundrx.keys()}
cobamp_model_boundrx_rev = {v[0]:k for k,v in cobamp_model_boundrx.items()}


outflows = set(chain(*cobamp_model_boundrx.values())) - set(media_boundrx.values())
efm_gapfill = {}

irrev_cobamp_model, irrev_cobamp_model_mapping = cobamp_model.make_irreversible()
input_drains = [k+'_flux_backwards' for k in outflows if k+'_flux_backwards' in irrev_cobamp_model.reaction_names]

gw = GapfillWrapper(irrev_cobamp_model)

with irrev_cobamp_model as m:
    orig_sol = irrev_cobamp_model.optimize({'biomass_human': 1})
    for k in input_drains: irrev_cobamp_model.set_reaction_bounds(k, lb=0, ub=0)
    media_sol = irrev_cobamp_model.optimize({'biomass_human': 1})


external_metabolites = set(m for m in model.metabolites if m.compartment == 's')
external_metabolites -= set(m for m in external_metabolites if 'pool' in m.name and 'NEFA' not in m.name)
real_compounds = [m for m in external_metabolites if m.formula is not None and len(m.formula) > 0]



## TODO: remove drains for compounds without formula (to eliminate metabolite pools, etc...)
if not os.path.exists(MODEL_DF_PATH_FULL_GAPFILL):
    efm_media_gapfill_result = {}
    for i, mrow in biomass_gapfill_models.iterrows():
        print(i)
        # mdict is any result from a context-specific model reconstruction (dict[str,bool])
        mdict = mrow.to_dict()

        # list of reactions to remove from the context-specific model
        to_remove = set([k for k,v in mdict.items() if not v]) - set(cobamp_model_boundrx_rev.keys())

        # generate a cobamp model
        gapfill_model = get_model_reader(model).to_cobamp_cbm('CPLEX')
        gapfill_model.remove_reactions(list(to_remove)) # remove reactions not in the context-specific model
        gapfill_model.set_reaction_bounds('HMR_10024', lb=0, ub=1000) # set biomass drain to irreversible
        # make an irreversible version of the model - this is to create split reactions so that only uptake reactions
        # are added as part of the gapfill approach
        irrev_cobamp_model, irrev_cobamp_model_mapping = gapfill_model.make_irreversible()

        # decide which reactions in the original model should be added. this is basically all boundary reactions except
        # those already defined in the growth medium. media_boundrx is a dict mapping external metabolite ids to their
        # respective boundary reactions
        exp_outflows = set(chain(*gapfill_model.get_boundary_reactions().values())) - set(media_boundrx.values())

        # select drains from exp_outflows but only the reverse split reaction (uptake)
        input_drains = [k + '_flux_backwards' for k in exp_outflows if
                        k + '_flux_backwards' in irrev_cobamp_model.reaction_names]

        # create a gapfill wrapper instance with the irreversible model
        gw = GapfillWrapper(irrev_cobamp_model)

        from troppo.methods.gapfill.efm import DEFAULT_CONFIG
        DEFAULT_CONFIG['BIG_M_VALUE'] = 1e6
        # gapfill!
        gapfill_res = gw.run(avbl_fluxes=input_drains, algorithm='efm',
                             ls_override={'produced':['temp001s']}, kshproperties=DEFAULT_CONFIG)

        for k in list(set(input_drains) - set(gapfill_res[0])):
            irrev_cobamp_model.set_reaction_bounds(k, lb=0, ub=1)

        irrev_cobamp_model.set_reaction_bounds('biomass_human', lb=1, ub=1000)

        final_medium = (set(gapfill_res[0]) | set(media_boundrx.values()))

        sol = irrev_cobamp_model.optimize({k:1 for k in set(input_drains) - final_medium}, minimize=True)
        srs = sol.to_series()[set(input_drains) - final_medium]
        final_medium_pfba = {k.split('_flux_backwards')[0] for k in srs[srs > 1e-8].index}


        efm_media_gapfill_result[i] = list(final_medium | final_medium_pfba)
        print('\t',len(efm_media_gapfill_result[i]))

    with open(MODEL_DF_FULL_GAPFILL_MEDIA, 'w') as f:
        string = json.JSONEncoder().encode(list({k: list(v) for k, v in efm_media_gapfill_result.items()}.items()))
        f.write(string)
else:
    with open(MODEL_DF_FULL_GAPFILL_MEDIA, 'r') as f:
        efm_media_gapfill_result = dict([(tuple(k),v) for k,v in json.JSONDecoder().decode(f.read())])
