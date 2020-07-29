import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])


import os
import pandas as pd
import numpy as np
import re
from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions, pfba
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.utilities.file_io import pickle_object
from cobamp.utilities.parallel import batch_run
from cobamp.wrappers.external_wrappers import get_model_reader
from projects.breast_mcf7.scripts.omics import growth_media
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer
from troppo.methods_wrappers import GapfillWrapper


from itertools import chain

ROOT_FOLDER = 'projects/breast_mcf7'
DATA_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/CCLE_expression_multiple_apprx.csv')
SAMPLE_INFO = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/sample_info_v2.csv')
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
NTHREADS = 12

## read the metabolic model and simplify it
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

cobamp_model = get_model_reader(model).to_cobamp_cbm('CPLEX')


model_df = pd.read_csv(os.path.join(CS_MODEL_DF_FOLDER,'cs_models_all_combinations.csv'), index_col=[0,1,2])
# model_df[['HMR_10023', 'HMR_10024', 'biomass_human']] = True

## gapfill for growth under open medium bounds
gapfill_results = {}
gw = GapfillWrapper(model)

exchanges = set(chain(*cobamp_model.get_boundary_reactions().values()))


set(cobamp_model.bounds)
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

pd.DataFrame(gapfill_results).T.to_csv(
    os.path.join(CS_MODEL_DF_FOLDER, 'cs_models_all_combinations_efm_gapfill_open_media.csv'))

biomass_gapfill_models = pd.read_csv(
    os.path.join(CS_MODEL_DF_FOLDER, 'cs_models_all_combinations_efm_gapfill_open_media.csv'))

#
media_metabolites = growth_media.get_medium_metabolites('MEM', 'iHuman_id')
# extra_metabolites = {'m02631s', 'm02039s', 'm02394s', 'm01822s',
#                      'm02982s', 'm01938s', 'm01401s', 'm01330s',
#                      'm02394s', 'm01361s', 'm01822s', 'm01361s',
#                      'm01570s', 'm03147s', 'm01385s', 'm01629s'}
#
# {k:model.metabolites.get_by_id(k).name for k in extra_metabolites}




extra_metabs = {'m10005s', 'm02560s'}
media_to_remove = {'m01628s'}
media_metabolite_ids = ((set(media_metabolites['iHuman_id'].dropna().unique())) - media_to_remove) & \
                       set([m.id for m in model.metabolites]) | extra_metabs

cobamp_model = get_model_reader(model).to_cobamp_cbm('CPLEX')
cobamp_model_boundrx = cobamp_model.get_boundary_reactions()
media_boundrx = {k: cobamp_model_boundrx[k][0] for k in media_metabolite_ids if k in cobamp_model_boundrx.keys()}
cobamp_model_boundrx_rev = {v[0]:k for k,v in cobamp_model_boundrx.items()}



# ### check fatty acids to add
# extra_fas = set()
# for fametab in model.reactions.get_by_id('HMR_10033').metabolites:
#     tentative = fametab.id[:-1]+'s'
#     if tentative in cobamp_model_boundrx:
#         extra_fas |= {tentative}
#
# media_boundrx.update({m:cobamp_model_boundrx[m][0] for m in extra_fas})

outflows = set(chain(*cobamp_model_boundrx.values())) - set(media_boundrx.values())
efm_gapfill = {}

irrev_cobamp_model, irrev_cobamp_model_mapping = cobamp_model.make_irreversible()
input_drains = [k+'_flux_backwards' for k in outflows if k+'_flux_backwards' in irrev_cobamp_model.reaction_names]

from troppo.methods_wrappers import GapfillWrapper
gw = GapfillWrapper(irrev_cobamp_model)

with irrev_cobamp_model as m:
    orig_sol = irrev_cobamp_model.optimize({'biomass_human': 1})
    for k in input_drains: irrev_cobamp_model.set_reaction_bounds(k, lb=0, ub=0)
    media_sol = irrev_cobamp_model.optimize({'biomass_human': 1})

from cobamp.algorithms.kshortest import *

ksh_cfg = KShortestProperties()
ksh_cfg[K_SHORTEST_MPROPERTY_METHOD] = K_SHORTEST_METHOD_ITERATE
ksh_cfg[K_SHORTEST_OPROPERTY_BIG_M_CONSTRAINTS] = True
ksh_cfg[K_SHORTEST_OPROPERTY_BIG_M_VALUE] = 1e6
ksh_cfg[K_SHORTEST_OPROPERTY_FORCE_NON_CANCELLATION] = True
ksh_cfg[K_SHORTEST_OPROPERTY_MAXSOLUTIONS] = 1

from troppo.methods_wrappers import GapfillWrapper

external_metabolites = set(m for m in model.metabolites if m.compartment == 's')
external_metabolites -= set(m for m in external_metabolites if 'pool' in m.name and 'NEFA' not in m.name)
real_compounds = [m for m in external_metabolites if m.formula is not None and len(m.formula) > 0]



## TODO: remove drains for compounds without formula (to eliminate metabolite pools, etc...)
efm_media_gapfill_result = {}
for i, mrow in biomass_gapfill_models.iterrows():
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

    # gapfill!
    gapfill_res = gw.run(avbl_fluxes=input_drains, algorithm='efm',
                         ls_override={'produced':['temp001s']})

    efm_media_gapfill_result[i] = gapfill_res[0] if len(gapfill_res) > 0 else []
    print(i, gapfill_res)

