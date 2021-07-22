import os, sys

from itertools import chain

import pandas as pd


if __name__ == '__main__':
    arg = [argval.split('=')[1] for argval in sys.argv[1:] if '-config=' in argval]
    mp = '-no-mp' not in sys.argv[1:]

    if len(arg) > 0:
        ARGS_PATH = arg[0]
        with open(ARGS_PATH, 'r') as f:
            config_dict = {j[0].strip(): j[1].strip() for j in [k.split('=') for k in f.readlines()]}
    else:
        root_folder = 'projects/breast_mcf7'
        config_dict = {
            'MODEL_PATH': '',
            'CS_MODEL_DF_FOLDER': os.path.join(root_folder, 'results/human1/reconstructions/mcf7_comparison'),
            'CS_MODEL_NAMES': 'cs_models_all_combinations',
            'INDEX_COLUMNS': '0,1,2,3',
            'PROTECTED_REACTION_LIST': 'biomass_human,HMR_10023,HMR_10024',
            'SOURCES_TO_ADD': '/home/vvieira/cobamp:/home/vvieira/human_ts_models:/home/vvieira/troppo:'
                              '/home/vvieira/cobamp/src:/home/vvieira/troppo/src:/home/vvieira/human_ts_models'
        }


    if 'SOURCES_TO_ADD' in config_dict:
        sources_to_add = config_dict['SOURCES_TO_ADD'].split(':')
        for source in sources_to_add:
            print('Adding source-code folder:', source)
        sys.path.extend(sources_to_add)

    from cobra.io import read_sbml_model
    from models import get_human1_model
    from cobamp.wrappers import get_model_reader
    from troppo.methods_wrappers import GapfillWrapper

    from cobra.util.solver import solvers
    from cobra.core.configuration import Configuration

    Configuration.solver = solvers['cplex']


    print('Configuration:')
    for k, v in config_dict.items():
        print("\t" + k, '=', v)

    model = read_sbml_model(config_dict['MODEL_PATH']) \
        if 'MODEL_PATH' in config_dict.keys() != '' \
        else get_human1_model()
    cobamp_model = get_model_reader(model).to_cobamp_cbm('CPLEX')
    CS_MODEL_DF_FOLDER, CS_MODEL_NAMES = (config_dict[k] for k in ['CS_MODEL_DF_FOLDER', 'CS_MODEL_NAMES'])

    BIOMASS_RX = config_dict['BIOMASS_RX'] if 'BIOMASS_RX' in config_dict.keys() else 'biomass_human'
    BIOMASS_MT = config_dict['BIOMASS_MT'] if 'BIOMASS_MT' in config_dict.keys() else 'temp001s'

    MODEL_DF_PATH = os.path.join(CS_MODEL_DF_FOLDER, CS_MODEL_NAMES)
    model_df = pd.read_csv(MODEL_DF_PATH, index_col=list(map(int,config_dict['INDEX_COLUMNS'].split(','))))

    if 'PROTECTED_REACTION_LIST' in config_dict.keys():
        protected_reactions = config_dict['PROTECTED_REACTION_LIST'].split(',')
        model_df[protected_reactions] = True
    else:
        protected_reactions = []

    name_to_append = ''
    if 'OVERRIDE_COMBINATIONS' in config_dict:
        model_df_params = pd.read_csv(config_dict['OVERRIDE_COMBINATIONS'], index_col=0)
        models_to_reconstruct = [tuple(k) for k in model_df_params.values.tolist()]
        model_df = model_df.loc[models_to_reconstruct]
        name_to_append = '.'.join(config_dict['OVERRIDE_COMBINATIONS'].split('/')[-1].split('.')[:-1])

    ## gapfill for growth under open medium bounds
    gw = GapfillWrapper(model)

    exchanges = set(chain(*cobamp_model.get_boundary_reactions().values()))
    gapfill_results = {}

    from troppo.methods.gapfill.efm import DEFAULT_CONFIG
    DEFAULT_CONFIG['TIMELIMIT'] = 600

    for i, row in model_df.T.items():
        activity = row.to_dict()
        if row.sum() <= len(protected_reactions):
            print('Model',i,'has no reactions.')
        elif row.sum() <= int(0.01 * len(cobamp_model.reaction_names)):
            print('Model',i,'is too small to gapfill within a reasonable amount of time')
        else:
            for k in exchanges:
                activity[k] = True
            inactive_rxs = set([k for k,v in activity.items() if not v])
            gapfill_results[i] = activity
            with cobamp_model as m:
                for rx in inactive_rxs:
                    m.set_reaction_bounds(rx, lb=0, ub=0)
                sol = m.optimize({BIOMASS_RX: 1}, False)
            if sol.status() != 'optimal' or sol.objective_value() < 1e-3:
                print('Model',i,'did not pass:',sol)
                gapfill_sol = gw.run(avbl_fluxes=list(inactive_rxs),
                                     algorithm='efm',
                                     ls_override={'produced':[BIOMASS_MT]},
                                     kshproperties=DEFAULT_CONFIG)
                if len(gapfill_sol) > 0:
                    for k in gapfill_sol[0]:
                        activity[k] = True
                    with cobamp_model as m:
                        for rx in [k for k,v in activity.items() if not v]:
                            m.set_reaction_bounds(rx, lb=0, ub=0)
                        gsol = m.optimize({BIOMASS_RX: 1}, False)
                        print('Solution after gapfilling:',gsol)
                else:
                    del gapfill_results[i]
            else:
                print('Model',i,'passed:',sol)
    biomass_gapfill_models = pd.DataFrame(gapfill_results).T

    MODEL_DF_PATH_BIOMASS_GAPFILL = os.path.join(CS_MODEL_DF_FOLDER, CS_MODEL_NAMES +
                                                 '_biomass_gapfill_'+name_to_append+'.csv')
    biomass_gapfill_models.to_csv(MODEL_DF_PATH_BIOMASS_GAPFILL)
