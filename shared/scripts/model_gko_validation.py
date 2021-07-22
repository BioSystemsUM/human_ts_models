import sys, os
from os.path import join as joinpath
from os import listdir

if __name__ == '__main__':
    arg = [argval.split('=')[1] for argval in sys.argv[1:] if '-config=' in argval]

    if len(arg) > 0:
        ARGS_PATH = arg[0]
        with open(ARGS_PATH, 'r') as f:
            config_dict = {j[0].strip(): j[1].strip() for j in [k.split('=') for k in f.readlines()]}
    else:
        root_folder = 'projects/ccle_models'
        config_dict = {
            'MODEL_PATH': 'shared/models/Human-GEM_latest_consistent.xml',
            'CS_MODEL_DF_FOLDER': os.path.join(root_folder, 'results/human1/reconstructions/ccle_mcf7_params'),
            'CS_MODEL_NAMES': 'best_params_run_mapall.csv_biomass_gapfill_gapfill_map',
            'OVERRIDE_COMBINATIONS': 'projects/ccle_models/configs/run_maps/CCLE_expression_best_scores_run_map_koval0.csv',
            'INDEX_COLUMNS': '0,1,2',
            'SOURCES_TO_ADD': '/home/vvieira/cobamp:/home/vvieira/human_ts_models:/home/vvieira/troppo:'
                              '/home/vvieira/cobamp/src:/home/vvieira/troppo/src:/home/vvieira/human_ts_models',
            'GROWTH_REACTION_NAME': 'biomass_human',
            'NTHREADS':'12',
            'SIM_FRAMEWORK': 'cobamp'
        }

    if 'SOURCES_TO_ADD' in config_dict:
        sources_to_add = config_dict['SOURCES_TO_ADD'].split(':')
        for source in sources_to_add:
            print('Adding source-code folder:', source)
        sys.path.extend(sources_to_add)

    import pandas as pd

    import cobra
    from cobra.io import read_sbml_model
    from cobra.util.solver import solvers
    cobra.Configuration.solver = solvers['cplex']

    from cobamp.wrappers import get_model_reader, ConstraintBasedModelSimulator
    from cobamp.wrappers.cobamp import cobamp_simulate, cobamp_simulation_result_function, cobamp_fba
    from models import get_human1_model
    from mewpy.problems import AbstractKOProblem
    from mewpy.utils.constants import ModelConstants
    from mewpy.utils.process import RayEvaluator
    from mewpy.simulation import set_default_solver

    set_default_solver('cplex')

    IND_COLS = list(map(int,config_dict['INDEX_COLUMNS'].split(',')))
    THREADS = int(config_dict['NTHREADS'])
    SIM_FRAMEWORK = config_dict['SIM_FRAMEWORK'] if 'SIM_FRAMEWORK' in config_dict else 'mewpy'

    model_df = pd.concat([pd.read_csv(joinpath(config_dict['CS_MODEL_DF_FOLDER'], fn), index_col=IND_COLS)
     for fn in [l for l in listdir(config_dict['CS_MODEL_DF_FOLDER']) if l.startswith(config_dict['CS_MODEL_NAMES'])]], axis=0)

    name_to_append = ''
    if 'OVERRIDE_COMBINATIONS' in config_dict:
        model_df_params = pd.read_csv(config_dict['OVERRIDE_COMBINATIONS'], index_col=0)
        models_to_reconstruct = [tuple(k) for k in model_df_params.values.tolist()]
        model_df = model_df.loc[models_to_reconstruct]
        name_to_append = '.'.join(config_dict['OVERRIDE_COMBINATIONS'].split('/')[-1].split('.')[:-1])

    OUTPUT_NAME = config_dict['CS_MODEL_NAMES'] + name_to_append + '_essentiality.csv'
    model = read_sbml_model(config_dict['MODEL_PATH']) if config_dict['MODEL_PATH'] != '' else get_human1_model()
    cobamp_model = get_model_reader(model, ttg_ratio=9999).to_cobamp_cbm('CPLEX')

    model_genes = cobamp_model.gpr.get_genes()

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

    ModelConstants.RESET_SOLVER = False

    class KOProblem(AbstractKOProblem):
        def __init__(self, model, fevaluation=None, **kwargs):
            super(KOProblem, self).__init__(
                model, fevaluation=fevaluation, **kwargs)

            self.sim_args = {}
            self.sim_args['objective'] = kwargs.get('objective', None)
            self.sim_args['method'] = kwargs.get('method', 'FBA')
            self.sim_args['maximize'] = kwargs.get('maximize', True)
            self.sim_args['reference'] = kwargs.get('reference', None)

        def _build_target_list(self):
            return []

        def evaluator(self, candidate_list, arg):
            result_list = []
            for const in candidate_list:
                res = self.simulator.simulate(constraints=const, **self.sim_args)
                result_list.append(res)
            return result_list


    # do not reset the solver


    # build a problem

    # lists of constraints to evaluate
    # multi processing evaluation
    if SIM_FRAMEWORK == 'mewpy':
        problem = KOProblem(model, fevaluation=[])
        rayeval = RayEvaluator(problem, THREADS)
    else:
        problem = None
        rayeval = None

    series_list = []
    for mname,row in model_df.iterrows():
        print('Model =',mname)
        model_inactives = {k: (0, 0) for k, v in row.to_dict().items() if not v}
        bound_changes = [{k:(0,0) for k in kos} for kos in kos_to_test]

        if SIM_FRAMEWORK == 'cobamp':
            with cobamp_model as cm:
                for rx in model_inactives: cm.set_reaction_bounds(rx, lb=0, ub=0)

                sol = cm.optimize({config_dict['GROWTH_REACTION_NAME']: 1}, False)
                default_sol_value = int(sol.status() == 'optimal') * sol.objective_value()

                simulator = ConstraintBasedModelSimulator(cm, cobamp_simulate, cobamp_simulation_result_function)
                ko_results = simulator.batch_simulate(cobamp_fba, bound_changes,
                                                      [{config_dict['GROWTH_REACTION_NAME']: 1}],
                                                      [False], None, THREADS)

                ko_sim_status, ko_sim_ofv, _ = zip(*ko_results)
        elif SIM_FRAMEWORK == 'mewpy':
            with model as m:
                for k,v in model_inactives.items(): m.reactions.get_by_id(k).bounds = v
                m.objective = config_dict['GROWTH_REACTION_NAME']
                sol = m.optimize()
                default_sol_value = int(sol.status == 'optimal') * sol.objective_value
            btt = [{k: (0, 0) for k in kos} for kos in kos_to_test]
            for d in btt:
                d.update(model_inactives)
            res = rayeval.evaluate(btt, None)
            ko_sim_status, ko_sim_ofv = zip(*[(str(sol.status).lower() == 'optimal',
                                               sol.fluxes[config_dict['GROWTH_REACTION_NAME']]) for sol in res])

        else:
            raise Exception('No valid simulation framework was found matching '+SIM_FRAMEWORK)

        info = [[list(map(int,ko_sim_status)), ko_sim_ofv],['status','objective_value']]
        single_model_ko_df = pd.DataFrame([pd.Series(l, index=kos_to_test,name=n) for l,n in zip(*info)]).T
        single_model_ko_series = single_model_ko_df['status']*single_model_ko_df['objective_value'].fillna(0)

        gko_series = single_model_ko_series.reindex([reaction_states_inactive[k] for k in model_genes])\
            .fillna(default_sol_value)
        gko_series.index = model_genes
        gko_series.name = mname

        series_list.append(gko_series)

    growth_per_ko_df = pd.DataFrame(series_list).to_csv(joinpath(config_dict['CS_MODEL_DF_FOLDER'],OUTPUT_NAME))