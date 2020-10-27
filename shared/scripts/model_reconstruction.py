import os, sys
import pandas as pd
import re
import resource

if __name__ == '__main__':
    arg = [argval.split(' =')[1] for argval in sys.argv[1:] if '-config=' in argval]
    mp = '-no-mp' not in sys.argv[1:]

    if len(arg) > 0:
        ARGS_PATH = arg[0]
        with open(ARGS_PATH, 'r') as f:
            config_dict = {j[0].strip(): j[1].strip() for j in [k.split('=') for k in f.readlines()]}
    else:
        root_folder = 'projects/breast_mcf7'
        config_dict = {
            'ROOT_FOLDER': 'projects/breast_mcf7',
            'MODEL_PATH': '',
            'DATA_PATH': os.path.join(root_folder,
                                      'data/ccle/DepMap Public 20Q1/CCLE_expression_ACH-000019_scores.csv'),
            'SAMPLE_INFO': os.path.join(root_folder, 'data/ccle/DepMap Public 20Q1/sample_info_v2.csv'),
            'CS_MODEL_DF_FOLDER': os.path.join(root_folder, 'results/human1/reconstructions/mcf7_comparison'),
            'CS_MODEL_NAMES': 'cs_models_all_combinations',
            'GENE_NAME_GRAB_PATTERN': 'ENSG[0-9]*',
            'PROTECTED_REACTION_LIST': 'biomass_human,HMR_10023,HMR_10024',
            'INDEX_COLUMNS': '0,1,2,3',
            'NTHREADS': '12',
            'SOURCES_TO_ADD': '/home/vvieira/cobamp:/home/vvieira/human_ts_models:/home/vvieira/troppo:'
                              '/home/vvieira/cobamp/src:/home/vvieira/troppo/src:/home/vvieira/human_ts_models'
        }

    print('Configuration:')
    for k, v in config_dict.items():
        print("\t" + k, '=', v)

    if 'SOURCES_TO_ADD' in config_dict:
        sources_to_add = config_dict['SOURCES_TO_ADD'].split(':')
        for source in sources_to_add:
            print('Adding source-code folder:', source)
        sys.path.extend(sources_to_add)

    from cobra.io import read_sbml_model
    from cobra.flux_analysis import find_blocked_reactions
    from urllib.request import urlretrieve

    from cobamp.utilities.parallel import batch_run, cpu_count
    from troppo.methods_wrappers import ReconstructionWrapper
    from troppo.omics.core import OmicsMeasurementSet, OmicsContainer

    from troppo.methods_wrappers import integration_strategy_map
    from troppo.omics.integration import MINSUM, MINMAX

    from itertools import product

    if not os.path.exists(config_dict['CS_MODEL_DF_FOLDER']): os.makedirs(config_dict['CS_MODEL_DF_FOLDER'])


    def get_human1_model():
        path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/model/Human-GEM.xml')
        model = read_sbml_model(path)
        model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
        blocked = find_blocked_reactions(model)
        model.remove_reactions(blocked)
        model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
        return model

    ## create a regex pattern to identify ensembl names in a string
    ## Create an identifier mapping object
    ## This is needed only if you need to convert genes from one nomenclature to another
    print('Reading data...')
    exp_data = pd.read_csv(config_dict['DATA_PATH'], index_col=list(map(int,config_dict['INDEX_COLUMNS'].split(','))))

    # Create an omics measurement set object with the dataframe components
    omics_mset = OmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values)

    # Keep the ensembl gene ID only
    if 'GENE_NAME_GRAB_PATTERN' in config_dict.keys():
        print('\t','Grabbing gene names from pattern')
        ensembl_patt = re.compile(config_dict['GENE_NAME_GRAB_PATTERN'])
        omics_mset.column_names = [ensembl_patt.findall(k)[0] for k in omics_mset.column_names]


    data_dicts = {'_'.join(map(str, k)) if isinstance(k, tuple) else k: v for k, v in omics_mset.data.T.to_dict().items()}

    del exp_data, omics_mset

    protected = list(map(lambda x: x.strip(),config_dict['PROTECTED_REACTION_LIST'].split(','))) \
        if ('PROTECTED_REACTION_LIST' in config_dict.keys()) else []

    #### parameter setup

    params = {
              'algorithms':['tinit','fastcore'],
              'strategies': {
                  'tinit': integration_strategy_map['adjusted_score'](protected),
                  'fastcore': integration_strategy_map['default_core'](0, protected)
              },
              'functions':{'minmax': MINMAX, 'minsum': MINSUM},
    }

    if 'OVERRIDE_COMBINATIONS' in config_dict:
        model_df_params = pd.read_csv(config_dict['OVERRIDE_COMBINATIONS'], index_col=0)
        models_to_reconstruct = [tuple(k) for k in model_df_params.values.tolist()]
    else:
        models_to_reconstruct = list(product(params['algorithms'], data_dicts.keys(), params['functions'].keys()))

    print('Reconstructing',len(models_to_reconstruct),'models...')

    for k in params['algorithms']:
        batch = [j for j in models_to_reconstruct if j[0] == k]
        print(len(batch),'models to reconstruct with',k)


    model = read_sbml_model(config_dict['MODEL_PATH']) \
        if 'MODEL_PATH' in config_dict.keys() != '' \
        else get_human1_model()


    rw = ReconstructionWrapper(model, ttg_ratio=9999)
    params['rw'] = rw

    NTHREADS = int(cpu_count() if 'NTHREADS' not in config_dict.keys() else config_dict['NTHREADS'])

    print('Using',NTHREADS,'threads.')
    print('Available CPUs:',cpu_count())

    def reconstruct_model(options, params):
        print('\tResource statistics before reconstruction with',options[0],':', resource.getrusage(resource.RUSAGE_SELF))
        alg, d, func = options
        data_dict, aofunc = d, params['functions'][func]

        oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x',
                                   data=data_dict, nomenclature='custom')

        params['rw'].run_from_omics(omics_data=oc_sample, algorithm=alg, and_or_funcs=aofunc,
                          integration_strategy=params['strategies'][alg], solver='CPLEX', raise_errors=False)

    safe_threads = {'tinit': max(1, NTHREADS // 16), 'fastcore': NTHREADS}

    print('Resource statistics before reconstruction:',resource.getrusage(resource.RUSAGE_SELF))
    reconstructions = {}

    for k, v in safe_threads.items():
        batch = [j for j in models_to_reconstruct if j[0] == k]

        if len(batch) > 0:
            alg, dd, intf = zip(*batch)
            ddicts = [data_dicts[i] for i in dd]
            batch_args = list(zip(alg, ddicts, intf))

            if mp:
                print('\tResource statistics before reconstruction with',k,':', resource.getrusage(resource.RUSAGE_SELF))
                reconstructions.update(dict(zip(batch,batch_run(reconstruct_model, batch_args, params, threads=v))))
            else:
                for barg in batch_args:
                    reconstructions[tuple(batch)] = reconstruct_model(barg, params)

    print('Writing models')
    pd.DataFrame.from_dict(reconstructions, orient='index').\
        to_csv(os.path.join(config_dict['CS_MODEL_DF_FOLDER'], config_dict['CS_MODEL_NAMES'] + '.csv'))