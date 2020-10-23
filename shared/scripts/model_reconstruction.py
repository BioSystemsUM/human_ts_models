import os, sys
import pandas as pd
import re


arg = [arg.split('=')[1] for arg in sys.argv[1] if '-config=' in arg]
if len(arg) > 0:
    ARGS_PATH = arg[0]
    with open(ARGS_PATH, 'r') as f:
        config_dict = {j[0]:j[1] for j in [k.split('=') for k in f.readlines()]}
else:
    root_folder = 'projects/breast_mcf7'
    config_dict = {
        'ROOT_FOLDER' : 'projects/breast_mcf7',
        'MODEL_PATH' : '',
        'DATA_PATH' : os.path.join(root_folder, 'data/ccle/DepMap Public 20Q1/CCLE_expression_ACH-000019_scores.csv'),
        'SAMPLE_INFO' : os.path.join(root_folder, 'data/ccle/DepMap Public 20Q1/sample_info_v2.csv'),
        'CS_MODEL_DF_FOLDER' : os.path.join(root_folder, 'results/human1/reconstructions/mcf7_comparison'),
        'CS_MODEL_NAMES' : 'cs_models_all_combinations',
        'GENE_NAME_GRAB_PATTERN' : 'ENSG[0-9]*',
        'PROTECTED_REACTION_LIST' : 'biomass_human,HMR_10023,HMR_10024',
        'INDEX_COLUMNS' : '0,1,2,3',
        'NTHREADS' : '12',
        'SOURCES_TO_ADD' : '/home/vvieira/cobamp:/home/vvieira/human_ts_models:/home/vvieira/troppo:'
                           '/home/vvieira/cobamp/src:/home/vvieira/troppo/src:/home/vvieira/human_ts_models'
    }

if 'SOURCES_TO_ADD' in config_dict:
    sys.path.extend(config_dict['SOURCES_TO_ADD'].split(':'))

from cobra.io import read_sbml_model
from cobamp.utilities.parallel import batch_run, cpu_count
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer
from projects.breast_mcf7.scripts.human1.models import get_human1_model

from troppo.methods_wrappers import integration_strategy_map
from troppo.omics.integration import MINSUM, MINMAX

from itertools import product

if not os.path.exists(config_dict['CS_MODEL_DF_FOLDER']): os.makedirs(config_dict['CS_MODEL_DF_FOLDER'])

## create a regex pattern to identify ensembl names in a string
## Create an identifier mapping object
## This is needed only if you need to convert genes from one nomenclature to another
mapping = IdentifierMapping('human_transcriptomics',
                            pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt',
                                        index_col=0, sep='\t'))


exp_data = pd.read_csv(config_dict['DATA_PATH'], index_col=list(map(int,config_dict['INDEX_COLUMNS'].split(','))))

# Create an omics measurement set object with the dataframe components
omics_mset = TypedOmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values, mapping)

# Keep the ensembl gene ID only
if 'GENE_NAME_GRAB_PATTERN' in config_dict.keys():
    ensembl_patt = re.compile(config_dict['GENE_NAME_GRAB_PATTERN'])
    omics_mset.column_names = [ensembl_patt.findall(k)[0] for k in omics_mset.column_names]

model = read_sbml_model(config_dict['MODEL_PATH']) \
    if 'MODEL_PATH' in config_dict.keys() != '' \
    else get_human1_model()

data_dicts = {'_'.join(map(str, k)): v for k, v in omics_mset.data.T.to_dict().items()}

protected = list(map(lambda x: x.strip(),config_dict['PROTECTED_REACTION_LIST'].split(','))) \
    if ('PROTECTED_REACTION_LIST' in config_dict.keys()) else []

#### parameter setup
rw = ReconstructionWrapper(model, ttg_ratio=9999)

params = {'rw': rw,
          'algorithms':['tinit','fastcore'],
          'strategies': {
              'tinit': integration_strategy_map['adjusted_score'](protected),
              'fastcore': integration_strategy_map['default_core'](0, protected)
          },
          'functions':{'minmax': MINMAX, 'minsum': MINSUM},
          'data': data_dicts
}

models_to_reconstruct = list(product(params['algorithms'], data_dicts.keys(), params['functions'].keys()))

NTHREADS = int(cpu_count() if 'NTHREADS' not in config_dict.keys() else config_dict['NTHREADS'])

def reconstruct_model(options, params):
    alg, d, func = options
    data_dict, aofunc = params['data'][d], params['functions'][func]

    oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x',
                               data=data_dict, nomenclature='custom')

    rw.run_from_omics(omics_data=oc_sample, algorithm=alg, and_or_funcs=aofunc,
                      integration_strategy=params['strategies'][alg], solver='CPLEX')

safe_threads = {'tinit': max(1, NTHREADS // 16), 'fastcore': NTHREADS}

reconstructions = {}
for k,v in safe_threads.items():
    batch = [j for j in [models_to_reconstruct] if j == k]
    reconstructions.update(dict(zip(batch,batch_run(reconstruct_model, batch, params, threads=v))))

pd.DataFrame.from_dict(reconstructions, orient='index').\
    to_csv(os.path.join(config_dict['CS_MODEL_DF_FOLDER'], config_dict['CS_MODEL_NAMES'] + '.csv'))