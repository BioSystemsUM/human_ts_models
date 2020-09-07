import os
import pandas as pd
import re

import sys


print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

from cobamp.utilities.parallel import batch_run
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer
from projects.breast_mcf7.scripts.human1.models import get_human1_model
from itertools import chain

ROOT_FOLDER = 'projects/breast_mcf7'
DATA_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/CCLE_expression_ACH-000019_scores.csv')
SAMPLE_INFO = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/sample_info_v2.csv')
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions/mcf7_comparison')

if not os.path.exists(CS_MODEL_DF_FOLDER): os.makedirs(CS_MODEL_DF_FOLDER)
CS_MODEL_NAMES = 'cs_models_all_combinations'
NTHREADS = 12

## create a regex pattern to identify ensembl names in a string
ensembl_patt = re.compile('ENSG[0-9]*')
## Create an identifier mapping object
## This is needed only if you need to convert genes from one nomenclature to another
mapping = IdentifierMapping('human_transcriptomics',
                  pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt',
                              index_col=0, sep = '\t'))
sample_info = pd.read_csv(SAMPLE_INFO, index_col=0)
depmap_to_ccle = {v:k for k,v in sample_info['CCLE_Name'].to_dict().items()}


exp_data = pd.read_csv(DATA_PATH, index_col=list(range(4)))

# Create an omics measurement set object with the dataframe components
omics_mset = TypedOmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values, mapping)

# Keep the ensembl gene ID only
omics_mset.column_names = [ensembl_patt.findall(k)[0] for k in omics_mset.column_names]


## read the metabolic model and simplify it
model = get_human1_model()

data_dicts = {'_'.join(map(str,k)):v for k,v in omics_mset.data.T.to_dict().items()}
## define the biomass reaction as protected so it doesn't get removed from the core
protected = ['biomass_human','HMR_10023','HMR_10024']

## define a set of functions to replace AND/OR when generating scores
funcs = {'minsum':(min, sum),
         'minmax':(min, max)}

## generate a dictionary with all combinations of parameters + and/or functions
runs = dict(chain(*[[((k,n),(aof, v)) for k,v in data_dicts.items()] for n, aof in funcs.items()]))

## instantiate the reconstruction wrapper

## when reconstructing models, you should use it unconstrained with open exchanges
## if your core set contains reactions that are incompatible with the media, the algorithm
## won't be able to reconstruct the model. for this reason, one should first reconstruct and then gapfill
rw = ReconstructionWrapper(model, ttg_ratio=9999)


## define a multiprocessing friendly function to reconstruct models with fastcore
def fastcore_reconstruction_func(score_tuple, params):
    aofx, data_dict = score_tuple
    oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x', data=data_dict, nomenclature='custom')
    rw = [params[k] for k in ['rw']][0]  # load parameters
    t = 0
    def integration_fx(data_map):
        return [[k for k, v in data_map.get_scores().items() if
                 (v is not None and v > t) or k in protected]]

    try:
        return rw.run_from_omics(omics_data=oc_sample, algorithm='fastcore', and_or_funcs=aofx,
                                 integration_strategy=('custom', [integration_fx]), solver='CPLEX')
    except Exception as e:
        print(e)
        return {r: False for r in rw.model_reader.r_ids}

def tinit_reconstruction_func(score_tuple, params):
    aofx, data_dict = score_tuple
    oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x', data=data_dict, nomenclature='custom')
    rw = [params[k] for k in ['rw']][0]  # load parameters
    try:
        def tinit_integration_fx(data_map):
            maxv = max([k for k in data_map.get_scores().values() if k is not None])
            scores = {k:(v/maxv if v < 0 else v) if v is not None else 0 for k, v in data_map.get_scores().items()}
            scores.update({x:max(scores.values()) for x in protected})
            return scores

        return rw.run_from_omics(omics_data=oc_sample, algorithm='tinit', and_or_funcs=aofx,
                                 integration_strategy=('custom', [tinit_integration_fx]), solver='CPLEX')
    except Exception as e:
        print(e)
        return {r: False for r in rw.model_reader.r_ids}



result_dicts = {}
labs, iters = zip(*runs.items())

tlabs = [tuple(['tinit']+list(l)) for l in labs]
output = batch_run(tinit_reconstruction_func, iters, {'rw': rw}, threads=min(len(runs), 3))
batch_tinit_res = dict(zip(tlabs, output))
result_dicts.update(batch_tinit_res)

flabs = [tuple(['fastcore']+list(l)) for l in labs]
output = batch_run(fastcore_reconstruction_func, iters, {'rw': rw}, threads=min(len(runs), NTHREADS))
batch_fastcore_res = dict(zip(flabs, output))
result_dicts.update(batch_fastcore_res)


pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(os.path.join(CS_MODEL_DF_FOLDER,CS_MODEL_NAMES+'.csv'))