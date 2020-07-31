import os
import pandas as pd
import re

import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions, pfba
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.utilities.file_io import pickle_object
from cobamp.utilities.parallel import batch_run
from cobamp.wrappers.external_wrappers import get_model_reader
from projects.breast_mcf7.scripts.omics import growth_media
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer

from numpy import log
from itertools import chain

ROOT_FOLDER = 'projects/breast_mcf7'
DATA_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/CCLE_expression_multiple_apprx.csv')
SAMPLE_INFO = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/sample_info_v2.csv')
CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions')

MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'
NTHREADS = 40

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
if not os.path.exists(MODEL_PATH):
    path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/modelFiles/xml/HumanGEM.xml')
    model = read_sbml_model(path)
    write_sbml_model(model, MODEL_PATH)
else:
    model = read_sbml_model(MODEL_PATH)
model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
blocked = find_blocked_reactions(model)
model.remove_reactions(blocked)


data_dicts = {'_'.join(map(str,k)):v for k,v in omics_mset.data.T.to_dict().items()}
## define the biomass reaction as protected so it doesn't get removed from the core
protected = ['biomass_human']

## define a set of functions to replace AND/OR when generating scores
funcs = {'minsum':((lambda x: min([max(y, 0) for y in x])), sum),
         'minmax':((lambda x: min([max(y, 0) for y in x])), max)}

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
            return {k:v if v is not None else 0 for k, v in data_map.get_scores().items()}
        return rw.run_from_omics(omics_data=oc_sample, algorithm='tinit', and_or_funcs=aofx,
                                 integration_strategy=('custom', [tinit_integration_fx]), solver='CPLEX')
    except Exception as e:
        print(e)
        return {r: False for r in rw.model_reader.r_ids}


result_dicts = {}
labs, iters = zip(*runs.items())
flabs = [tuple(['fastcore']+list(l)) for l in labs]
output = batch_run(fastcore_reconstruction_func, iters, {'rw': rw}, threads=min(len(runs), NTHREADS))
batch_fastcore_res = dict(zip(flabs, output))
result_dicts.update(batch_fastcore_res)

tlabs = [tuple(['tinit']+list(l)) for l in labs]
output = batch_run(tinit_reconstruction_func, iters, {'rw': rw}, threads=min(len(runs), 1))
batch_tinit_res = dict(zip(tlabs, output))
result_dicts.update(batch_tinit_res)

pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(os.path.join(CS_MODEL_DF_FOLDER,'cs_models_all_combinations.csv'))