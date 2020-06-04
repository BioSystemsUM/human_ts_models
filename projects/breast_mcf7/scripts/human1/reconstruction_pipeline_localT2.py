import os
import pandas as pd
import re

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.utilities.parallel import batch_run
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer

from numpy import arange, log2, log, inf
from itertools import product, chain



## setting paths
ROOT_FOLDER = 'projects/breast_mcf7'
DATA_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/CCLE_expression_full.csv')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'

## create a regex pattern to identify ensembl names in a string
ensembl_patt = re.compile('ENSG[0-9]*')
## Create an identifier mapping object
## This is needed only if you need to convert genes from one nomenclature to another
mapping = IdentifierMapping('human_transcriptomics',
                  pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt',
                              index_col=0, sep = '\t'))
## Read the expression data csv
exp_data = pd.read_csv(DATA_PATH, index_col=0)

# Create an omics measurement set object with the dataframe components
omics_mset = TypedOmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values, mapping)
# Keep the ensembl gene ID only
omics_mset.column_names = [ensembl_patt.findall(k)[0] for k in omics_mset.column_names]

## Preprocessing and conversion to scores
## Define a step value to walk across
step = 0.05

## Generate different percentiles to obtain the respective values
quantiles = arange(0.1, 0.9+step, step)
gene_quantiles = omics_mset.data.quantile(quantiles)

## Generate dictionaries with the parameter options we'll need in the reconstruction part
## Minimum expression threshold
def local_threshold(x, q1df, q2df):
    return x.clip(q1df, q2df)

def process_samples(x, t, a):
    return a*log2(1+(x/t)).fillna(0).replace({inf: 10*log(2)})

lower_quantiles = omics_mset.data.quantile([0.1, 0.25, 0.5])
upper_quantiles = omics_mset.data.quantile([0.5, 0.75, 0.9])
avg_genes = omics_mset.data.mean()

lq_opts = {str(round(k, 2)):v for k,v in lower_quantiles.iterrows()}
uq_opts = {str(round(k, 2)):v for k,v in upper_quantiles.iterrows()}


## Assign a sample to reconstruct
sample_id = 'ACH-000019'
sample_df = omics_mset.data.loc[sample_id,:]
threshold_combinations = list(product(*[lq_opts.keys(),uq_opts.keys()]))

data_dicts = {(k,v): process_samples(sample_df,local_threshold(avg_genes, lq_opts[k], uq_opts[v]),5).to_dict()
              for k,v in threshold_combinations}



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

## define the biomass reaction as protected so it doesn't get removed from the core
protected = ['biomass_human']

## instantiate the reconstruction wrapped
rw = ReconstructionWrapper(model, ttg_ratio=9999)

## define a set of functions to replace AND/OR when generating scores
funcs = {'minsum':((lambda x: min([max(y, 0) for y in x])), sum),
         'minmax':((lambda x: min([max(y, 0) for y in x])), max)}

## generate a dictionary with all combinations of parameters + and/or functions
runs = dict(chain(*[[((k,n),(aof, v)) for k,v in data_dicts.items()] for n, aof in funcs.items()]))

## define a multiprocessing friendly function to reconstruct models with fastcore
def fastcore_reconstruction_func(score_tuple, params):
    aofx, data_dict = score_tuple
    oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x', data=data_dict, nomenclature='custom')
    rw = [params[k] for k in ['rw']][0]  # load parameters
    t = 5*log(2)
    try:
        def integration_fx(data_map):
            return [[k for k, v in data_map.get_scores().items() if
                     (v is not None and v > t) or k in protected]]
        return rw.run_from_omics(omics_container=oc_sample, algorithm='fastcore', and_or_funcs=aofx,
                                 integration_strategy=('custom', [integration_fx]), solver='CPLEX')
    except Exception as e:
        print(e)
        return {r: False for r in rw.model_reader.r_ids}


# to run a single model...
# define three parameters:
# min_threshold: minimum expression value required for all genes to be considered possibly active
# certain_threshold: minimum expression value required to consider a gene as certainly active
# local_quantile_threshold: a dictionary with threshold values for each gene. genes with expression between
# min_threshold and certain_threshold are considered active if their threshold is above the one specified in this
# dictionary

# to reconstruct multiple models, use this part
result_dicts = {}
labs, iters = zip(*runs.items())
output = batch_run(fastcore_reconstruction_func, iters, {'rw': rw}, threads=min(len(runs), 12))
batch_fastcore_res = dict(zip(labs, output))
result_dicts.update(batch_fastcore_res)

CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions')
pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(os.path.join(CS_MODEL_DF_FOLDER,'cs_models_lt2.csv'))