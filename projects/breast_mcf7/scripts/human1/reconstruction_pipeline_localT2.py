import os
import pandas as pd
import re

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.utilities.parallel import batch_run
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer
import matplotlib.pyplot as plt
from numpy import arange
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
min_exp_options = {str(round(quantiles[i], 2)): gene_quantiles.loc[quantiles[i],:].mean() for i in range(4)}

## Certain expression threshold
cert_exp_options = {str(round(quantiles[i], 2)): gene_quantiles.loc[quantiles[i],:].mean()
                    for i in range(len(quantiles)//2, len(quantiles))}

## Gene-wise quantiles
local_quantile_options = {str(round(quantiles[i], 2)):gene_quantiles.loc[quantiles[i],:]
                          for i in range((len(quantiles)//2)+1)}

## Assign a sample to reconstruct
sample_id = 'ACH-000019'
sample_df = exp_data.loc[[sample_id],:]

## Define a function to take a dictionary with genes mapped to their expression value and a set of parameters,
## returning another dictionary with processed values according to the following rules:
## x >= certexp? convert to x/certexp : same as fold change
## minexp <= x < certexp? convert to (x-local_quantile)/q, but bound the value to stay between -1 and 1
## x < minexp? convert to ((x-minexp)/minexp) -1

def process_values(exp_dict, minexp, certexp, local_quantiles):
    def funcv(x, q):
        return x/certexp if x >= certexp else min(max((x-q)/q, -1),1) if x >= minexp else -1-((x-minexp)/minexp)
    return {k:funcv(v, local_quantiles[k]) for k,v in exp_dict.items()}


# generate all possible parameter combinations
threshold_combinations = list(product(*[min_exp_options.keys(),cert_exp_options.keys(),local_quantile_options.keys()]))
# convert the data from pandas series to dict
sample_dict = omics_mset.data.loc[sample_id,:].to_dict()

# generate the processed expression scores according to each parameter value
data_dicts = {(a,b,c):process_values(sample_dict,
                             min_exp_options[a], cert_exp_options[b], local_quantile_options[c])
              for a,b,c in threshold_combinations}


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
    try:
        def integration_fx(data_map):
            return [[k for k, v in data_map.get_scores().items() if
                     (v is not None and v > 0) or k in protected]]
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

# min_threshold = min_exp_options['0.1'] # average value of the 10th percentile for all genes
# certain_threshold = cert_exp_options['0.5'] # average value of the 50th percentile for all genes
# local_quantile_threshold = local_quantile_options['0.25'] # 25th percentile values for each gene
#
# sample_data_dict = process_values(sample_dict, min_threshold, certain_threshold, local_quantile_threshold)
# and_or_functions = funcs['minsum']
#
# cs_model = fastcore_reconstruction_func((and_or_functions, sample_data_dict), {'rw': rw})


# to reconstruct multiple models, use this part snippet
result_dicts = {}
labs, iters = zip(*runs.items())
output = batch_run(fastcore_reconstruction_func, iters, {'rw': rw}, threads=min(len(runs), 12))
batch_fastcore_res = dict(zip(labs, output))
result_dicts.update(batch_fastcore_res)

CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions')
pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(os.path.join(CS_MODEL_DF_FOLDER,'cs_models.csv'))
model_df = pd.DataFrame.from_dict(result_dicts, orient='index')


## save result_dicts somewhere