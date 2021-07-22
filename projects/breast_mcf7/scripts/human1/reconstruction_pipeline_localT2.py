import os
import pandas as pd
import re

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions, gapfill, pfba
from cobra.io import read_sbml_model, write_sbml_model

from cobamp.algorithms.kshortest import KShortestProperties, K_SHORTEST_OPROPERTY_BIG_M_CONSTRAINTS
from cobamp.utilities.file_io import pickle_object
from cobamp.utilities.parallel import batch_run
from cobamp.wrappers.external_wrappers import get_model_reader
from projects.breast_mcf7.scripts.omics import growth_media
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import TypedOmicsMeasurementSet, IdentifierMapping, OmicsContainer

from numpy import arange, log2, log, inf
from itertools import product, chain

## setting paths
ROOT_FOLDER = 'projects/breast_mcf7'
DATA_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/CCLE_expression_full.csv')
PROT_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/Proteomics/protein_quant_current_normalized.csv.gz')
SAMPLE_INFO = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/sample_info_v2.csv')
MODEL_PATH = 'projects/breast_mcf7/support/Human-GEM.xml'

## create a regex pattern to identify ensembl names in a string
ensembl_patt = re.compile('ENSG[0-9]*')
## Create an identifier mapping object
## This is needed only if you need to convert genes from one nomenclature to another
mapping = IdentifierMapping('human_transcriptomics',
                  pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt',
                              index_col=0, sep = '\t'))
sample_info = pd.read_csv(SAMPLE_INFO, index_col=0)
depmap_to_ccle = {v:k for k,v in sample_info['CCLE_Name'].to_dict().items()}


exp_data = pd.read_csv(DATA_PATH, index_col=0)

# Create an omics measurement set object with the dataframe components
omics_mset = TypedOmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values, mapping)

# Keep the ensembl gene ID only
omics_mset.column_names = [ensembl_patt.findall(k)[0] for k in omics_mset.column_names]

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


#
# lqp_opts = {str(round(k, 2)):v for k,v in prot_mset.data.quantile([0.1, 0.25, 0.5]).iterrows()}
# uqp_opts = {str(round(k, 2)):v for k,v in prot_mset.data.quantile([0.5, 0.75, 0.9]).iterrows()}
# avg_prots = prot_mset.data.mean()
#
# ## Assign a sample to reconstruct
sample_id = 'ACH-000019'
sample_df = omics_mset.data.loc[sample_id,:]
# prot_sample_df = prot_mset.data.loc[sample_id,:]

threshold_combinations = list(product(*[lq_opts.keys(),uq_opts.keys()]))

data_dict_rna = {(k,v): process_samples(sample_df,local_threshold(avg_genes, lq_opts[k], uq_opts[v]),5).to_dict()
              for k,v in threshold_combinations}
#
# data_dict_prot = {(k,v): process_samples(prot_sample_df,local_threshold(avg_prots, lqp_opts[k], uqp_opts[v]),5).to_dict()
#               for k,v in threshold_combinations}

data_dicts = {k:{i:j for i,j in data_dict_rna[k].items()} for k in threshold_combinations}

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
    t = 5*log(2)
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
    def_val = max(oc_sample.get_Data().values())
    t = 5*log(2)
    try:
        def tinit_integration_fx(data_map):
            vtransform = lambda x,t: x/t if x >= t else (x-t)/t
            return {k:(def_val if k in protected else vtransform(v,t) if v is not None else 0)
                    for k, v in data_map.get_scores().items()}
        return rw.run_from_omics(omics_data=oc_sample, algorithm='tinit', and_or_funcs=aofx,
                                 integration_strategy=('custom', [tinit_integration_fx]), solver='CPLEX')
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

output_tinit = batch_run(tinit_reconstruction_func, iters, {'rw': rw}, threads=1)
result_dicts.update(dict(zip(labs, output)))
result_dicts.update(batch_fastcore_res)

CS_MODEL_DF_FOLDER = os.path.join(ROOT_FOLDER, 'results/human1/reconstructions')
pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(os.path.join(CS_MODEL_DF_FOLDER,'cs_models_lt2_prot.csv'))

model_df = pd.read_csv(os.path.join(CS_MODEL_DF_FOLDER,'cs_models_lt2.csv'), index_col=[0,1])


media_metabolites = growth_media.get_medium_metabolites('MEM', 'iHuman_id')
# extra_metabolites = {'m02631s', 'm02039s', 'm02394s', 'm01822s',
#                      'm02982s', 'm01938s', 'm01401s', 'm01330s',
#                      'm02394s', 'm01361s', 'm01822s', 'm01361s',
#                      'm01570s', 'm03147s', 'm01385s', 'm01629s'}
#
# {k:model.metabolites.get_by_id(k).name for k in extra_metabolites}
extra_metabs = {'m10005s'}
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
ksh_cfg[K_SHORTEST_OPROPERTY_BIG_M_VALUE] = 1e3
ksh_cfg[K_SHORTEST_OPROPERTY_FORCE_NON_CANCELLATION] = True
ksh_cfg[K_SHORTEST_OPROPERTY_MAXSOLUTIONS] = 1

from troppo.methods_wrappers import GapfillWrapper

## TODO: remove drains for compounds without formula (to eliminate metabolite pools, etc...)
efm_gapfill_result = {}
for i, mrow in model_df.iterrows():
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

    efm_gapfill_result[i] = gapfill_res[0]
    print(i, gapfill_res[0])


pfba_gapfill = {}
for i, mrow in model_df.iterrows():
    mdict = mrow.to_dict()
    with model as m:
        # m.remove_reactions([k for k,v in mdict.items() if not v and (k in cobamp_model.reaction_names) and (k not in outflows | set(media_boundrx.values()))])
        # for r in outflows:
        #     m.reactions.get_by_id(r).bounds = (-1000, 1000)
        # for r,v in media_boundrx.items():
        #     try:
        #         m.reactions.get_by_id(v).bounds = (-1000, 1000)
        #     except:
        #         print(v,'not found')
        for k,v in mdict.items():
            if k not in cobamp_model_boundrx_rev.keys() and not v:
                m.reactions.get_by_id(k).bounds = (0, 0)
        for k in outflows:
            m.reactions.get_by_id(k).bounds = (-0.001, 1000)

        rvars = [m.reactions.get_by_id(r).reverse_variable for r in outflows]
        m.objective = sum(rvars)
        m.objective_direction = 'min'
        m.reactions.biomass_human.bounds = (0.01, 1000)
        try:
            sol = m.optimize()

            m.objective = 'biomass_human'
            m.objective_direction = 'max'

            ofl_flx = sol.fluxes[outflows]
            for k in outflows:
                m.reactions.get_by_id(k).bounds = (0, 1000)
            print('\tadded',ofl_flx[ofl_flx < 0].index)
            for k in ofl_flx[ofl_flx < 0].index:
                m.reactions.get_by_id(k).bounds = (-0.001, 1000)
            m.reactions.biomass_human.bounds = (0, 1000)
            pfba_sol = pfba(m)


            pfba_gapfill[i] = {'fixed_gaps': ofl_flx[ofl_flx < 0].index,
                               'bounds': {m:m.bounds for m in m.reactions}}
        except:
            pass

pickle_object(pfba_gapfill, os.path.join(CS_MODEL_DF_FOLDER, 'cs_models_lt2_prot_fba_gapfill.csv'))
print(model.summary(pfba_sol, names=True, threshold=0.000001))

