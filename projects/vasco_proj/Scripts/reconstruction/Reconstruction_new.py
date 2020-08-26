from troppo.omics.readers.generic import TabularReader
from cobamp.wrappers.external_wrappers import get_model_reader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO
from cobra.io import read_sbml_model, load_matlab_model  # ler o modelo com  o cobra
from cobra.io import write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions, flux_variability_analysis, find_essential_reactions
from cobra.medium.minimal_medium import minimal_medium
from cobra.flux_analysis import pfba
from projects.vasco_proj.src.Paths import *

from cobamp.utilities.parallel import batch_run
import pandas as pd
from json import JSONEncoder, JSONDecoder

import matplotlib.pyplot as plt
import seaborn as sns
from numpy import log, linspace, array
import re


if __name__ == '__main__':
    # helper functions
    # patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')  # find .{number} references
    patt = re.compile('_AT[0-9]{1}')  # find .{number} references
    replace_alt_transcripts = lambda x: patt.sub('', x)  # replace .{number} with nothing

    # paths to necessary files
    # this is still hardcoded

    # Context-specific model reconstruction #

    # model preprocessing only if the model isn't loaded already
    # this consists of:
    # - removing artificial sinks and drug modules
    # - removing blocked reactions

    if not os.path.exists(CMODL_PATH):
        model_consistent = read_sbml_model(MODEL_PATH)
        model_consistent.remove_reactions([r for r in model_consistent.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
        blocked_reactions = find_blocked_reactions(model_consistent)
        model_consistent.remove_reactions(blocked_reactions, remove_orphans=True)
        write_sbml_model(model_consistent, CMODL_PATH)  # write a model file if it doesn't exist
    else:
        model_consistent = read_sbml_model(CMODL_PATH)

    # files=data=pd.read_excel(TRANSCRIPTOMICS_PATH + ("Proc/Compiled.xlsx"), header=0, index_col=[0,1])
    files=os.listdir(TRANSCRIPTOMICS_PATH + "Proc")
    f = pd.read_excel(TRANSCRIPTOMICS_PATH + ("Proc/Comp_z.xlsx"), header=0, index_col=[0,1])

    # read the csv on DATA_PATH as a transcriptomics data set and get OmicsContainers, one per sample
    # ocs = TabularReader(path_or_df=file, nomenclature='entrez_id', omics_type='transcriptomics').to_containers()


    ocs = TabularReader(path_or_df=f, nomenclature='entrez_id', omics_type='transcriptomics').to_containers()
    oc_sample = [oc for oc in ocs if oc.get_Condition() == ('Supdata_proc_z_score.xlsx',1)][0]

    # oc_sample = [oc for oc in ocs if oc.get_Condition() == (0)][0]


    # create a reconstruction wrapper instance that loads the consistent model
    # GPRs are assumed to be in DNF form if the ratio between tokens and unique genes is above ttg_ratio
    # the helper function replace_alt_transcripts is used here to replace .{number} with nothing

    rw = ReconstructionWrapper(model_consistent, ttg_ratio=9999, gpr_gene_parse_function=replace_alt_transcripts)


    # gene score threshold to determine core reactions
    # on this dataset, genes were normalized using the following rule : 5*log(1+(expression/mean))
    # thus, if the expression is average, the final score will be 5*log(2)
    # we intend to keep reactions above this threshold as core
    #sample1
    # ts = linspace(0.25, 1, 4)
    # sample2
    ts = linspace(-1.2, 0, 5)
    # sample3
    # ts = linspace(-1.2, 0.5, 6)

    # since we will be running multiple samples, we can generalize a model reconstruction as a function
    # this function should take two arguments
    # the first will be the variable between samples (in this case, a different omics container)
    # the second should be a dictionary with static variables needed to build the model:
    # - threshold
    # - reconstruction wrapper
    xxx = None
    
    def fastcore_reconstruction_func(t, params):
        rw = [params[k] for k in ['rw']][0]  # load parameters
        try:
            # if no errors appear, call the run_from_omics method passing the omics_container,
            # algorithm string, integration strategy (how the core is determined) and a solver
            # for fastcore, a threshold-based integration strategy retrieves core reactions if the score
            # is above the threshold t
            def integration_fx(data_map):
                return [[k for k, v in data_map.get_scores().items() if
                         (v is not None and v > t) or k in ['BIOMASS_maintenance']]]

            return rw.run_from_omics(omics_container=oc_sample, algorithm='fastcore', and_or_funcs=(min, sum),
                                     integration_strategy=('custom', [integration_fx]), solver='CPLEX')
        except Exception as e:
            # the result from run_from_omics is a dict mapping reaction ids and a boolean flag - True if
            # the reaction is in the model or false otherwise
            # in case an error arises, assume all reactions are False
            print(e)
            return {r: False for r in rw.model_reader.r_ids}


    def tinit_reconstruction_func(t, params):
        rw = [params[k] for k in ['rw']][0]  # load parameters
        try:
            # if no errors appear, call the run_from_omics method passing the omics_container,
            # algorithm string, integration strategy (how the core is determined) and a solver
            # for fastcore, a threshold-based integration strategy retrieves core reactions if the score
            # is above the threshold t
            protected = ['biomass_reaction']

            def score_apply(data_map):
                dm = {k: 0 if v is None else (min(v, 10 * log(2)) - t) if k not in protected else (10 * log(2)) for k, v
                      in data_map.items()}
                return dm

            return rw.run_from_omics(omics_container=oc_sample, algorithm='tinit', and_or_funcs=(min, sum),
                                     integration_strategy=('continuous', score_apply), solver='CPLEX')
        except Exception as e:
            # the result from run_from_omics is a dict mapping reaction ids and a boolean flag - True if
            # the reaction is in the model or false otherwise
            # in case an error arises, assume all reactions are False
            print(e)
            # raise e
            return {r: False for r in rw.model_reader.r_ids}


    # parallel reconstruction can be achieved with the batch_run function that takes in 4 key arguments:
    # - the function to apply multiple times - reconstruction_func
    # - a list of objects for which we want to apply the function - ocs (list of omics_containers)
    # - a dictionary with the static params - containing a 't' and 'rw' entry (see reconstruction_func)
    # - an integer value specifying the amount of parallel processes to be run


    # fastcore_res=fastcore_reconstruction_func(ts * 10 * log(2), {'rw': rw})
    # batch_fastcore_res = batch_run(fastcore_reconstruction_func, ts * 10 * log(2), {'rw': rw}, threads=min(len(ts), 12))

    # tinit_res=tinit_reconstruction_func([5*log(2)], {'rw': rw})

    # batch_tinit_res = batch_run(tinit_reconstruction_func, [5*log(2)], {'rw': rw}, threads=1)
    # batch_tinit_res = [tinit_reconstruction_func(5 * log(2), {'rw': rw})]
    # create a dict mapping tuple of algorithm and condition informations: e.g. ('fastcore','MCF7')
    # to each result from fastcore
    # fastcore_res_dict = dict(zip([('fastcore', str(i)) for i in range(len(ts))], fastcore_res))
    # tinit_res_dict = dict(zip([('tinit', str(i)) for i in range(3)], tinit_res))

    # result_dicts = {}
    # result_dicts.update(fastcore_res_dict)
    # result_dicts.update(tinit_res_dict)
    fastcore_res_dict = {}

    for t in ts:
        fastcore_res = fastcore_reconstruction_func(t * 10 * log(2), {'rw': rw})
        fastcore_res_dict["fastcore" + str(t)] = fastcore_res

    pd.DataFrame.from_dict(fastcore_res_dict,orient="index").to_excel('C:/Users/vasco/PycharmProjects/human_ts_models/projects/vasco_proj/Data/Transcriptomics/Proc/r3d_compact_mcf7_fastcore_sample_1.xlsx')
    # write these results as a dataframe for future reference
    # pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(CS_MODEL_DF_PATH)
    #
    # result_dicts = pd.read_csv(CS_MODEL_DF_PATH, index_col=[0, 1]).T.to_dict()

