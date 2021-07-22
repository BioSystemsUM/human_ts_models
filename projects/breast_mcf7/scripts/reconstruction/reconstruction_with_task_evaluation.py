from cobamp.wrappers.external_wrappers import get_model_reader
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO
from cobra.io import read_sbml_model, load_matlab_model
from cobra.io import write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions,flux_variability_analysis, find_essential_reactions
from cobra.medium.minimal_medium import minimal_medium
from cobra.flux_analysis import pfba
from projects.breast_mcf7.src.filepaths import *
DATA_PATH = os.path.join(ROOT_FOLDER, 'data/ccle/DepMap Public 20Q1/CCLE_expression_full.csv')

from cobamp.utilities.parallel import batch_run
import pandas as pd
from json import JSONEncoder, JSONDecoder

import matplotlib.pyplot as plt
import seaborn as sns
from numpy import log, linspace, array
import re

if __name__ == '__main__':

	# helper functions
	#patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')  # find .{number} references
	patt = re.compile('_AT[0-9]{1}')  # find .{number} references
	replace_alt_transcripts = lambda x: patt.sub('', x)  # replace .{number} with nothing


	# paths to necessary files
	# this is still hardcoded

	with open(MEDIA_DICT_PATH, 'r') as f:
		media_metabolites = JSONDecoder().decode(f.read())

	# Context-specific model reconstruction #

	# model preprocessing only if the model isn't loaded already
	# this consists of:
	# - removing artificial sinks and drug modules
	# - removing blocked reactions
	if not os.path.exists(CMODL_PATH):
		model_consistent = read_sbml_model(MODEL_PATH)
		# model_consistent.remove_reactions(
		# 	[r for r in model_consistent.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
		# blocked_reactions = find_blocked_reactions(model_consistent)
		# model_consistent.remove_reactions(blocked_reactions, remove_orphans=True)
		write_sbml_model(model_consistent, CMODL_PATH)  # write a model file if it doesn't exist
	else:
		model_consistent = read_sbml_model(CMODL_PATH)

	model_consistent_media = model_consistent.copy()

	all_boundary = set([r.id for r in model_consistent_media.boundary])
	extra_metabs = {'o2[e]','ala_L[e]','gln_L[e]'}
	rem_metabs = {'alagln[e]', 'cl[e]'}

	uptake_drains = set(['EX_' + k.replace('_','__') for k in set(media_metabolites['MEM']) - rem_metabs | extra_metabs])
	uptake_drains_fix = set([k.replace('[e]','_e') for k in uptake_drains])

	for k in all_boundary - uptake_drains_fix:
		if k[:3] != 'DM_':
			model_consistent_media.reactions.get_by_id(k).bounds = (0, 1000)
		else:
			model_consistent_media.reactions.get_by_id(k).bounds = (-1000, 1000)
	for k in uptake_drains_fix:
		model_consistent_media.reactions.get_by_id(k).bounds = (-1000, 1000)

	model_consistent_media.remove_reactions(find_blocked_reactions(model_consistent_media))
	write_sbml_model(model_consistent_media, CMODL_MEDIA_PATH)


	# read the csv on DATA_PATH as a transcriptomics data set and get OmicsContainers, one per sample
	data = pd.read_csv(DATA_PATH, index_col=0)
	cols = [k.split('(')[1][:-1] for k in data.columns]
	data.columns = cols
	tab_rdr = TabularReader(path_or_df=data, nomenclature='ensembl_gene_id', omics_type='transcriptomics', cache_df=True)
	ocs = tab_rdr.to_containers()
	oc_sample = [oc for oc in ocs if oc.get_Condition() == 'ACH-000019'][0]
	oc_sample.convertIds('entrez_id')

	# create a reconstruction wrapper instance that loads the consistent model
	# GPRs are assumed to be in DNF form if the ratio between tokens and unique genes is above ttg_ratio
	# the helper function replace_alt_transcripts is used here to replace .{number} with nothing
	rw = ReconstructionWrapper(model_consistent_media, ttg_ratio=9999, gpr_gene_parse_function=replace_alt_transcripts)
	# gene score threshold to determine core reactions
	# on this dataset, genes were normalized using the following rule : 5*log(1+(expression/mean))
	# thus, if the expression is average, the final score will be 5*log(2)
	# we intend to keep reactions above this threshold as core
	ts = linspace(0.25, 0.75, 11)


	# since we will be running multiple samples, we can generalize a model reconstruction as a function
	# this function should take two arguments
	# the first will be the variable between samples (in this case, a different omics container)
	# the second should be a dictionary with static variables needed to build the model:
	# - threshold
	# - reconstruction wrapper
	def fastcore_reconstruction_func(t, params):
		rw = [params[k] for k in ['rw']][0]  # load parameters
		try:
			# if no errors appear, call the run_from_omics method passing the omics_container,
			# algorithm string, integration strategy (how the core is determined) and a solver
			# for fastcore, a threshold-based integration strategy retrieves core reactions if the score
			# is above the threshold t
			def integration_fx(data_map):
				return [[k for k, v in data_map.get_scores().items() if
				         (v is not None and v > t) or k in ['BIOMASS_maintenance'] + list(uptake_drains_fix)]]

			return rw.run_from_omics(omics_data=oc_sample, algorithm='fastcore', and_or_funcs=(min, sum),
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
			protected = list(uptake_drains_fix) + ['biomass_reaction']

			def score_apply(data_map):
				dm = {k: 0 if v is None else (min(v, 20*log(2)) - t) if k not in protected else (20 * log(2)) for k, v in data_map.items()}
				return dm

			return rw.run_from_omics(omics_data=oc_sample, algorithm='tinit', and_or_funcs=(min, sum),
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
	batch_fastcore_res = batch_run(fastcore_reconstruction_func, ts*2*log(2), {'rw': rw}, threads=min(len(ts), 12))

	#batch_tinit_res = batch_run(tinit_reconstruction_func, [5*log(2)], {'rw': rw}, threads=1)
	batch_tinit_res = [tinit_reconstruction_func(log(2), {'rw': rw})]
	# create a dict mapping tuple of algorithm and condition informations: e.g. ('fastcore','MCF7')
	# to each result from fastcore
	fastcore_res_dict = dict(zip([('fastcore', str(i)) for i in range(len(ts))], batch_fastcore_res))
	tinit_res_dict = dict(zip([('tinit', str(i)) for i in range(3)], batch_tinit_res))

	result_dicts = {}
	result_dicts.update(fastcore_res_dict)
	result_dicts.update(tinit_res_dict)

	# write these results as a dataframe for future reference
	pd.DataFrame.from_dict(result_dicts, orient='index').to_csv(CS_MODEL_DF_PATH)

	result_dicts = pd.read_csv(CS_MODEL_DF_PATH, index_col=[0, 1]).T.to_dict()
	# Task evaluation #

	# read the original model to avoid compatibility issues with the tasks (e.g. missing metabolites from the block)
	task_model = model_consistent_media.copy()
	all_metabs = {m.id for m in task_model.metabolites}
	# parse tasks from a previously existing JSON
	# the supplied file contains tasks adapted from the publication of Richelle et. al, 2019
	special_cases = {'ala_B[c]': 'ala_B_c', 'Tyr_ggn[c]':'Tyr_ggn_c'}

	if not os.path.exists(SUBTASKS_PATH):
		task_list_orig = [t for t in JSONTaskIO().read_task(TASKS_PATH)]
		task_list = []

		for task in task_list_orig:
			n1 = lambda k: k.replace('_','__').replace('[','_').replace(']','')
			n2 = lambda k: k.replace('[','_').replace(']','')

			task.inflow_dict = {n1(k) if n1(k) in all_metabs else n2(k) if n2(k) in all_metabs else special_cases[k] if k in special_cases else k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in
			                    task.inflow_dict.items()}
			task.outflow_dict = {n1(k) if n1(k) in all_metabs else n2(k) if n2(k) in all_metabs else special_cases[k] if k in special_cases else k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}
			missing_components = (set(task.inflow_dict.keys()) | set(task.outflow_dict.keys())) - all_metabs
			if len(missing_components) > 0:
				print('Task metabolites missing for', task,'-',missing_components)
			else:
				task_list.append(task)

		for task in task_list:
			task.mandatory_activity = []

		JSONTaskIO().write_task(SUBTASKS_PATH, task_list)
	else:
		task_list = JSONTaskIO().read_task(SUBTASKS_PATH)


	# tasks should be evaluated without open boundary reactions. We can easily close them on the COBRA model
	for k in task_model.boundary:
		k.knock_out()

	# get the names of all reactions in the model - this will be useful further on
	all_reactions = set([r.id for r in task_model.reactions])

	# since we have multiple models, we need to evaluate the 210 tasks for each model (54 samples)
	# first, we create a structure to hold all of these results - a dictionary
	task_eval_results = {}
	# for each k (tuple with algorithm and sample information) and result (dict with reaction presences)...
	for k, result in result_dicts.items():
		# using with statements to change the COBRA model temporarily
		# this is done to knock-out reaction not appearing the FASTCORE result
		with task_model as context_specific_model:
			protected = set([k for k, v in result.items() if v])  # get reactions included in the sample-specific model
			to_remove = all_reactions - protected  # get reactions except the protected ones
			for rid in to_remove:
				context_specific_model.reactions.get_by_id(rid).knock_out()  # knock-out reactions not in the model

			# create a task evaluator instance with the context specific model and the supplied task list and solver
			task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')

			# get task names (for future reference)
			task_names = task_eval.tasks

			# use the batch_function from the TaskEvaluator class (takes the name of a loaded task, a params
			# dictionary with the task evaluator associated to the 'tev' key) and set the amount of threads to be used
			batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=10)
		# each element in the list of results in batch_res_tasks is a tuple of length 3 with the following:
		# 0 - boolean flag representing the task evaluation
		# 1 - Solution instance used to evaluate the task
		# 2 - A dictionary with reactions supposed to be active mapped to True/False according to that criterion

		# keep only items 0 and 2 of the task result - we don't need the flux distribution
		task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
		print(k, len(protected), len([v for k, v in task_csm_res.items() if v[0]]), 'tasks completed.')
		# assign this dictionary to it's sample on the master results dictionary
		task_eval_results[k] = task_csm_res

	# save these results for later analysis as a JSON file
	with open(TASK_RESULTS_PATH, 'w') as f:
		f.write(JSONEncoder().encode([(k, v) for k, v in task_eval_results.items()]))