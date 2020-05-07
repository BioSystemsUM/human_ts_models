from json import JSONDecoder
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy as np
from projects.breast_mcf7.src.filepaths import *
from troppo.tasks.task_io import JSONTaskIO
import matplotlib.pyplot as plt

if __name__ == '__main__':
	import os
	## paths to necessary files
	## this is still hardcoded

	task_list = [t for t in JSONTaskIO().read_task(TASKS_PATH)]

	with open(TASK_RESULTS_PATH, 'r') as f:
		task_eval_results_slim = JSONDecoder().decode(f.read())

	fastcore_res_df = pd.read_csv(CS_MODEL_DF_PATH, index_col=[0, 1])
	fastcore_res_dict = fastcore_res_df.T.to_dict()

	k, v = zip(*task_eval_results_slim)
	task_eval_results_slim = dict(zip([tuple(i) for i in k], v))

	task_full_name_mapper = {t.name: str(t.should_fail)[0] + '-' + t.name + ' ' + t.annotations['description'] for t
							 in task_list}
	task_subsystem_mapper = {t.name: t.annotations['subsystem'] for t in task_list}

	task_res_df_orig = pd.DataFrame.from_dict(
		{k: {tk: tv[0] for tk, tv in v.items()} for k, v in task_eval_results_slim.items()}, orient='index').rename(
		columns=task_full_name_mapper)
	task_df_subsystem_series = pd.Series({fn:task_subsystem_mapper[t] for t,fn in task_full_name_mapper.items()})

	task_res_df = task_res_df_orig.copy()

	df_full_plot = task_res_df_orig.T
	df_full_plot['subsystem'] = task_df_subsystem_series
	df_full_plot = df_full_plot.groupby('subsystem').mean()




	def heatmap(df_from_dict, plot_path, figsize=(40, 20)):
		# ord_samp = sch.fcluster(samp_linkg, 0.7*max(samp_linkg[:,2]), 'distance')
		# ord_task = sch.fcluster(task_linkg, 0.7*max(task_linkg[:,2]), 'distance')
		df_from_dict_plot = df_from_dict
		plt.figure(figsize=figsize)
		plt.subplots_adjust(hspace=0.001, wspace=0.001)
		plt.gca().set_aspect(1)
		plt.pcolormesh(df_from_dict_plot)
		plt.yticks(np.arange(0.5, len(df_from_dict_plot.index), 1), [t for t in df_from_dict_plot.index])
		plt.xticks(np.arange(0.5, len(df_from_dict_plot.columns), 1), [k for k in df_from_dict_plot.columns],
				   rotation=90)
		plt.savefig(plot_path)


	PLOT_FOLDER = os.path.join(ROOT_FOLDER, 'results/plots')
	if not os.path.exists(PLOT_FOLDER): os.makedirs(PLOT_FOLDER)

	heatmap(task_res_df_orig, os.path.join(PLOT_FOLDER, 'CCLE_mcf7_fastcore_orig_minsum.png'), figsize=(30, 20))

	heatmap(df_full_plot, os.path.join(PLOT_FOLDER, 'CCLE_mcf7_fastcore_full_subsys_minsum.png'),
	        figsize=(15, 15))

