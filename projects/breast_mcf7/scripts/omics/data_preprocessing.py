import sys; print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

from troppo.omics.id_converter import idConverter
import pandas as pd
from numpy import array, log2, inf, apply_along_axis, power, sum, log, fromfunction, where
from functools import reduce


def dict_merge(ld):
	return reduce(lambda x,y: dict(x, **y), ld)

def list_concat(ll):
	return reduce(lambda x,y: [j if i == None else i for i,j in zip(x,y)], ll)

def get_mappable_df(exp_df_orig, final_nomenclature="entrez_id"):
	gene_ids = list(zip(*[[t.replace('(', '').replace(')', '') for t in f.split(' ')] for f in exp_df_orig.columns]))
	conversion_dict = dict_merge(
		[idConverter(gene_id_list, nomenc, final_nomenclature) for gene_id_list, nomenc in
		 zip(gene_ids, ["symbol", "ensembl_gene_id"])])

	converted_ids = [[conversion_dict[g] if g in conversion_dict.keys() else None for g in gene_id_list] for
					 gene_id_list in gene_ids]

	final_cols = array(list_concat(converted_ids))
	exp_df_mappable = exp_df_orig.loc[:, final_cols != None]
	return exp_df_mappable.rename(
		columns=dict(zip(exp_df_mappable.columns, final_cols[final_cols != None])))


def preprocessing_function(exp_df, scaling_factor=5):
	qs = [0.25, 0.75]
	quantiles = exp_df.quantile(qs, axis='index').T
	mean = exp_df.mean()
	threshold = mean.copy()
	q1r, q3r = (threshold < quantiles[qs[0]], threshold > quantiles[qs[1]])
	threshold[q1r] = quantiles.loc[q1r, qs[0]]
	threshold[q3r] = quantiles.loc[q3r, qs[1]]
	exp_df_thr = scaling_factor * log2(1 + (exp_df/threshold).fillna(0))
	return exp_df_thr

if __name__ == '__main__':
	from itertools import product

	EXPFILE_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_full.csv'
	EXPFINL_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_full_entrez.csv'
	EXPPP_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_log2pp.csv'
	EXPCOMB_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_multiple_apprx.csv'

	exp_df_ccle = pd.read_csv(EXPFILE_PATH, index_col=0)

	from numpy import arange

	qvalues = [0.1, 0.25, 0.5, 0.75, 0.9]
	quantiles = exp_df_ccle.quantile(qvalues)

	sample = 'ACH-000019'
	# global strategy
	sample_series = exp_df_ccle.loc[sample, :]

	global_thresholds = quantiles.T.apply(lambda x: x.mean())
	maxexp = log(exp_df_ccle.max().max())
	global_dicts = {('global',i,None,None):(sample_series/g).apply(log).clip(-maxexp,maxexp).to_dict()
	                for i,g in enumerate(global_thresholds)}
	# threshold for each is 0

	# local with 2 states
	# first parameter - global threshold for activity
	# second parameter - local threshold for activity
	param_combinations = list(product(*[range(len(quantiles))]*2))
	local1_dicts = {}
	for k,v in param_combinations:
		gt, lt = global_thresholds.iloc[k], quantiles.iloc[v,:]
		activity = (sample_series/gt).apply(log).clip(-maxexp, maxexp)
		gt_active = activity >= 0
		activity[gt_active] = (log(sample_series[gt_active] / lt[gt_active])*maxexp).clip(-maxexp, maxexp)
		local1_dicts[('local1',k,None,v)] = activity.to_dict()
	# threshold for activity is 0

	#local with 3 states
	local2_dicts = {}
	global_lt2_params = list(product(range(len(qvalues)),
	                                 list(zip(*where(fromfunction(lambda i,j: i < j, [len(qvalues)]*2))))))

	for v, k in global_lt2_params:
		gtl, gtu, lt = global_thresholds.iloc[k[0]], global_thresholds.iloc[k[1]], quantiles.iloc[v,:]
		upp_activity = (1+(sample_series/gtu).apply(log)).clip(-maxexp, 1+maxexp)
		gtu_inactive = upp_activity < 1
		low_activity = (sample_series[gtu_inactive]/gtl).apply(log).clip(-maxexp, maxexp)
		gtl_maybes = low_activity > 0
		activity_maybe = (sample_series[gtu_inactive][gtl_maybes]/lt[gtu_inactive][gtl_maybes]).\
			apply(log).clip(-maxexp, maxexp)
		upp_activity[activity_maybe.index] = activity_maybe.clip(-1, 1)

		local2_dicts[('local2',k[0],k[1],v)] = upp_activity.to_dict()


	pd.concat(pd.DataFrame.from_dict(d).T for d in [global_dicts, local1_dicts, local2_dicts]).to_csv(EXPCOMB_PATH)
	# exp_df_ccle = get_mappable_df(exp_df_ccle, 'ensembl_gene_id')
	# exp_df_ccle.to_csv(EXPFINL_PATH)
	# exp_df_ccle_pp = preprocessing_function(exp_df_ccle)
	# exp_df_ccle_pp.to_csv(EXPPP_PATH)
	#
	# MTB_DATA = 'projects/breast_mcf7/data/ccle/CCLE 2019/CCLE_metabolomics_20190502.csv'
	#
	# mtb_df_ccle = pd.read_csv(MTB_DATA, index_col=0).reset_index(drop=True).set_index('DepMap_ID')
	# mtb_min_max = (mtb_df_ccle - mtb_df_ccle.min()) / (mtb_df_ccle.max() - mtb_df_ccle.min())
	#
	# import matplotlib.pyplot as plt
	# exp_df_ccle.mean().hist(bins=100)
	# plt.show()
	# exp_df_ccle_pp.fillna(0).replace(inf, 0).mean().hist(bins=100)
	# plt.show()

# import scipy.stats as st
	#
	# DISTRIBUTIONS = [st.norm, st.t, st.lognorm, st.cauchy]
	#
	# fits = {}
	# for dist in DISTRIBUTIONS: fits[dist.name] = apply_along_axis(dist.fit, 0, mtb_df_ccle)
	#
	# errors = {}
	# for dist_name, params in fits.items():
	# 	dist = getattr(st, dist_name)
	# 	pdf_gen = array(list(map(lambda x: dist.pdf(x[0], *x[1][:-2], loc=x[1][-2], scale=x[1][-1]),
	# 	                   zip(mtb_df_ccle.values, fits[dist_name].T))))
	# 	errors[dist_name] = sum(power(mtb_df_ccle.values - pdf_gen, 2.0), axis=1)
	#
	# fits_df = pd.concat([pd.DataFrame(v).assign(idx=mtb_df_ccle.columns, dist=k) for k,v in fits.items()]).set_index(['idx','dist'])
	# sse_df = pd.DataFrame(errors, index=mtb_df_ccle.columns)
