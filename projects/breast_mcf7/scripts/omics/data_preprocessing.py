from troppo.omics.id_converter import idConverter
import pandas as pd
from numpy import array, log2, inf
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

	EXPFILE_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_full.csv'
	EXPFINL_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_full_entrez.csv'
	EXPPP_PATH = 'projects/breast_mcf7/data/ccle/DepMap Public 20Q1/CCLE_expression_log2pp.csv'

	exp_df_ccle = pd.read_csv(EXPFILE_PATH, index_col=0)
	exp_df_ccle = get_mappable_df(exp_df_ccle, 'entrez_id')
	exp_df_ccle.to_csv(EXPFINL_PATH)
	exp_df_ccle_pp = preprocessing_function(exp_df_ccle)
	exp_df_ccle_pp.to_csv(EXPPP_PATH)