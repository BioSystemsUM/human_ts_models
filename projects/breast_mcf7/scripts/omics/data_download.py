from shared.src.data.ccle_tools import CCLEDataDownloader

ccle_downloader = CCLEDataDownloader(local_dir='human_ts_models/projects/breast_mcf7/data/ccle/')

def latest_expression(md):
	return (md['fileName'] == 'CCLE_expression_full.csv') & (md['releaseName'] == 'DepMap Public 20Q1')

def latest_info(md):
	return (md['fileName'] == 'sample_info_v2.csv') & (md['releaseName'] == 'DepMap Public 20Q1')

ccle_downloader.download_datasets(latest_expression)
ccle_downloader.download_datasets(latest_info)

gtex_median_url = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8' \
'_RNASeQCv1.1.9_gene_median_tpm.gct.gz'

