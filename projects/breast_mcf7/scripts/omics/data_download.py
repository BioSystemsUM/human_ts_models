from shared.src.data.ccle_tools import CCLEDataDownloader

ccle_downloader = CCLEDataDownloader(local_dir='projects/breast_mcf7/data/ccle/')

def latest_expression(md):
	return (md['fileName'] == 'CCLE_expression_full.csv') & (md['releaseName'] == 'DepMap Public 20Q1')

def latest_info(md):
	return (md['fileName'] == 'sample_info_v2.csv') & (md['releaseName'] == 'DepMap Public 20Q1')

# def latest_metabolomics(md):
# 	return md['fileName'] == 'CCLE_metabolomics_20190502.csv'

def latest_gene_dep(md):
	return (md['fileName'] == 'Achilles_gene_effect.csv') & (md['releaseName'] == 'DepMap Public 20Q1')


def latest_proteomics(md):
	return (md['fileName'] == 'protein_quant_current_normalized.csv')

# ccle_downloader.download_datasets(latest_metabolomics)
ccle_downloader.download_datasets(latest_expression)
ccle_downloader.download_datasets(latest_info)
ccle_downloader.download_datasets(latest_gene_dep)
ccle_downloader.download_datasets(latest_proteomics)

gtex_median_url = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8' \
'_RNASeQCv1.1.9_gene_median_tpm.gct.gz'