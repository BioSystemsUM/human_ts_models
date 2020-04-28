import requests
import json
import pandas as pd
import os
import numpy as np


def download_file(url, name):
	with open(name,'w') as f:
		f.write(requests.get(url).content.decode('utf-8'))

class CCLEDataDownloader(object):
	INTERNAL_MD_FILE_NAME = 'ccle_file_repo.csv'
	def __init__(self, local_dir=None):
		if local_dir == None:
			self.ldir = os.getcwd()
		else:
			self.ldir = local_dir
		metadata_file_name = os.path.join(self.ldir,self.INTERNAL_MD_FILE_NAME)
		if os.path.exists(metadata_file_name):
			self.file_md = pd.read_csv(metadata_file_name)
		else:
			self.file_md = self.refresh_file_md()
			self.file_md = self.file_md[~self.file_md['downloadUrl'].isnull()]
			self.file_md['filePath'] = None
			self.update_local_md_file()
			# self.file_md.to_csv(os.path.join(self.ldir, self.INTERNAL_MD_FILE_NAME))

	def update_local_md_file(self):
		os.makedirs(self.ldir)
		self.file_md.to_csv(os.path.join(self.ldir, self.INTERNAL_MD_FILE_NAME), index='file')


	def refresh_file_md(self):
		response = requests.get('https://depmap.org/portal/download/api/downloads')
		depmap_md = json.loads(response.content.decode('utf-8'))
		return pd.DataFrame(depmap_md['table'])

	def download_datasets(self, mask_function):
		subset = self.file_md[mask_function(self.file_md)]
		for k in subset.iterrows():
			file_url, file_name, release = k[1][['downloadUrl','fileName','releaseName']].tolist()

			final_folder = os.path.join(self.ldir,release)
			final_path = os.path.join(final_folder, file_name)

			if (not os.path.exists(final_path)):
				if (not os.path.exists(final_folder)):
					os.makedirs(final_folder)
				print('Downloading',file_name,'(',k[1]['size'],')','to',final_folder)
				download_file(file_url, final_path)
				print(file_name, 'downloaded and saved on', final_path)
				self.file_md.loc[file_name, 'filePath'] = final_path

			else:
				if k[1]['filePath'] != np.nan:
					self.file_md.loc[file_name,'filePath'] = final_path
					print(file_name,'already present on',final_path)
			self.update_local_md_file()

	def get_dataframes(self, mask_function):
		df_dict = {}
		self.download_datasets(mask_function)
		subset = self.file_md[mask_function(self.file_md)]
		for k in subset.iterrows():
			file_name, local_path = k[0], k[1]['filePath']
			try:
				df_dict[file_name] = pd.read_csv(local_path)
			except:
				print('Could not read',file_name)
		return df_dict

if __name__ == '__main__':

	ccle_downloader = CCLEDataDownloader(local_dir='human_ts_models/projects/breast_mcf7/data/ccle/')
	release_id = 'DepMap Public 19Q3'
	mask = lambda x: ((x['terms'] == 'ccle') & (x['release'] == release_id) & (x.index == 'CCLE_expression_full.csv'))
	data = ccle_downloader.get_dataframes(mask)
	expression_data = data['CCLE_expression_full.csv'].set_index(data['CCLE_expression_full.csv'].columns[0])