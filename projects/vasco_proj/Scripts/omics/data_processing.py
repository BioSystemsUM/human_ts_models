# data processing script
# a set of functions to help the reading and pre-processing of the dataset
import os
import mygene
import pandas as pd
from projects.vasco_proj.src.Paths import *


# types: "Proteomics","Metabolomics","Transcriptomics"

def get_files(type):
    if type == 'Metabolomics':
        files = [i for i in os.listdir(METABOLOMICS_PATH) if i.endswith(".xlsx")]
        path = METABOLOMICS_PATH
    if type == 'Proteomics':
        files = [i for i in os.listdir(PROTEOMICS_PATH) if i.endswith(".xlsx")]
        path = PROTEOMICS_PATH
    if type == 'Transcriptomics':
        files = [i for i in os.listdir(TRANSCRIPTOMICS_PATH) if i.endswith(".xlsx")]
        path = TRANSCRIPTOMICS_PATH
    return (files, path)



def mapping_convertion(exp_file_original, final_nomenclature="entrezgene", species="human", sheet_name="Sheet1"):
    file = pd.read_excel(exp_file_original)
    col_convert = file.iloc[:,0]
    xli = []
    for i in col_convert:
        xli.append(str(i))
    mg = mygene.MyGeneInfo()
    converted_list = mg.querymany(xli, scopes='symbol', fields=final_nomenclature,
                                  species=species)  # objeto do tipo mygene com uma lista de dicionários

    final_dict_conversion = {k['query']: k['entrezgene'] for k in converted_list if (k['query'] in xli and 'entrezgene' in k.keys())}
    return final_dict_conversion

def map(exp_file,dict):
    for i in exp_file.columns:
        if i not in dict.keys():
            exp_file= exp_file.drop(str(i),1)
    data_f = exp_file.rename(columns=dict)
    return data_f

# zero for rows - axis
# one for columns - axis

# removes specific columns or rows

def data_select(file, gene_names_index, data_index=(), transpose=False):
    gene_names = file.iloc[:, gene_names_index]
    data = file.iloc[:, data_index[0]:data_index[1]]
    data_file = pd.concat([gene_names, data], 1)
    if transpose:
        data_file.T
    return file


def data_cut(exp_file, numbers=[], axis=1):
        if axis == 1:
            exo_file.drop(exp_file.columns[numbers], axis=axis, inplace=True)
        else:
            exp_file.drop(exp_file.row[numbers], axis=axis, inplace=True)
        return exp_file


# zero for rows - axis
# one for columns - axis
# two for both

def spliter(exp_file, axis=None, row=(), column=()):
    if axis == 0:
        file1 = exp_file[row[0]:row[1], :]
        file2 = exp_file[row[1]:, :]
    if axis == 1:
        file1 = exp_file[:, column[0]:column[1]]
        file2 = exp_file[:, column[1]:]
    return file1, file2


# a: an array like object containing data
# axis: the axis along which to calculate the z-scores. Default is 0.
# ddof: degrees of freedom correction in the calculation of the standard deviation. Default is 0.
# nan_policy: how to handle when input contains nan. Default is propagate, which returns nan. ‘raise’ throws an error and ‘omit’ performs calculations ignoring nan values.

# DataFrame.apply(func, axis=0, raw=False, result_type=None, args=(), **kwds
# DataFrame.dropna().apply()

# methods for the pre_processing funtion
# "maximum_absolute"
# "min_max"
# "z_score"
# "robust_scaling"
methods=['maximum_absolute','min_max','z_score','robust_scaling']

def pre_processing(dataframe, method):
    if method == "maximum_absolute":
        df_scaled = dataframe.apply(maximum_absolute_scaling, axis=1)
        return (df_scaled,method)
    if method == "min_max":
        df_scaled = dataframe.apply(min_max_scaling, axis=1)
        return (df_scaled, method)
    if method == "z_score":
        df_scaled = dataframe.apply(z_score, axis=1)
        return (df_scaled, method)
    if method == "robust_scaling":
        df_scaled = dataframe.apply(robust_scaling, axis=1)
        return (df_scaled, method)
    else:
        print("unknown method")


# The maximum absolute scaling
# The maximum absolute scaling rescales each feature between -1 and 1 by dividing every observation by its maximum absolute value.
# apply the maximum absolute scaling in Pandas using the .abs() and .max() methods

# def maximum_absolute_scaling(dataframe):
#     for column in dataframe.columns:
#         dataframe[column] = dataframe[column] / dataframe[column].abs().max()
#     return dataframe

def maximum_absolute_scaling(column):
    column = column / column.abs().max()  # ver dataframe.apply(lambda x: x.abs().max(), )
    return column


# The min-max approach (often called normalization) rescales the feature to a fixed range of [0,1] by
# subtracting the minimum value of the feature and then dividing by the range.

# def min_max_scaling(df):
#     for column in df_norm.columns:
#         df_norm[column] = (df_norm[column] - df_norm[column].min()) / (df_norm[column].max() - df_norm[column].min())
#     return column

def min_max_scaling(column):
        column = (column - column.min()) / (column.max() - column.min())
    return column

# The z-score method (often called standardization) transforms the data into a distribution with a mean of 0 and a standard deviation of 1.
# Each standardized value is computed by subtracting the mean of the corresponding feature and then dividing by the standard deviation.

# def z_score(df):
#     for column in df_std.columns:
#         df_std[column] = (df_std[column] - df_std[column].mean()) / df_std[column].std()
#     return df_std

def z_score(column):
    column = (column - column.mean()) / column.std()
    return column

# This method comes in handy when working with data sets that contain many outliers because it uses statistics
# that are robust to outliers (median and interquartile range), in contrast with the previous scalers, which use statistics that are highly affected by
# outliers such as the maximum, the minimum, the mean, and the standard deviation.


# def robust_scaling(df):
#     # apply robust scaling
#     for column in df_robust.columns:
#         df_robust[column] = (df_robust[column] - df_robust[column].median()) / (
#                     df_robust[column].quantile(0.75) - df_robust[column].quantile(0.25))
#     return df_robust

def robust_scaling(column):
    column = (column - column.median()) / (column.quantile(0.75) - column.quantile(0.25))
    return column

data_types = ["Transcriptomics", "Proteomics", "Metabolomics"]
datadic = {"Transcriptomics": None, "Proteomics": None, "Metabolomics": None}


if __name__ == '__main__':
    #process transcriptomic data
    data_type="Transcriptomics"
    files=get_files(data_type)
    j = files[0][0]
    path = str(files[1]) + str(j)
    mapping_dic = mapping_convertion(path)
    data_f=pd.read_excel(path , index_col=0)
    data_f=data_f.T
    data_f.rename(columns=data_f.iloc[0])
    #drop unwanted columns an rows
    data_f=data_f.drop(index=["Full Name","Fold Change"],axis=0)
    data_f.rename(columns=data_f.iloc[0])
    #mapping entrez
    data_f=map(data_f,mapping_dic)
    data_f = data_f.drop(data_f.columns[0], axis=1)
    for i in methods:
        res=pre_processing(data_f,i)
        res[0].to_excel(files[1] + str("proc/") + "10.1371journal.pone.0135426_proc_" + res[1] + ".xlsx", index=False)

    #file 2
    files
    j=files[0][1]
    path = str(files[1]) + str(j)
    mapping_dic = mapping_convertion(path)
    data_f=pd.read_excel(path , index_col=0)
    data_f=data_f.T
    data_f.rename(columns=data_f.iloc[0])
    #mapping entrez
    data_f=map(data_f,mapping_dic)
    data_f = data_f.drop(data_f.columns[0], axis=1)
    for i in methods:
        res=pre_processing(data_f,i)
        res[0].to_excel(files[1] + str("proc/") + "Supdata_proc_" + res[1] + ".xlsx", index=False)


    #join datasets

    files = [i for i in os.listdir(TRANSCRIPTOMICS_PATH+str("proc/")) if i.endswith(".xlsx")]
    omic_data={i: pd.read_excel(TRANSCRIPTOMICS_PATH+str("proc/")+str(i)) for i in files}
    compiled=pd.concat(omic_data, join="outer")
    compiled.to_excel(TRANSCRIPTOMICS_PATH+str("proc/") + "Compiled.xlsx", index=True)

    # #criar o dicionario de multiindexação
    mult=pd.MultiIndex.from_frame(compiled)