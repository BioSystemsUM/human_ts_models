### Modules:
from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import ast
import cobra
from functools import reduce
from EssTasksGenesMetbFunct import subsetCellTypes

### Functions:
def FVA(BaseDir, modelMedAdap):
    '''
    - does and saves Flux Variability Analysis (FVA) results on reconstructed models
    :param BaseDir: basic directory
    :param modelMedAdap: geneic metabolic model adapted for medium composition
    '''
    path = os.path.join(BaseDir, 'results/Reconst/RecModelBounds.tab')
    modelBds = pd.read_csv(path, sep='\t', index_col=[0])
    modelBds = modelBds.applymap(lambda x: ast.literal_eval(x))
    pathDir = os.path.join(BaseDir, 'results/Simulations/FVA')
    if not os.path.exists(pathDir):
        os.makedirs(pathDir)
    for mdN in modelBds:
        mdBds = modelBds[mdN]
        with modelMedAdap as md:
            for rcid, bd in zip(mdBds.index, mdBds):
                md.reactions.get_by_id(rcid).bounds = bd
            # do FVA for biomass from 0 to 100% with 10% intervals:
            frcL = list()
            for frc in np.arange(0, 1.1, 0.1):
                fva = cobra.flux_analysis.flux_variability_analysis(md, md.reactions, fraction_of_optimum=frc)  # dataframe where reactions are rows and columns are min and max FVA fluxes
                fva['mean'] = fva.mean(axis=1)  # add column with average flux (mean of min and max FVA fluxes)
                fva.index.name = 'Reaction_ID'  # add row index title
                fva.columns = map(lambda x: '_'.join([str(frc), x]), fva.columns)
                frcL.append(fva)
            mdDf = reduce(lambda x,y: x.join(y, how='outer'), frcL)
            mdDf.to_csv(os.path.join(pathDir, mdN), sep='\t')

def clustermap(data, path, csize, rsize, fn, xlabSiz, ylabSiz, hmypos, hmhfactor, hmxpos, hmwfactor, rdengx, cdengy, rdengwfactor, cdenghfactor, showvalues, rowdeng=True, coldeng=True, cbarx=0.4, cbarh=0.05):
    '''
    - plots clustermap of values where columns are tissues and rows are features,
    :param data: dataframe to plot
    :param path: path to final plot
    :param csize: factor to increase width of graph
    :param rsize: factor to increase heigh of graph
    :param fn: factor to increase letter size overal
    :param xlabSiz: fontsize of column labels
    :param ylabSiz: fontsize of row labels
    :param hmypos: value to add to current heatmap position to move it vertically
    :param hmhfactor: factor by which current heatmap height is multiplied to expand or shrink heatmap vertically
    :param hmxpos: value to add to current heatmap position to move it horizontally
    :param hmwfactor: factor by which current heatmap width is multiplied to expand or shrink heatmap horizontally
    :param rdengx: value to add to current row dendogram position to move dendrogram horizontally
    :param cdengy: value to add to current column dendogram position to move dendrogram vertically
    :param rdengwfactor: factor by which current row dendogram width is multiplied to expand or shrink dendrogram horizontally
    :param cdenghfactor: factor by which current column dendogram width is multiplied to expand or shrink dendrogram vertically
    :param showvalues: shows values in each cell of clustermap (True/False) or if we give a dataframe with values adds those values
    :param rowdeng: show row dendrogram (True/False) (it can still cluster even if dendrogram not shown)
    :param coldeng: show row dendrogram (True/False) (it can still cluster even if dendrogram not shown)
    :param cbarx: x position of color bar
    :param cbarh: height of color bar
    '''
    nc = data.shape[1]
    nr = data.shape[0]
    fcl = fn / nc
    sns.set(font_scale=fcl)
    colrs = sns.color_palette("dark:white_r", as_cmap=True)
    res = sns.clustermap(data, cmap=colrs, linecolor='black', xticklabels=True, yticklabels=True, linewidths=0.1, figsize=(csize * nc, rsize * nr), dendrogram_ratio=(0.2, 0.05), cbar_pos=(cbarx, 0.9, 0.01, cbarh), cbar_kws={'ticks': [int(data.min().min()), int(0), int(data.max().max())]}, annot=showvalues, fmt='d')
    res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_ymajorticklabels(), size=ylabSiz)
    res.ax_heatmap.set_xticklabels(res.ax_heatmap.get_xmajorticklabels(), size=xlabSiz, rotation=90)
    hm = res.ax_heatmap.get_position()
    col = res.ax_col_dendrogram.get_position()
    row = res.ax_row_dendrogram.get_position()
    res.ax_heatmap.set_position([hm.x0 + hmxpos, hm.y0 + hmypos, hm.width * hmwfactor, hm.height * hmhfactor])
    res.ax_row_dendrogram.set_position([row.x0 + rdengx, row.y0 + hmypos, row.width * rdengwfactor, row.height * hmhfactor])
    res.ax_col_dendrogram.set_position([col.x0 + hmxpos, col.y0 + cdengy, col.width * hmwfactor, col.height * cdenghfactor])
    if not rowdeng:
        res.ax_row_dendrogram.set_visible(False)
    if not coldeng:
        res.ax_col_dendrogram.set_visible(False)
    #plt.show()
    res.savefig(path, dpi=300)
    plt.close('all')

def onerowmap(data, path, csize, rsize, fn, xlabSiz, ylabSiz, hmypos, hmhfactor, hmxpos, hmwfactor, showvalues):
    '''
    - heatmap with one row
    - same parameters as function above
    '''
    nc = data.shape[1]
    nr = data.shape[0]
    fcl = fn / nc
    sns.set(font_scale=fcl)
    colrs = sns.color_palette("dark:white_r", as_cmap=True)
    res = sns.clustermap(data, cmap=colrs, linecolor='black', xticklabels=True, yticklabels=True, linewidths=0.1, figsize=(csize * nc, rsize * nr), row_cluster=False, col_cluster=False, cbar_pos=(0.4, 0.8, 0.01, 0.2), cbar_kws={'ticks': [int(data.min().min()), int(0), int(data.max().max())]}, annot=showvalues, fmt='d', rasterized=True)
    res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_ymajorticklabels(), size=ylabSiz, rotation=0)
    res.ax_heatmap.set_xticklabels(res.ax_heatmap.get_xmajorticklabels(), size=xlabSiz, rotation=90)
    hm = res.ax_heatmap.get_position()
    res.ax_heatmap.set_position([hm.x0 + hmxpos, hm.y0 + hmypos, hm.width * hmwfactor, hm.height * hmhfactor])
    #plt.show()
    res.savefig(path, dpi=300)
    plt.close('all')

def MiRsTFs(BaseDir, model, StdD, TissueCode):
    '''
    - saves tables with info on genes inversely/directly correlated with growth
    - saves tables and heatmaps with info on transcription factors(TFs) knockouts and miRNAs (miRs) that may potentially target genes directly correlated with biomass
    :param BaseDir: base directory
    :param model: generic model
    :param StdD: dict with tissue name and corresponding models
    :param TissueCode: dict with tissue code vs tissue name
    '''
    # get dataframe with reactions on rows and models on columns, with 'True' when reaction flux is inversely correlated with biomass percentage:
    D = dict()
    pathD = os.path.join(BaseDir, 'results/Simulations/FVA')
    pathF = os.path.join(BaseDir, 'results/Correlation')
    if not os.path.exists(pathF):
        os.makedirs(pathF)
    for f in os.listdir(pathD):
        df = pd.read_csv(os.path.join(pathD, f), sep='\t', index_col=0)
        # FVA of 100% biomass (= 1.0) must be excluded. The solver doesn't use precise values of biomass flux
        # (it rounds the value - cause of decimal to hexadecimal conversion)
        # for 100% biomass if solver rounds flux to a higher number, it will give a solution that when used afterwards is infeasible
        # cause in fact that value is higher than 100% biomass:
        df.drop(['1.0_minimum', '1.0_maximum', '1.0_mean'], axis=1, inplace=True)
        meanLbs = [c for c in df.columns if 'mean' in c] # column names corresponding to mean of fluxes
        df = df[meanLbs]
        y = pd.Series([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]) # perc of biomass
        y.index = df.columns
        D[f] = df.apply(lambda x: x.corr(y, method='pearson'), axis=1)
    mdf = pd.DataFrame.from_dict(D)
    inv = mdf < -0.7 # reactions with inverse correlation with biomass for each model
    dir = mdf > 0.7 # reactions correlated with biomass for each model
    # get similar dataframes as above but for genes corresponding to reactions inversely correlated with biomass:
    hgnc_df = pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt', index_col=0, sep='\t')  # gene id conversion dataframe
    hgnc_df = hgnc_df[['ensembl_gene_id', 'symbol']]
    hgnc_df.index = hgnc_df['ensembl_gene_id']
    hgnc_df.drop('ensembl_gene_id', axis=1, inplace=True)
    def selectCorr(df):
        AllGenes = [g.id for g in model.genes]  # list of genes in universal model
        genesD = dict()
        for md in df:
            gLst = list()
            for rid, v in zip(df.index, df[md]):
                if v:
                    rc = model.reactions.get_by_id(rid)
                    for g in rc.genes:
                        gLst.append(g.id)
            gLst = np.unique(gLst)
            mdGenesD = {gid: True if gid in gLst else False for gid in AllGenes}
            genesD[md] = mdGenesD  # dict with model id vs inner dict with gene id and True/False
        dfG = pd.DataFrame.from_dict(genesD)
        convD = hgnc_df['symbol'].to_dict()
        dfG.drop([g for g in dfG.index if g not in convD.keys()], inplace=True) # exclude ensembl genes that do not have corresponding symbol
        dfG.index = [convD[g] for g in dfG.index]
        inCSCnotCC, inCSCinCC, inCCnotCSC = subsetCellTypes(mdF = dfG, StdD=StdD, TissueCode=TissueCode)
        return inCSCnotCC, inCSCinCC, inCCnotCSC
    inCSCnotCC, inCSCinCC, inCCnotCSC = selectCorr(inv)
    def organize(df):
        return pd.DataFrame({ts: [', '.join([g for g, bool in v.items() if bool])] for ts, v in df.to_dict().items()}).T
    inCSCnotCC = organize(inCSCnotCC)
    inCSCnotCC.to_csv(os.path.join(pathF, 'Inverse_inCSCnotCC.tab'), sep='\t', header=False)
    inCSCinCC = organize(inCSCinCC)
    inCSCinCC.to_csv(os.path.join(pathF, 'Inverse_inCSCinCC.tab'), sep='\t', header=False)
    ## get similar dataframe as above but for genes corresponding to reactions positively correlated with biomass:
    inCSCnotCC, inCSCinCC, inCCnotCSC = selectCorr(dir)
    inCSCnotCC = organize(inCSCnotCC)
    inCSCnotCC.to_csv(os.path.join(pathF, 'Direct_inCSCnotCC.tab'), sep='\t', header=False)
    inCSCinCC = organize(inCSCinCC)
    inCSCinCC.to_csv(os.path.join(pathF, 'Direct_inCSCinCC.tab'), sep='\t', header=False)
    ## find MiRs that affect genes positively correlated with biomass:
    hsa_MTI = pd.read_excel('https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/8.0/hsa_MTI.xlsx') # old link: 'http://mirtarbase.cuhk.edu.cn/cache/download/8.0/hsa_MTI.xlsx'
    hgnc_df.reset_index(inplace=True)
    frD = {n: set(g['Target Gene'].dropna()) for n,g in hsa_MTI.groupby('miRNA')} # dict with miR name vs gene targets in database
    ## all methods of multitest p-value adjustment for hypergeometric test were tested and no miR was significant (adj p-val < 0.05) in any tissue
    ## so instead decided to count those mirs with more targets and use top 10:
    def topmiR(df, top, n):
        mirs = set(frD.keys()) # all miRs in database
        tgD = df.iloc[:,0].T.to_dict() # dict with tissue vs genes
        # get dict with tissue vs dict, inner dict is miR vs number of miR targets that exist in a tissue:
        intsD = dict()
        for ts,gs in tgD.items(): # for each tissue
            tsg = set(gs.split(', ')) # tissue genes
            intD = dict()
            for mir in mirs: # for each miR
                mg = frD[mir] # mir genes
                ints = tsg.intersection(mg)
                if ints:
                    intD[mir] = len(ints)
            intsD[ts] = intD
        # dataframe with miRs in rows, tissues as columns and values are number of overlaped targets:
        tsDf = pd.DataFrame(intsD)
        # exclude mirs with less than n targets:
        tsDf = tsDf[tsDf >= n]
        # find miRs with top 10 number of targets for each tissue,
        # it might be more than 10 miRs per tissue, as some mirs might have the same number of targets:
        tsD = dict()
        for ts in tsDf:
            tgN = np.sort(tsDf[ts].dropna().unique())[::-1][:top]
            tsD[ts] = tsDf[ts].apply(lambda x: int(x) if x in tgN else 0)
        topMirdf = pd.DataFrame(tsD)
        return topMirdf[topMirdf.any(axis=1)]  # exclude rows always 0/False
    top_inCSCnotCC = topmiR(df = inCSCnotCC, top = 10, n=10) # top 10 gives 30 rows
    clustermap(data = top_inCSCnotCC, path = os.path.join(pathF, 'miR_Direct_inCSCnotCC.svg'), csize = 1.4, rsize = 0.3, fn = 6, xlabSiz = 10, ylabSiz= 9, hmypos = 0.15, hmhfactor = 0.8, hmxpos = -0.145, hmwfactor = 0.3, rdengx = 0, cdengy = -0.014, rdengwfactor = 0.2, cdenghfactor = 0.95, showvalues = True, rowdeng=False, coldeng=False)
    top_inCSCinCC = topmiR(df = inCSCinCC, top = 10, n=10) # top 10 gives 35 rows
    clustermap(data = top_inCSCinCC, path = os.path.join(pathF, 'miR_Direct_inCSCinCC.svg'), csize = 1.4, rsize = 0.3, fn = 6, xlabSiz = 10, ylabSiz= 9, hmypos = 0.165, hmhfactor = 0.8, hmxpos = -0.145, hmwfactor = 0.3, rdengx = 0, cdengy = -0.014, rdengwfactor = 0.2, cdenghfactor = 0.4, showvalues = True, rowdeng=False, coldeng=False)
    ## find TFs that when ko downregulate genes positively correlated with biomass:
    # knockTF database has gene expression profiles when TFs are knock out:
    kTF = pd.read_csv('http://www.licpathway.net/KnockTF/download_FC/differential%20expression%20of%20genes%20in%20all%20datasets.txt', sep='\t')
    # get dict with TF vs target gene which expression decreases when TF is ko for TFs in database:
    kTF = kTF[kTF['Log2FC'] < -1.5]  # select targets which gene expression decreases when TF is ko
    frD = {n: set(g['Gene'].dropna()) for n, g in kTF.groupby('TF')}
    # hypergeometric test (https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458):
    M = len({t for v in frD.values() for t in v})  # population size: number of all genes activated by TFs
    TFs = set(frD.keys())  # all TFs in database
    def hypfct(df):
        tgD = df.iloc[:, 0].T.to_dict()  # dict with tissue vs genes
        tsD = dict()
        intsD = dict()
        for ts, gs in tgD.items():  # for each tissue
            tsg = set(gs.split(', '))  # tissue genes
            N = len(tsg)  # sample size: number of genes in each tissue
            tfD = dict()
            intD = dict()
            for tf in TFs:  # for each TF
                mg = frD[tf]  # TF genes
                n = len(mg)  # number of successes in the population
                X = len(tsg.intersection(mg))  # number of successes in the sample: number of genes from a tissue inactivated by a TF ko
                pval = hypergeom.sf(X - 1, M, n, N)
                int = tsg.intersection(mg)  # actual genes from a tissue inactivated by a TF ko
                tfD[tf] = pval
                intD[tf] = int
            tsD[ts] = tfD
            intsD[ts] = intD
        return tsD, intsD
    tsD, intsD = hypfct(df=inCSCinCC)
    tsDf = pd.DataFrame(tsD)  # dataframe with TFs in rows, tissues as columns and values are p-values
    tsDf = tsDf.apply(lambda x: multi.multipletests(x, method='fdr_bh')[1], axis=0)  # adjust p-values for multiple-testing
    tsDf[tsDf >= 0.05] = 1  # all non-significant TFs will get a pvalue of 1
    tsDf = tsDf.applymap(lambda x: -np.log(x))
    toKeep = tsDf.any(axis=1)  # all rows that have at least 1 value dif from 0
    tsDf = tsDf[toKeep]
    tsDfint = pd.DataFrame(intsD).applymap(len)  # dataframe with TFs in rows, tissues as columns and values are number of overlapped targets
    tsDfint = tsDfint[toKeep]
    clustermap(data=tsDf, path=os.path.join(pathF, 'TFko_Direct_inCSCinCC.svg'), csize=1.4, rsize=0.5, fn=6, xlabSiz=10,
               ylabSiz=9, hmypos=0.35, hmhfactor=0.5, hmxpos=-0.15, hmwfactor=0.3, rdengx=0, cdengy=-0.014,
               rdengwfactor=0.5, cdenghfactor=0.95, showvalues=tsDfint, rowdeng=False, coldeng=False, cbarx=0.35, cbarh=0.08)
    tsD, intsD = hypfct(df=inCSCnotCC)
    tsDf = pd.DataFrame(tsD)  # dataframe with TFs in rows, tissues as columns and values are p-values
    tsDf = tsDf.apply(lambda x: multi.multipletests(x, method='fdr_bh')[1], axis=0)  # adjust p-values for multipletesting.
    tsDf[tsDf >= 0.05] = 1  # all non-significant TFs will get a pvalue of 1
    tsDf = tsDf.applymap(lambda x: -np.log(x))
    toKeep = tsDf.any(axis=1)  # all rows that have at least 1 value dif from 0
    tsDf = tsDf[toKeep]
    tsDfint = pd.DataFrame(intsD).applymap(len)  # dataframe with TFs in rows, tissues as columns and values are number of overlapped targets
    tsDfint = tsDfint[toKeep]
    onerowmap(data = tsDf, path=os.path.join(pathF, 'TFko_Direct_inCSCnotCC.svg'), csize=1.4, rsize=0.8, fn=6, xlabSiz=10, ylabSiz=9, hmypos=0.3, hmhfactor=1, hmxpos=-0.15, hmwfactor=0.4, showvalues=tsDfint)
