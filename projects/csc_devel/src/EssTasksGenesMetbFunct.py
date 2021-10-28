### Modules:
import os
import pandas as pd
import ast
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from cobamp.wrappers.external_wrappers import get_model_reader
from LethalityEvalFunct import InactReactFromGeneKo
from multiprocessing import cpu_count
from LethalityEvalFunct import batch_ko_optimize
from troppo.tasks.core import TaskEvaluator
from collections import Counter
from statsAndSimFunct import medianMdSameStd, applyclustermap, subsetCellTypes
from functools import reduce
import xml.etree.ElementTree as ET


### Functions:
def PtGeneTargets(BaseDir, modelMedAdap, task_list):
    '''
    - identify genes that are essential for reconstructed model growth and do not affect biomass or essential tasks in generic medium adapted model
    :param BaseDir: base directory
    :param modelMedAdap: generic model adapted for medium composition
    :param task_list: list of 57 essential tasks
    :return:
    '''
    p = os.path.join(BaseDir, 'support/Essentiality')
    NTHREADS = cpu_count() - cpu_count() / 4
    ## save essential genes and tasks completed with gene kos for generic metabolic model adapted for medium composition (if not already done):
    if not os.path.exists(os.path.join(p, 'TasksworkingWithGeneKo.tab')):
        # get cobamp model from cobrapy model:
        cobamp_model = get_model_reader(modelMedAdap, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
        # all genes in cobamp model:
        model_genes = cobamp_model.gpr.get_genes()
        # get biomass flux (objective value) on wild-type:
        sol = cobamp_model.optimize({'biomass_human': 1}, False)
        ofv = sol.objective_value()
        # get inactive reactions when a gene is knockout in generic model adapted for medium composition:
        kos_to_test, reaction_states_inactive = InactReactFromGeneKo(cobamp_model, model_genes)
        unique_ko_dict_list = [{i:(0,0) for i in k} for k in kos_to_test] # list of dicts. each dict represents a gene ko and contains reaction:(0,0)
        # ko combinations of inactive reactions and get simulated lethality score (obj. value of ko model / obj. value of wild-type model):
        bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context={}, objective={'biomass_human': 1}, objective_sense=False, threads=NTHREADS)
        # create dict where key is comb of reactions ko and value is score representing if that combination was essential (0) or affected biomass (<1) or not (>=1):
        opt_res = dict(zip(kos_to_test, bkopt))
        of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in reaction_states_inactive.items()}
        # if else condition  guarantees that genes
        # which don't have associated inactive reaction when ko get a score of 1 - meaning they don't affect biomass
        # threshold for considering a gene lethal (when biomass o model with ko < 0.1% of biomass of wild-type model):
        t = 0.001*ofv
        # dict where essential genes (lethal genes) in generic medium adapted model are 'True':
        essGen = {g: [True] if s < t else [False] for g, s in of_res.items()}
        # save dataframe with lethal genes in generic medium adapted model:
        if not os.path.exists(p):
            os.makedirs(p)
        pd.DataFrame.from_dict(essGen).to_csv(os.path.join(p, 'EssGenesGenericAdapMd.tab'), sep='\t', index=False)
        # get a task evaluator object for generic medium adapted model:
        modelMedAdap_copy = modelMedAdap.copy()
        for r in modelMedAdap_copy.boundary: # close boundary reactions to test tasks
            r.knock_out()
        cobamp_model_copy = get_model_reader(modelMedAdap_copy, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
        tev = TaskEvaluator(model=cobamp_model_copy, tasks=task_list)
        # get dataframe for a generic medium adapted model, each row one gene ko, and each column one task:
        ko_task_evaluation = tev.batch_evaluate(bound_changes=unique_ko_dict_list, threads=int(NTHREADS), output_sol=False, mp_batch_size=50000)
        taskGen = pd.DataFrame([{k: v[0] for k, v in ko_task_evaluation[i].items()} for i in range(len(unique_ko_dict_list))])
        rct2gene = {r:g for g,r in reaction_states_inactive.items()} # dict with inactive reactions as keys and gene ko as values
        taskGen.index = [rct2gene[el] for el in kos_to_test]
        # save dataframe with task completion (True when task complete) info for generic medium adapted model:
        taskGen.to_csv(os.path.join(p, 'TasksworkingWithGeneKo.tab'), sep='\t')
    ## exclude genes essential for growth or essential tasks in generic model from kos to test in reconstructed models:
    essGen = pd.read_csv(os.path.join(p, 'EssGenesGenericAdapMd.tab'), sep='\t') # loads essential/lethal genes in generic adapted model
    taskGen = pd.read_csv(os.path.join(p, 'TasksworkingWithGeneKo.tab'), sep='\t', index_col=[0]) # loads tasks working when a gene is ko in generic adapted model
    # get cobamp model from cobrapy model:
    cobamp_model = get_model_reader(modelMedAdap, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
    # all genes in cobamp model:
    model_genes = cobamp_model.gpr.get_genes()
    # get inactive reactions when a gene is knockout in generic model adapted for medium composition:
    kos_to_test, reaction_states_inactive = InactReactFromGeneKo(cobamp_model, model_genes)
    # exclude genes essential for growth or essential tasks in generic model from kos to test in reconstructed model:
    exc = [g for g in essGen if essGen[g][0]] # essential genes for growth in generic md
    exc = exc + list(taskGen.index[taskGen.sum(axis=1) != taskGen.shape[1]]) # essential genes for tasks in generic md
    exc = list(set(exc)) # get unique genes
    excr = [reaction_states_inactive[g] for g in exc] # groups of ko reactions corresponding to gene kos to be excluded
    kos_to_test = [ko for ko in kos_to_test if ko not in excr]
    # reconstruct models and test groups of reaction kos corresponding to each gene:
    path = os.path.join(BaseDir, 'results/Reconst/RecModelBounds.tab')
    modelBds = pd.read_csv(path, sep='\t', index_col=[0])
    RecGeneD = dict()
    for mdN in modelBds:
        mdBds = modelBds[mdN].to_dict()
        bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context=mdBds, objective={'biomass_human': 1}, objective_sense=False, threads=NTHREADS)
        with cobamp_model as m:
            for crx, Bd in mdBds.items():
                Bd = ast.literal_eval(Bd)  # values are as str so eval is needed
                m.set_reaction_bounds(crx, lb=Bd[0], ub=Bd[1])  # change reaction bounds
            solrec = m.optimize({'biomass_human': 1}, False)
            ofvrec = solrec.objective_value()
        # create dict where key is comb of reactions ko and value is score representing if that combination was essential (0) or affected biomass (<1) or not (>=1):
        opt_res = dict(zip(kos_to_test, bkopt))
        of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in reaction_states_inactive.items()}
        # if else condition  guarantees that genes
        # which don't have associated inactive reaction when ko or that are excluded genes (see above) get a score of 1 - meaning they don't affect biomass
        # threshold for considering a gene lethal (when biomass o model with ko < 0.1% of biomass of reconstructed wild-type model):
        t = 0.001*ofvrec
        # dict where essential genes (lethal genes) in reconstructed model that do Not affect biomass and tasks on gneric model are 'True':
        essGen = {g: True if s < t else False for g, s in of_res.items()}
        RecGeneD[mdN] = essGen
    # save dataframe with potential target genes for reconstructed models:
    p = os.path.join(BaseDir, 'results/Essentiality')
    if not os.path.exists(p):
        os.makedirs(p)
    pd.DataFrame.from_dict(RecGeneD).to_csv(os.path.join(p, 'PtEssGenesRec.tab'), sep='\t')

def metb2RcBds(modelMedAdap):
    '''
    - get dictionary with metabolite names as keys and corresponding reactions which bounds have to be changed to, to simulate metabolite ko
    :param modelMedAdap: generic model adapted for medium composition
    '''
    # get all metabolite names (common between same metabolites in different compartments):
    mtbN = {m.name for m in modelMedAdap.metabolites}
    # get dict mapping metabolite name to list of metabolite ids:
    mtbNm2IdD = dict()
    for name in mtbN:
        mtLst = list()
        for m in modelMedAdap.metabolites:
            if m.name == name:
                mtLst.append(m.id)
        mtbNm2IdD[name] = mtLst
    # get dict where key is metabolite name and value is tuple of tuples,
    # outer tuple contains tupples, where each inner tuple represents a reaction where that metabolite is a reagent/substrate
    # inner tuple is (metaboliteID+compartment combination, lb, ub).
    # lb and ub are bounds that reaction has to have so that the effect of metabolite in model is abrogated:
    mtb2RcD = dict()
    for k, v in mtbNm2IdD.items(): # for each metabolite
        mtbT = list()
        for mid in v: # for each metabolite + compartment combination
            m = modelMedAdap.metabolites.get_by_id(mid)
            for r in m.reactions: # for each reaction of metabolite+compartment combination
                if m in r.reactants:
                    rs = (r.id, r.lower_bound, 0.0)
                elif m in r.products:
                    rs = (r.id, 0.0, r.upper_bound)
                mtbT.append(rs)
        # sometimes the same reaction has two/more combinations of metabolite+compartment of same metabolite,
        # one combination in reagents and the other in the products. in those cases, that reaction appears twice/more in 'mtbT',
        # one time with bounds (-1000, 0) and the other with (0, 1000).
        # code bellow, removes the duplicated reaction and replaces by same reaction with bounds (0,0):
        dupRc = [k for k,v in Counter([r[0] for r in mtbT]).items() if v > 1] # find duplc reactions
        mtbT = [rcT for rcT in mtbT if rcT[0] not in dupRc] + [(dp, 0.0, 0.0) for dp in dupRc]
        # add tuple with reactions where metabolite is reagent to dict where key is metabolite name:
        mtb2RcD[k] = tuple(mtbT)
    return mtb2RcD

def PtMtbTargets(BaseDir, modelMedAdap, task_list):
    '''
    - identify metabolites that are essential for reconstructed model growth and do not affect biomass or essential tasks in generic medium adapted model
    :param BaseDir: base directory
    :param modelMedAdap: generic model adapted for medium composition
    :param task_list: list of 57 essential tasks
    '''
    p = os.path.join(BaseDir, 'support/Essentiality')
    #NTHREADS = cpu_count() - cpu_count() / 4
    NTHREADS=10
    ## save essential genes and tasks completed with gene kos for generic metabolic model adapted for medium composition (if not already done):
    if not os.path.exists(os.path.join(p, 'TasksworkingWithMetbKo.tab')):
        # get cobamp model from cobrapy model:
        cobamp_model = get_model_reader(modelMedAdap, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
        # get biomass flux (objective value) on wild-type:
        sol = cobamp_model.optimize({'biomass_human': 1}, False)
        ofv = sol.objective_value()
        # get list of reactions+reaction bounds when a metabolite is knockout in generic model adapted for medium composition:
        mtb2RcD = metb2RcBds(modelMedAdap) # dict with metabolite names as keys and corresponding reactions which bounds have to be changed to, to simulate metabolite ko
        kos_to_test = list(mtb2RcD.values())  # list of tuples, each outer tuple in list corresponds to one metabolite, each inner tuple to a inactive reaction
        unique_ko_dict_list = [{i[0]:i[1:] for i in k} for k in kos_to_test] # list of dicts. each dict represents a metabolite ko and contains reaction bounds to change
        # ko a metabolite and get simulated lethality score (obj. value of ko model / obj. value of wild-type model):
        bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context={}, objective={'biomass_human': 1}, objective_sense=False, threads=NTHREADS, EssMtb=True)
        # create dict where key is comb of reactions to change bounds and value is score representing if that combination was essential (0) or affected biomass (<1) or not (>=1):
        opt_res = dict(zip(kos_to_test, bkopt))
        of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in mtb2RcD.items()}
        # above if else condition guarantees that metabolites
        # which don't have associated reaction bounds to change when ko get a score of 1 - meaning they don't affect biomass
        # threshold for considering a metabolite lethal (when biomass o model with ko < 0.1% of biomass of wild-type model):
        t = 0.001*ofv
        # dict where essential metabolites (lethal metabolites) in generic medium adapted model are 'True':
        essMtb = {g: [True] if s < t else [False] for g, s in of_res.items()}
        # save dataframe with lethal metabolites in generic medium adapted model:
        if not os.path.exists(p):
            os.makedirs(p)
        pd.DataFrame.from_dict(essMtb).to_csv(os.path.join(p, 'EssMtbsGenericAdapMd.tab'), sep='\t', index=False)
        # get a task evaluator object for generic medium adapted model:
        modelMedAdap_copy = modelMedAdap.copy()
        for r in modelMedAdap_copy.boundary: # close boundary reactions to test tasks
            r.knock_out()
        cobamp_model_copy = get_model_reader(modelMedAdap_copy, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
        tev = TaskEvaluator(model=cobamp_model_copy, tasks=task_list)
        # get dataframe for a generic medium adapted model, each row one metabolite ko, and each column one task:
        ko_task_evaluation = tev.batch_evaluate(bound_changes=unique_ko_dict_list, threads=int(NTHREADS), output_sol=False, mp_batch_size=50000)
        taskMetb = pd.DataFrame([{k: v[0] for k, v in ko_task_evaluation[i].items()} for i in range(len(unique_ko_dict_list))])
        rct2metb = {r: g for g,r in mtb2RcD.items()} # dict with reactions with changed bounds as keys and metabolite ko as values
        taskMetb.index = [rct2metb[el] for el in kos_to_test]
        # save dataframe with task completion (True when task complete) info for generic medium adapted model:
        taskMetb.to_csv(os.path.join(p, 'TasksworkingWithMtbKo.tab'), sep='\t')
    ## exclude metabolites essential for growth or essential tasks in generic model from kos to test in reconstructed models:
    essMtb = pd.read_csv(os.path.join(p, 'EssMtbsGenericAdapMd.tab'), sep='\t') # loads essential/lethal metabolites in generic adapted model
    taskMtb = pd.read_csv(os.path.join(p, 'TasksworkingWithMtbKo.tab'), sep='\t', index_col=[0]) # loads tasks working when a metabolite is ko in generic adapted model
    # get cobamp model from cobrapy model:
    cobamp_model = get_model_reader(modelMedAdap, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
    # get list of reactions+reaction bounds when a metabolite is knockout in generic model adapted for medium composition:
    mtb2RcD = metb2RcBds(modelMedAdap)
    kos_to_test = list(mtb2RcD.values())  # list of tuples, each outer tuple in list corresponds to one metabolite, each inner tuple to a inactive reaction
    unique_ko_dict_list = [{i[0]: i[1:] for i in k} for k in kos_to_test]  # list of dicts. each dict represents a metabolite ko and contains reaction bounds to change
    # exclude metabolites essential for growth or essential tasks in generic model from kos to test in reconstructed model:
    exc = [g for g in essMtb if essMtb[g][0]] # essential metabolites for growth in generic md
    exc = exc + list(taskMtb.index[taskMtb.sum(axis=1) != taskMtb.shape[1]]) # essential metabolites for tasks in generic md
    exc = list(set(exc)) # get unique metabolites
    excr = [mtb2RcD[g] for g in exc] # groups of changed reactions corresponding to metabolite kos to be excluded
    kos_to_test = [ko for ko in kos_to_test if ko not in excr]
    # reconstruct models and test groups of afected reactions corresponding to each metabolite:
    path = os.path.join(BaseDir, 'results/Reconst/RecModelBounds.tab')
    modelBds = pd.read_csv(path, sep='\t', index_col=[0])
    RecMtbD = dict()
    for mdN in modelBds:
        mdBds = modelBds[mdN].to_dict()
        bkopt = batch_ko_optimize(model=cobamp_model, kos=kos_to_test, context=mdBds, objective={'biomass_human': 1}, objective_sense=False, threads=NTHREADS, EssMtb=True)
        with cobamp_model as m:
            for crx, Bd in mdBds.items():
                Bd = ast.literal_eval(Bd)  # values are as str so eval is needed
                m.set_reaction_bounds(crx, lb=Bd[0], ub=Bd[1])  # change reaction bounds
            solrec = m.optimize({'biomass_human': 1}, False)
            ofvrec = solrec.objective_value()
        # create dict where key is comb of reactions changed and value is score representing if that combination was essential (0) or affected biomass (<1) or not (>=1):
        opt_res = dict(zip(kos_to_test, bkopt))
        of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in mtb2RcD.items()}
        # above if else condition guarantees that metabolitesinCSCnotCC
        # which don't have associated reaction bounds to change when ko get a score of 1 - meaning they don't affect biomass
        # threshold for considering a metabolite lethal (when biomass o model with ko < 0.1% of biomass of wild-type model):
        t = 0.001*ofvrec
        # dict where essential metabolites (lethal metabolites) in reconstructed model that do Not affect biomass and tasks on gneric model are 'True':
        essMtb = {g: True if s < t else False for g, s in of_res.items()}
        RecMtbD[mdN] = essMtb
    # save dataframe with potential target metabolites for reconstructed models:
    p = os.path.join(BaseDir, 'results/Essentiality')
    if not os.path.exists(p):
        os.makedirs(p)
    pd.DataFrame.from_dict(RecMtbD).to_csv(os.path.join(p, 'PtEssMtbRec.tab'), sep='\t')

def heatmap(data, path, csize, rsize, fn, xlabSiz, ylabSiz, hmypos, hmhfactor, hmxpos, hmwfactor):
    '''
    - plots heatmap of booleans where columns are tissues and rows are features,
    with no cluster in rows (rows are sorted from those in higher number of tissues to less number of tissues)
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
    '''
    data = data.loc[data.sum(axis=1).sort_values(ascending=False).index]
    nc = data.shape[1]
    nr = data.shape[0]
    fcl = fn / nc
    sns.set(font_scale=fcl)
    colrs = sns.color_palette("dark:white_r", as_cmap=True)
    res = sns.clustermap(data, cmap=colrs, linecolor='black', xticklabels=True, yticklabels=True,
                         linewidths=0.1, figsize=(csize * nc, rsize * nr), row_cluster=False, col_cluster=True, dendrogram_ratio=(0.001),
                         cbar_pos=(0.9, 0.92, 0.02, 0.05))
    res.cax.set_visible(False)
    res.ax_row_dendrogram.set_visible(False)
    res.ax_col_dendrogram.set_visible(False)
    res.ax_heatmap.set_xticklabels(res.ax_heatmap.get_xmajorticklabels(), size=xlabSiz, rotation=90)
    res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_ymajorticklabels(), size=ylabSiz)
    hm = res.ax_heatmap.get_position()
    res.ax_heatmap.set_position([hm.x0 + hmxpos, hm.y0 + hmypos, hm.width * hmwfactor, hm.height * hmhfactor])
    #plt.tight_layout()
    #plt.show()
    #matplotlib.spines.Spine(res.ax_heatmap,'Line-like', matplotlib.spines.Spine.get_path, **{'color':'blue'})
    res.savefig(path, dpi=300)
    plt.close('all')

def EssGeneAnal(BaseDir, StdD, TissueCode):
    '''
    - creates heatmaps/tables with genes essential for one cell type and not the other
    :param BaseDir: basic directory
    :param StdD: dict with tissue name and corresponding models to use
    :param TissueCode: dict with tissue code vs tissue name
    '''
    # get dict with ensembl gene id vs symbol:
    hgnc_df = pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt', index_col=0, sep='\t')  # gene id conversion dataframe
    hgnc_df = hgnc_df[['ensembl_gene_id', 'symbol']]
    hgnc_df.index = hgnc_df['ensembl_gene_id']
    hgnc_df.drop('ensembl_gene_id', axis=1, inplace=True)
    convD = hgnc_df['symbol'].to_dict()
    # load dataframe with gene vs reconstructed model with 'Trues' when gene is essential in that model:
    mdDfIn = pd.read_csv(os.path.join(BaseDir, 'results/Essentiality', 'PtEssGenesRec.tab'), sep='\t', index_col=[0])
    mdDfIn = mdDfIn[mdDfIn.sum(axis=1) != 0] # exclude genes with all False, for all models
    mdDfIn.index = [convD[i] for i in mdDfIn.index]
    inCSCnotCC, inCSCinCC, inCCnotCSC = subsetCellTypes(mdF=mdDfIn, StdD=StdD, TissueCode=TissueCode)
    # exclude 'SLC22A5' as essential gene. This gene codes for OCTN2, a cytoplasmatic transporter of carnitine.
    # Since our models have all exchange reactions closed except the ones corresponding to entry of medium components,
    # we wouldn't expect models to have flux through the uptake of carnitine (as there is no carnitine in extracelular compartment).
    # the reason why 'SLC22A5' was simulated as essential was that it was incorrectly annotated to the transport reaction of carnitine in-between
    # the mitochondria and cytoplasm. the situation was signaled to the Human1 developers and fixed.:
    inCSCnotCC.drop(index='SLC22A5', inplace=True)
    inCSCinCC.drop(index='SLC22A5', inplace=True)
    # do clustermap without dendrograms:
    path = os.path.join(BaseDir, 'results/Essentiality/EssGen_inCSCnotCC.svg')
    heatmap(data=inCSCnotCC, path=path, csize=1*1.5, rsize=0.5*1.5, fn=18, xlabSiz=50, ylabSiz=32, hmypos=0.12, hmhfactor=0.95*0.9, hmxpos=0, hmwfactor=0.9)
    path = os.path.join(BaseDir, 'results/Essentiality/EssGen_inCSCinCC.svg')
    heatmap(data=inCSCinCC, path=path, csize=0.5, rsize=0.3, fn=18, xlabSiz=12, ylabSiz=10, hmypos=0.025, hmhfactor=1, hmxpos = 0, hmwfactor = 1)

def addinfo(subsD, mtb2RcD, dfm, ts, modelMedAdap, activeRcts):
    '''
    - adds corresponding subsystem, reactions and genes info to metabolite names
    :param subsD: dict with reaction id vs subsystems
    :param mtb2RcD: dict with metabolite vs tuple where elements are reactions that use that metabolite as substrate.
                    each element is a tuple: (reaction id, lb, ub) where lb and ub are bounds that simulate ko of metabolite
    :param dfm: dataframe with metabolite in rows and tissues in columns
    :param ts: tissue name
    :param modelMedAdap: generic model adapted for medium composition
    :param activeRcts: dataframe with active reactions (with flux) for each tissue
    '''
    s = pd.Series(dfm.index[dfm[ts]]) # metabolites
    actR = list(activeRcts.index[activeRcts[ts]]) # list with reactions active in a tissue
    dfa = pd.DataFrame.from_dict({mt: [', '.join(set([subsD[r[0]] for r in mtb2RcD[mt]]))] for mt in s}).T
    dfa.columns = ['Subsystems']
    dfb = pd.DataFrame.from_dict({mt: [', '.join(set(r[0] for r in mtb2RcD[mt] if r[0] in actR))] for mt in s}).T
    dfb.columns = ['Reactions']
    dfc = pd.DataFrame.from_dict({mt: [', '.join({g.id for r in mtb2RcD[mt] if (r[0] in actR) for g in modelMedAdap.reactions.get_by_id(r[0]).genes})] for mt in s}).T
    dfc.columns = ['Genes']
    return dfa.join(dfb, how='outer').join(dfc, how='outer')

def EssMtbAnal(BaseDir, StdD, TissueCode):
    '''
    - creates heatmaps/tables with metabolites essential for one cell type and not the other
    :param BaseDir: basic directory
    :param StdD: dict with tissue name and corresponding models to use
    :param TissueCode: dict with tissue code vs tissue name
    '''
    # load dataframe with metabolite vs reconstructed model with 'Trues' when metabolite is essential in that model:
    mdDfIn = pd.read_csv(os.path.join(BaseDir, 'results/Essentiality', 'PtEssMtbRec.tab'), sep='\t', index_col=[0])
    mdDfIn = mdDfIn[mdDfIn.sum(axis=1) != 0] # exclude genes with all False, for all models
    inCSCnotCC, inCSCinCC, inCCnotCSC = subsetCellTypes(mdF=mdDfIn, StdD=StdD, TissueCode=TissueCode)
    # do clustermap without dendrograms:
    path = os.path.join(BaseDir, 'results/Essentiality/EssMtb_inCSCnotCC.png')
    heatmap(data=inCSCnotCC, path=path, csize=1.7, rsize=0.55, fn=18, xlabSiz=28, ylabSiz=24, hmypos= 0.025, hmhfactor=0.95, hmxpos = 0, hmwfactor = 0.8)
    sort_inCSCnotCC = inCSCnotCC.loc[inCSCnotCC.sum(axis=1).sort_values(ascending= False).index]
    sort_inCSCnotCC.to_csv(os.path.join(BaseDir, 'results/Essentiality/EssMtb_inCSCnotCCOrder.tab'),sep='\t')
    path = os.path.join(BaseDir, 'results/Essentiality/EssMtb_inCSCinCC.png')
    heatmap(data=inCSCinCC, path=path, csize=1.7, rsize=0.55, fn=18, xlabSiz=28, ylabSiz=24, hmypos=0.025, hmhfactor=0.95, hmxpos = 0, hmwfactor = 0.79)
    sort_inCSCinCC = inCSCinCC.loc[inCSCinCC.sum(axis=1).sort_values(ascending=False).index]
    sort_inCSCinCC.to_csv(os.path.join(BaseDir, 'results/Essentiality/EssMtb_inCSCinCCOrder.tab'), sep='\t')

def antimetabolites(BaseDir, modelMedAdap, StdD, TissueCode):
    '''
    - Identify Antimetabolites that may abrogate the effect of essential metabolites
    :param BaseDir: base directory
    :param modelMedAdap: generic model adapted for medium composition
    :param StdD: dict with tissue name and corresponding models to use
    :param TissueCode: dict with tissue code vs tissue name
    '''
    base = ET.parse(os.path.join(BaseDir, 'data/drugbank_database/full_database.xml'))
    root = base.getroot()
    # get list of antimetabolites in drugbank database:
    lk = '{http://www.drugbank.ca}'
    AntLst = list()
    for drug in root:
        for cat in drug.find(lk+'categories').findall(lk+'category'):
            if 'Antimetabolites' in cat.find(lk+'category').text:
                AntLst.append(drug.find(lk+'name').text)
    # get dictionary with antimetabolite id in drugbank as key and corresponding targets, enzymes, carriers, transporters as values:
    mtbD = dict()
    for m in AntLst:
        for drug in root.findall(lk+'drug'):
            name = drug.find(lk+'name').text
            if name == m: # if drug is an antimetabolite
                mtLst = list()
                for e in drug.find(lk+'enzymes'): # for each enzyme corresponding to an antimetabolite
                    for pep in e.findall(lk+'polypeptide'): # for each peptide corresponding to a enzyme
                        f = pep.find(lk + 'external-identifiers').find(lk + 'external-identifier').find(lk + 'identifier').text
                        if 'HGNC' in f: # if enzyme is in human
                            mtLst.append(f)
                for tg in drug.find(lk+'targets'): # for each target corresponding to an antimetabolite
                    for pep in tg.findall(lk+'polypeptide'): # for each peptide corresponding to a target
                        f = pep.find(lk + 'external-identifiers').find(lk + 'external-identifier').find(lk + 'identifier').text
                        if 'HGNC' in f:  # if target is in human
                            mtLst.append(f)
                for c in drug.find(lk+'carriers'): # for each carrier corresponding to an antimetabolite
                    for pep in c.findall(lk+'polypeptide'): # for each peptide corresponding to a carrier
                        f = pep.find(lk + 'external-identifiers').find(lk + 'external-identifier').find(lk + 'identifier').text
                        if 'HGNC' in f:  # if carrier is in human
                            mtLst.append(f)
                for t in drug.find(lk+'transporters'): # for each transporter corresponding to an antimetabolite
                    for pep in t.findall(lk+'polypeptide'): # for each peptide corresponding to a transporter
                        f = pep.find(lk + 'external-identifiers').find(lk + 'external-identifier').find(lk + 'identifier').text
                        if 'HGNC' in f:  # if transporter is in human
                            mtLst.append(f)
                if set(mtLst): # if set is not empty
                    mtbD[m] = set(mtLst)
    # target vs antimetabolite dict:
    TAntM = {el:k for k,v in mtbD.items() for el in v}
    # prepare to convert essencial metabolites' gene ensembl id to HGNC gene id:
    hgnc_df = pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt', sep='\t')  # gene id conversion dataframe
    hgnc_df = hgnc_df[['ensembl_gene_id', 'hgnc_id']]
    hgnc_df.index = hgnc_df['ensembl_gene_id']
    hgnc_df.drop('ensembl_gene_id', axis=1, inplace=True)
    convD = hgnc_df['hgnc_id'].to_dict() # dict with ensembl gene id vs hgnc gene id
    # get dataframe with essential metabolite vs anti-metabolite:
    def EssMtbAntD(df):
        EssD = {k: [convD[g] for g in v.split(', ') if g in convD.keys()] for k,v in df['Genes'].to_dict().items()} # dict with essential metabolite vs list of hgnc gene ids corresponding to reactions afected by that metabolite
        D = {mt: {TAntM[g] for g in gLst if g in TAntM.keys()} for mt, gLst in EssD.items()} # dict with essential metabolite vs antimetabolite
        df = pd.DataFrame.from_dict({k: [', '.join(v)] for k,v in D.items()}).T
        df.columns = ['Antimetabolite']
        df.index.name = 'Essential Metabolite'
        return df
    ## convert essential metabolite to antimetabolite in inCSCnotCC and inCSCinCC dataframes:
    # load dataframe with metabolite vs reconstructed model with 'Trues' when metabolite is essential in that model:
    mdDfIn = pd.read_csv(os.path.join(BaseDir, 'results/Essentiality', 'PtEssMtbRec.tab'), sep='\t', index_col=[0])
    mdDfIn = mdDfIn[mdDfIn.sum(axis=1) != 0] # exclude genes with all False, for all models
    inCSCnotCC, inCSCinCC, inCCnotCSC = subsetCellTypes(mdF=mdDfIn, StdD=StdD, TissueCode=TissueCode)
    # add subsystem, reactions, genes info:
    subsD = pd.read_csv(os.path.join(BaseDir, 'support/Pathways/subsytemInfo.tab'), sep='\t', index_col=0)['Subsystem'].to_dict()
    mtb2RcD = metb2RcBds(modelMedAdap)  # dict with metabolite names as keys and corresponding reactions which bounds have to be changed to, to simulate metabolite ko
    # identify active reactions (with flux) in each model:
    path = os.path.join(BaseDir, 'results/Reconst/RecModelBounds.tab')
    modelBds = pd.read_csv(path, sep='\t', index_col=[0])
    modelBds = modelBds.applymap(lambda x: ast.literal_eval(x))
    actRD = dict()
    for mdN in modelBds:
        mdBds = modelBds[mdN]
        with modelMedAdap as gmd:
            for rcid, bd in zip(mdBds.index, mdBds):
                if rcid in [r.id for r in gmd.reactions]:
                    gmd.reactions.get_by_id(rcid).bounds = bd
            opt = gmd.optimize()
            flx = opt.fluxes
            actRD[mdN] = ((flx < -0.5) | (flx > 0.5)).to_dict()
    actRDf = pd.DataFrame.from_dict(actRD)
    ActRcinCSCnotCC, ActRcinCSCinCC, ActRcinCCnotCSC = subsetCellTypes(mdF=actRDf, StdD=StdD, TissueCode=TissueCode)
    # dataframe with antimetabolites, corresponding essential metabolites and essentiality in each tissue:
    def getDf(dfm, activeRcts):
        # get dict with tissue vs dataframe with essential metabolites in that tissue and corresponding antimetabolites:
        d = {ts: EssMtbAntD(addinfo(subsD, mtb2RcD, dfm, ts, modelMedAdap, activeRcts)) for ts in dfm}
        # dataframe with antimetabolites, corresponding essential metabolites and essentiality in each tissue:
        for ts, df in d.items():
            df[ts] = df.shape[0]*[True]
            df.reset_index(inplace=True)
        ff = reduce(lambda x,y: pd.merge(x,y, how='outer', on=['Essential Metabolite', 'Antimetabolite']), d.values()).fillna(False)
        ff = ff[ff['Antimetabolite'] != ''] # exclude rows with no antimetabolite (just essential metabolite)
        return ff
    ff1 = getDf(inCSCnotCC, ActRcinCSCnotCC)
    ff1.to_csv(os.path.join(BaseDir, 'results/Essentiality/AntMtb_inCSCnotCC.tab'), sep='\t', index=None)
    ff2 = getDf(inCSCinCC, ActRcinCSCinCC)
    ff2.to_csv(os.path.join(BaseDir, 'results/Essentiality/AntMtb_inCSCinCC.tab'), sep='\t', index=None)
    # heatmap with anti-metabolites (no essential metabolites) and potential activity (True) or not in each tissue:
    def onlyAntm(f):
        antm = set([mt for el in f['Antimetabolite'].apply(lambda x: x.split(', ')) for mt in el])
        D = dict()
        for ps in range(2, f.shape[1]): # for each tissue
            tsL = list()
            for antL, r in zip(f['Antimetabolite'], f.iloc[:,ps]):
                if r: # if antimetabolite potencially works
                    tsL.extend(antL.split(', ')) # save antimetabolite name
            D[f.columns[ps]] =  {el: True if el in set(tsL) else False for el in antm}
        return pd.DataFrame.from_dict(D)
    inCSCnotCC = onlyAntm(ff1)
    inCSCinCC = onlyAntm(ff2)
    path = os.path.join(BaseDir, 'results/Essentiality/AntMtbOnly_inCSCnotCC.svg')
    heatmap(data=inCSCnotCC, path=path, csize=0.5, rsize=0.3, fn=18, xlabSiz=12, ylabSiz=10, hmypos=-0.2, hmhfactor=2.1, hmxpos=0, hmwfactor=1)
    path = os.path.join(BaseDir, 'results/Essentiality/AntMtbOnly_inCSCinCC.svg')
    heatmap(data=inCSCinCC, path=path, csize=0.5, rsize=0.75, fn=18, xlabSiz=12, ylabSiz=10, hmypos=0.025, hmhfactor=0.4, hmxpos=0, hmwfactor=1)

