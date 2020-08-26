from projects.vasco_proj.src.Paths import *
from cobra.io import read_sbml_model, load_matlab_model  # ler o modelo com  o cobra
from cobra.io import write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions, flux_variability_analysis, find_essential_reactions
from cobra.medium.minimal_medium import minimal_medium
from cobra.flux_analysis import pfba
from cobra.util.solver import linear_reaction_coefficients
from projects.vasco_proj.src.Paths import *
import os
import cobra.test
import pandas as pd

#load the model

recontruction=pd.read_excel('C:/Users/vasco/PycharmProjects/human_ts_models/projects/vasco_proj/Data/Transcriptomics/Proc/r3d_compact_mcf7_fastcore_sample_2.xlsx',header=0, index_col=0)
model=cobra.io.read_sbml_model(CMODL_PATH)

recontruction[(recontruction == False).any(axis=1)]
line1=recontruction.loc['fastcore-1.2']

fal_index=list(line1.index[line1==False])

df[(df == 'something1').any(axis=1)]





model.remove_reactions(reactions=fal_index,remove_orphans=True)
model.objective='BIOMASS_maintenance'
biomass_rxn = model.reactions.get_by_id('BIOMASS_maintenance')
model.reactions.get_by_id('BIOMASS_maintenance').bounds
coef_obj=linear_reaction_coefficients(model)

#FBA
solution=model.optimize()
obj_value=solution.objective_value()
fluxes=solution.fluxes
sol_frame=solution.to_frame()
pfba_solution = cobra.flux_analysis.pfba(model)

#FVA

flux_variability_analysis(model,fraction_of_optimum=0.9)