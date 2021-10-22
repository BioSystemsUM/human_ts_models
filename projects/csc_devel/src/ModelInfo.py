### Modules:
import os
import pandas as pd
import cobra
from collections import OrderedDict
import re
from cobra.io import read_sbml_model
from medium import AdjustExReac

### Functions:
def createiHumanInfo (model, modelInfoPath):
    '''
    - Creates a excel file with info on iHuman model with a tab for reactions, a tab for metabolites and a tab for genes
    :param model: basic recon metabolic model
    :param modelInfoPath: path to file where to save model info
    :return: saves model info into a file
    '''
    # create dataframe with iHuman model info on reactions:
    reactions = pd.DataFrame(OrderedDict((('ID', pd.Series([r.id for r in model.reactions])),
                                          ('NAME', pd.Series([r.name for r in model.reactions])),
                                          ('REACTION', pd.Series([r.reaction for r in model.reactions])),
                                          ('COMPARTMENT', pd.Series([r.compartments for r in model.reactions])),
                                          ('SUBSYSTEM', pd.Series([r.subsystem for r in model.reactions])),
                                          ('LOWER_BOUND', pd.Series([r.lower_bound for r in model.reactions])),
                                          ('UPPER_BOUND', pd.Series([r.upper_bound for r in model.reactions])),
                                          ('GENE_REACTION_RULE', pd.Series([r.gene_reaction_rule for r in model.reactions])),
                                          ('EC-CODE', pd.Series([v for r in model.reactions for k,v in r.annotation.items() if k=='ec-code'])),
                                          ('SBO', pd.Series([v for r in model.reactions for k,v in r.annotation.items() if k=='sbo'])),
                                          ('BIGG.REACTION', pd.Series([v for r in model.reactions for k,v in r.annotation.items() if k=='bigg.reaction'])),
                                          ('KEGG.REACTION', pd.Series([v for r in model.reactions for k,v in r.annotation.items() if k=='kegg.reaction'])),
                                          ('METANETX.REACTION', pd.Series([v for r in model.reactions for k,v in r.annotation.items() if k=='metanetx.reaction']))
                                           )))
    reactions.COMPARTMENT = reactions.COMPARTMENT.apply(lambda x: re.sub("}", ")", re.sub("{", "(", str(x)))) # formats column 'COMPARTMENT'
    # create dataframe with recon iHuman model info on metabolites:
    metabolites = pd.DataFrame(OrderedDict((('ID', pd.Series([m.id for m in model.metabolites])),
                                            ('NAME', pd.Series([m.name for m in model.metabolites])),
                                            ('FORMULA', pd.Series([m.formula for m in model.metabolites])),
                                            ('COMPARTMENT', pd.Series([m.compartment for m in model.metabolites]))
                                            )))
    # create dataframe with iHuman model info on genes:
    genes = pd.DataFrame(OrderedDict((('ID', pd.Series([g.id for g in model.genes])),
                                      ('REACTIONS', pd.Series([g.reactions for g in model.genes]))
                                       )))
    writer = pd.ExcelWriter(modelInfoPath, engine='xlsxwriter')
    reactions.to_excel(writer, sheet_name='REACTIONS', index=False) # create tab with reactions info
                                                                    # index=False means don't add row names
    metabolites.to_excel(writer, sheet_name='METABOLITES', index=False) # create tab with metabolites info
    genes.to_excel(writer, sheet_name='GENES', index=False) # create tab with genes info
    writer.save()

def preprocessiHuman(modelPathSBML, modelProcPath, modelProcInfoPath, modelMedAdapPath, modelInfoPath, MediumPath, BiomassID):
    '''
    - saves info (reactions, genes, metabolites) on generic model, pre-processed model and model adapted for medium composition on .xlsx files
    - removes blocked reactions, except those that are exchange reactions related with medium metabolites
    - gives model adapted for medium composition
    :param modelPathSBML: link to where generic iHuman(human1) model is
    :param modelProcPath: path to pre-processed metabolic model
    :param modelMedAdapPath: path to metabolic model adapted for medium composition
    :param modelInfoPath: path to file with info on generic model
    :param MediumPath: path to output file with medium composition to be tested
    :param BiomassID: id of biomass reaction
    :return: returns generic model and model adapted for medium composition; and saves pre-processed/medium adapted/generic model info
    '''
    if not os.path.exists(modelMedAdapPath):
        model = read_sbml_model(modelPathSBML) # load model in SBML format cause .mat format does not have boundaries
        model.objective = BiomassID # define biomass production as objective
        createiHumanInfo(model, modelInfoPath)
        MedReac = set(pd.read_excel(MediumPath).Reaction_ID)  # boundaries for medium components
        blocked_reactions = cobra.flux_analysis.variability.find_blocked_reactions(model) # identify blocked reactions
        Toremove = list(set(blocked_reactions) - MedReac)
        model.remove_reactions(Toremove, remove_orphans=True) # remove blocked reactions except boundaries for medium components
        createiHumanInfo(model, modelProcInfoPath)  # save info on processed model to a xlsx
        cobra.io.write_sbml_model(model, modelProcPath)
        # close input drains except those of medium - which have specific bounds:
        modelMedAdap = model.copy()
        modelMedAdap = AdjustExReac(modelMedAdap, MediumPath)
        if modelMedAdap.optimize().status != 'optimal' or modelMedAdap.optimize().objective_value <= 0:
            print('No biomass is produced or model infeasible. Change medium composition')
        else:
            cobra.io.write_sbml_model(modelMedAdap, modelMedAdapPath)
    else:
        model = cobra.io.read_sbml_model(modelProcPath)
        modelMedAdap = cobra.io.read_sbml_model(modelMedAdapPath)
    return model, modelMedAdap

