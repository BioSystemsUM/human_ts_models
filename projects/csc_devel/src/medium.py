### Modules:
import os
import pandas as pd

### Functions:
def CloseInfluxAllEx(model):
    '''
    - close influx of all exchange reactions
    :param model: generic metabolic model
    :return: model with all schange reactions within (0,1000)
    '''
    for r in model.exchanges:
            r.lower_bound = 0
    return model

def AddMediumBounds(MediumBoundsPath, model):
    '''
    - Adds medium bounds to exchange reactions associated with experimental medium compounds
    :param MediumBoundsPath: path to file containing bounds of medium exchange reactions
    :param model: generic metabolic model
    :return: model with flux constrains on exchange reactions of a medium
    '''
    MediumBoundsTable = pd.read_excel(MediumBoundsPath)
    MediumBoundsDic = MediumBoundsTable.set_index('Reaction_ID').T.to_dict()
    for k, v in MediumBoundsDic.items():
        LB = float(v['Lower_bound'])
        UB = float(v['Upper_bound'])
        model.reactions.get_by_id(k).bounds = (LB, UB)
    return model

def AdjustExReac (model, MediumBoundsPath):
    '''
    - adjusts exchange reactions for medium composition: influx of all exchange reactions is closed except for medium compounds
    :param model: generic metabolic model
    :param MediumBoundsPath: path to file containing bounds of medium exchange reactions
    :return: model with flux constrains on exchange reactions of a medium and remaining exchange reactions closed (0,1000)
    '''
    modelAdj = model.copy()
    CloseInfluxAllEx(modelAdj)
    AddMediumBounds(MediumBoundsPath, modelAdj)
    return modelAdj

def findExtraMediumMetabolites(model, HamsPath, BaseMediumPath, MediumPath):
    '''
    - when a model cannot grow in a medium, this funct identifies metabolites (of HAMs medium) that have to be added to first medium in order for model to grow
    :param model: generic metabolic model
    :param HamsPath: path to file with HAMs medium composition
    :param BaseMediumPath: path to file with composition of medium we want to test
    :param MediumPath: path to file where results will be saved
    :return: creates file with original medium composition plus metabolites of HAMs medium that are needed for model to grow
    '''
    if not os.path.exists(MediumPath):
        hamsTable = pd.read_excel(HamsPath)
        MediumTable = pd.read_excel(BaseMediumPath)
        extra = set(hamsTable.Reaction_ID) - set(MediumTable.Reaction_ID) # metabolites present in Ham's medium but not on other medium
        # add basic medium bounds plus extra metabolites from Ham's medium:
        testModel = model.copy()
        testModel = AdjustExReac(testModel, BaseMediumPath)
        for rid in extra:
            testModel.reactions.get_by_id(rid).bounds = (float(-1000), float(1000))
        # find 'extra' metabolites that are not essential and can be left out of medium formulation:
        extraFinal = list()
        for rid in extra:
            testModel.reactions.get_by_id(rid).bounds = (float(0),float(1000))
            opt = testModel.optimize()
            if opt.status != 'optimal' or opt.objective_value <= 0:
                extraFinal.append(rid)
            testModel.reactions.get_by_id(rid).bounds = (float(-1000), float(1000))
            # note: we should manually check if metabolites of 'extraFinal' reactions are biologically acceptable
        for rid in extraFinal:
            extraR = hamsTable.loc[hamsTable.Reaction_ID == rid]
            MediumTable = MediumTable.append(extraR, ignore_index=True, sort=False)
        writer = pd.ExcelWriter(MediumPath, engine='xlsxwriter')
        MediumTable.to_excel(writer, sheet_name='Bounds', index=False)
        writer.save()