import sys, os

from troppo.tasks.task_io import ExcelTaskIO

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/vvieira/cobamp', '/home/vvieira/human_ts_models', '/home/vvieira/troppo',
                 '/home/vvieira/cobamp/src', '/home/vvieira/troppo/src', '/home/vvieira/human_ts_models'])

from urllib.request import urlretrieve

from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model


ROOT_FOLDER = 'projects/breast_mcf7'
MODEL_PATH = os.path.join(ROOT_FOLDER, 'support/Human-GEM.xml')

def get_human1_model():
    if not os.path.exists(MODEL_PATH):
        path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/modelFiles/xml/HumanGEM.xml')
        model = read_sbml_model(path)
        write_sbml_model(model, MODEL_PATH)
    else:
        model = read_sbml_model(MODEL_PATH)
    model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
    blocked = find_blocked_reactions(model)
    model.remove_reactions(blocked)
    model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
    return model

def get_human1_essential_tasks(model):
    URL = 'https://github.com/SysBioChalmers/Human-GEM/raw/master/data/metabolicTasks/metabolicTasks_Essential.xlsx'
    path, _ = urlretrieve(URL)
    task_list = ExcelTaskIO().read_task(path)

    metab_map = {m.name + '[' + m.compartment + ']': m.id if m.compartment != 'x' else m.id[:-1] + 's' for m in
                 model.metabolites}
    metab_map['NEFA blood pool in[x]'] = 'm02560s'
    replace_func = lambda x: metab_map[x] if x in metab_map.keys() else x

    for t in task_list: t.id_replace(replace_func)  # only works for inflow_dict/outflow_dict
    for t in task_list:
        t.reaction_dict = {k: [{metab_map[i]: j for i, j in v[0].items() if i in metab_map.keys()}, v[1]]
                           for k, v in t.reaction_dict.items()}

    task_list[27].outflow_dict.update(task_list[27].inflow_dict)
    del task_list[27].outflow_dict['ALLMETSIN[s]']
    return task_list
