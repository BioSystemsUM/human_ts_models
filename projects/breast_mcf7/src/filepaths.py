import os

# paths to necessary files
# this is still hardcoded
ROOT_FOLDER = 'projects/breast_mcf7'
SHR_FOLDER = 'shared/'

MODEL_PATH = os.path.join(SHR_FOLDER, 'models/Recon3D/Recon3D.xml')  # SBML model path
DATA_PATH = os.path.join(ROOT_FOLDER,
                         'data/ccle/DepMap Public 20Q1/CCLE_expression_log2pp.csv')  # tabular csv with expression
TASKS_PATH = os.path.join(SHR_FOLDER, 'task_sets/nl2019_tasks_r3d_compact.json')  # JSON containing metabolic tasks
SUBTASKS_PATH = os.path.join(ROOT_FOLDER, 'support/nl2019_tasks_r3d_bigg_compact.json')  # JSON containing metabolic tasks
CMODL_PATH = os.path.join(ROOT_FOLDER,
                          'support/Recon3D_bigg_consistent.xml')  # Consistent SBML model path - recommended!
CMODL_MEDIA_PATH = os.path.join(ROOT_FOLDER,
                                'support/Recon3D_bigg_MEM_consistent.xml')  # Consistent SBML model path - recommended!

TASK_RESULTS_PATH = os.path.join(ROOT_FOLDER, 'results',
                                 'r3d_compact_task_results_mcf7_minsum.json')  # task evaluation
CS_MODEL_DF_PATH = os.path.join(ROOT_FOLDER, 'results',
                                'r3d_compact_mcf7_fastcore_minsum.csv')  # context-specific models extracted from algos
MEDIA_DICT_PATH = os.path.join(ROOT_FOLDER, 'support/r3d_media_metabolites.json')
