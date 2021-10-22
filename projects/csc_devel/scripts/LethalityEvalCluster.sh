#!/bin/bash
source activate py356 # acttivates conda environment
source /home/tbarata/CSCs/human_ts_models/projects/csc_devel/src/paths.py # path to paths.py file. It is hardcoded.
script=$projFld/$BaseDir/'scripts/LethalityEvalCluster.py'
python $script $pos $T
