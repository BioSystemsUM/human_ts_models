#!/usr/bin/env bash
source /home/tbarata/CSCs/human_ts_models/projects/csc_devel/src/paths.py # path to paths.py file. It is hardcoded.
dir=$projFld/$BaseDir/'support/Reconst/ReconstOut'
scriptDir=$projFld/$BaseDir/'scripts'
script=$scriptDir/'StdGapfill.sh'
cd $dir
for fname in *; do qsub -d $projFld -q day -l nodes=1:ppn=24,walltime=24:00:00 -v FileName=$fname $script; done


