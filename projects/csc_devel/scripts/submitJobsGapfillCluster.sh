#!/usr/bin/env bash
source /home/tbarata/CSCs/human_ts_models/projects/csc_devel/src/paths.py # path to paths.py file. It is hardcoded.
dir=$projFld/$BaseDir/'support/LethalityEval/Reac4Reconst/ReconstOut'
scriptDir=$projFld/$BaseDir/'scripts'
script=$scriptDir/'GapfillCluster.sh'
cd $dir
for stdN in $dir
do cd $stdN
   for fname in *; do qsub -d $projFld -q day -l nodes=1:ppn=24,walltime=24:00:00 -v T=$stdN,FileName=$fname $script; done
done


