#!/usr/bin/env bash
source /home/tbarata/CSCs/human_ts_models/projects/csc_devel/src/paths.py # path to paths.py file. It is hardcoded.
sizeFile=$projFld/$BaseDir/'support/Reconst/NumbModels.tab'
scriptDir=$projFld/$BaseDir/'scripts'
script=$scriptDir/'ReconstStdCluster.sh'
size=$(cat $sizeFile)
seqs=$(seq 0 $(echo $(( $size - 1 ))))
for i in $seqs; do qsub -d $projFld -q day -l nodes=1:ppn=24,walltime=24:00:00 -v pos=$i $script; done








