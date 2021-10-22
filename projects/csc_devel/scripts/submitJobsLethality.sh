#!/usr/bin/env bash
source /home/tbarata/CSCs/human_ts_models/projects/csc_devel/src/paths.py # path to paths.py file. It is hardcoded.
scriptDir=$projFld/$BaseDir/'scripts'
script=$scriptDir/'LethalityEvalCluster.sh'
dir=$projFld/$BaseDir/'support/TestGeneScores/models2Test'
cd $dir
for std in $(ls | grep .*.pkl)
do File=$(echo $std | rev | cut -d'.' -f2- | rev)
   sizeFile=$projFld/$BaseDir/'support/TestGeneScores/models2Test/NumbModels'$File'.tab'
   size=$(cat $sizeFile)
   seqs=$(seq 0 $(echo $(( $size - 1 ))))
   for i in $seqs; do qsub -d $projFld -q day -l nodes=1:ppn=24,walltime=24:00:00 -v pos=$i,T=$File $script; done
done







