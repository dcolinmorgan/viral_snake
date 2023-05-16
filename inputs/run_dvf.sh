#!/bin/sh

#PBS -l nodes=1:ppn=48
#PBS -l mem=500g
#PBS -l walltime=96:00:00
#PBS -N dvf
#PBS -q cgsd 
#PBS -e .dvfe
#PBS -o .dvfo
#qsub run/epilepsy/dvf.sh

# module load python3 cuda11.0/toolkit/11.0.3 cudnn8.0-cuda11.0/8.0.5.39

# module purge
ml miniconda3/4.12.0
module add parallel
source activate dvf
# source activate mypy3
# source activate myGPU


## assemble and run metaphlan
## qsub run/epilepsy/pbs_assemble
## qsub run/gcn/pbs_metaphlan


## get deepvirfinder scores and use that to filter input into prodigal
# assembledir='assemble2'
# contigs=$(ls assemble/*/*.fa)

function dvfcont {
# for i in $contigs
# do
i=$1
# python ~/run/virus/DeepVirFinder/dvf.py -i $i -l 200 -c 12
python ~/run/virus/DeepVirFinder/dvf.py -i $i -l 500 -c 12
# python ~/run/virus/DeepVirFinder/dvf.py -i $i -l 1000 -c 12

# done
}
export -f dvfcont

parallel dvfcont ::: assemble/*/final.contigs.fa


## dont really need genes since proteins will blast agains NCBI viral DB
# prodigal -i assemble/GMOREA0142-DNA/final_contigs.fna  -o assemble/GMOREA0142-DNA/prod.genes -a assemble/GMOREA0142-DNA/prod.prots  -p anon 
#### run this below ####
# python run/epilepsy/run_prodigal.py 

# blast NCBI refseq  /groups/cgsd/hbyao/arg/ncbi/2020-06-11.1  ncbiamr-06-11.1.zip  pbs_hmm