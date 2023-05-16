#!/bin/sh

#PBS -l nodes=1:ppn=48
#PBS -l mem=100g
#PBS -l walltime=96:00:00
#PBS -N bpz
#PBS -q cgsd 
#PBS -e .bpze
#PBS -o .bpzo
#qsub run/epilepsy/run_blastp.sh

# module load python3 cuda11.0/toolkit/11.0.3 cudnn8.0-cuda11.0/8.0.5.39

# module purge
ml miniconda3/4.12.0
module add parallel
# source activate dvf
source activate mypy3
# source activate myGPU




contigs=$(ls assemble/*/prod.prots)

function blastpN {
# for i in $contigs
# do
# blastp $i 
# i=$1
i=$(eval "echo "$1" | cut -d / -f1,2")

blastp -query $1 -db /groups/cgsd/dcmorgan/virusdb/refseq_viral_db -out $i/viral_hits.blast -outfmt 6
# done
}
export -f blastpN

parallel blastpN ::: assemble/*/prod.prots

