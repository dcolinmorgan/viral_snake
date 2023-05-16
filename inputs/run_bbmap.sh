#!/bin/sh

#PBS -l nodes=1:ppn=48
#PBS -l mem=100g
#PBS -l walltime=96:00:00
#PBS -N bbm
#PBS -q cgsd 
#PBS -e .bberr
#PBS -o .bbout
###qsub run/epilepsy/run_bbmap.sh

ml miniconda3/4.12.0
module add parallel
# ml bedtools
ml bowtie2
# source activate dvf
source activate mypy3


# index contigs
# > map raw reads (before megahit) to contigs -- samtools/ bedtools
# >> match with blast output


# cd /groups/cgsd/shengxu/epilepsy/LauG_Metagenomics_CPOS-221215-MHWK-15822a/LauG_Metagenomics_CPOS-221215-MHWK-15822a/clean_data/
# https://github.com/edamame-course/Metagenome/blob/master/2018-06-29-counting-abundance-with-mapped-reads.md

function m_bbmap {
# for i in *
# do
i=$(eval "echo "$1" | cut -d / -f9 |cut -d _ -f1")
j=$(eval "echo "$1" | cut -d _ -f1,2,3,4,5,6")
tr ',' '\n' < ~/assemble/${i}/blast_contigs.fa > ~/assemble/${i}/blast_contigs.faa
bowtie2-build ~/assemble/${i}/blast_contigs.faa ~/assemble/${i}/build_v
bowtie2 -x ~/assemble/${i}/build_v -1 ${j}_clean_1.fastq -2 ${j}_clean_2.fastq -S ~/assemble/${i}/eg2.sam


samtools view -bS ~/assemble/${i}/eg2.sam > ~/assemble/${i}/aln.bam
samtools sort ~/assemble/${i}/aln.bam -o ~/assemble/${i}/sorted_aln.bam
# samtools index ~/assemble/${i}/sorted_aln.bam
# samtools view -c -F 0 ~/assemble/${i}/sorted_aln.bam
samtools idxstats ~/assemble/${i}/sorted_aln.bam > ~/assemble/${i}/idxstats.txt
# git clone https://github.com/metajinomics/mapping_tools.git
# samtools depth  ~/assemble/${i}/aln.bam  > ~/assemble/${i}/depth.txt #|  awk '{sum+=$3} END { print "Average = ",sum/NR}'



}
export -f m_bbmap

parallel m_bbmap ::: /groups/cgsd/shengxu/epilepsy/LauG_Metagenomics_CPOS-221215-MHWK-15822a/LauG_Metagenomics_CPOS-221215-MHWK-15822a/clean_data/*



for subdir in *; do cp $subdir/idxstats.txt $subdir.idxstats.txt; done;
find . -size 0 -delete
python2 ~/run/epilepsy/mapping_tools/get_count_table.py *.idxstats.txt > contig_counts.tsv

# must need to index blast_contifs.fa file ## could use bowtie instead

# i=$(eval "echo "$1" | cut -d / -f9 |cut -d _ -f1")
# tr ',' '\n' < ~/assemble/GMOREA0160-DNA/blast_contigs.fa > ~/assemble/GMOREA0160-DNA/blast_contigs.faa
# bowtie2-build ~/assemble/GMOREA0160-DNA/blast_contigs.faa ~/assemble/GMOREA0160-DNA/build_v
# bowtie2 -x ~/assemble/GMOREA0168-DNA/build_v -1 ${1}_clean_1.fastq -2 ${1}_clean_2.fastq -S ~/assemble/${i}/eg2.sam
# samtools view -bS ~/assemble/${i}/eg2.sam > ~/assemble/${i}/aln.bam
# samtools depth  ~/assemble/${i}/aln.bam  > ~/assemble/${i}/depth.txt