
#PBS -l nodes=1:ppn=40
#PBS -l mem=200gb
#PBS -N mpa4_2
#PBS -l walltime=444:00:00
#PBS -q cgsd
#PBS -e .mpa4e
#PBS -o .mpa4o
###qsub run/epilepsy/run_mpa.sh


##cd $PBS_O_WORKDIR

##cd /groups/cgsd/shengxu/CRC_data
##cd /groups/cgsd/shengxu/ARG_data/fastq/ ##no_duplicated

module load miniconda3
source activate mypy3
ml gcc
ml python3
ml bowtie2

# cd /groups/cgsd/shengxu/epilepsy/LauG_Metagenomics_CPOS-221215-MHWK-15822a/LauG_Metagenomics_CPOS-221215-MHWK-15822a/
# for i in $(cat idlist)
# do
#     metaphlan no_duplicated/${i}_good_out_R1.fastq,no_duplicated/${i}_good_out_R2.fastq --bowtie2out ~/epi_metaphlan_out/${i}.bowtie2.bz2 --nproc 40 --input_type fastq -o m_out/${i}.txt -t rel_ab_w_read_stats -s ~/epi_sam_out/${i}.sam
# done


cd /groups/cgsd/shengxu/epilepsy/LauG_Metagenomics_CPOS-221215-MHWK-15822a/LauG_Metagenomics_CPOS-221215-MHWK-15822a/old-data/

for i in $(cat idlist_2)
do
    metaphlan no_duplicated/${i}_good_out_R1.fastq,no_duplicated/${i}_good_out_R2.fastq --bowtie2out ~/epi_metaphlan_out/${i}.bowtie2.bz2 --nproc 40 --input_type fastq -o ~/assemble/${i}/mpa4_out.txt -t rel_ab_w_read_stats -s ~/epi_sam_out/${i}.sam
done