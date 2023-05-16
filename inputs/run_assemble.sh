#PBS -l nodes=1:ppn=24
#PBS -l mem=40gb
#PBS -N asmbl
#PBS -l walltime=24:00:00
#PBS -q cgsd
#PBS -e .asmble
#PBS -o .asmblo
###qsub run/epilepsy/run_assemble.sh

# cd /groups/cgsd/shengxu/epilepsy/LauG_Metagenomics_CPOS-221215-MHWK-15822a/LauG_Metagenomics_CPOS-221215-MHWK-15822a/

# module load megahit
# #for sample in `cat $idlist`
# for i in $(cat idlist_2)
# do
#     megahit -1 clean_data/${i}_clean_1.fastq -2 clean_data/${i}_clean_2.fastq -t 48 -o ~/assemble/$i
# done



# cd /groups/cgsd/shengxu/epilepsy/LauG_Metagenomics_CPOS-221215-MHWK-15822a/LauG_Metagenomics_CPOS-221215-MHWK-15822a/old-data/

# module load megahit
# #for sample in `cat $idlist`
# for i in $(cat idlist)
# do
#     megahit -1 clean_data/${i}_clean_1.fastq -2 clean_data/${i}_clean_2.fastq -t 48 -o ~/assemble/$i
# done



# seq=$(ls /groups/cgsd/shengxu/CRC_data/no_duplicated/*good_out_R1*)
seq=$(ls /groups/cgsd/shengxu/CRC_data/no_duplicated/*good_out_R1* | cut -d / -f7| cut -d _ -f1)

module load megahit
#for sample in `cat $idlist`
for i in $(eval "echo "$seq" | cut -d _ -f1")

do
    
    megahit -1 /groups/cgsd/shengxu/CRC_data/no_duplicated/${i}_good_out_R1.fastq -2 /groups/cgsd/shengxu/CRC_data/no_duplicated/${i}_good_out_R2.fastq -t 48 -o /groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/${i}
done