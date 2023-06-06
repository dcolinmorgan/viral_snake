# esearch -db sra -query PRJNA544527
# esearch -db sra -query PRJNA544527 | esummary | xtract -pattern DocumentSummary -element Biosample,Sample@acc,Experiment@acc,Run@acc,Platform@instrument_model > PRJNA544527-info.tsv
# head -n 2 PRJNA544527-info.tsv > PRJNA544527-info-subset.tsv

ml sratoolkit/3.0.0

for accession in $(cut -f 2 PRJNA544527-metagen.tsv)
do

    printf "\n  Working on: ${accession}\n\n"
    fasterq-dump --split-files ${accession}

done