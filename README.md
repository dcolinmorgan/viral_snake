# viral_snake : a parallelized, system-agnostic snakemake pipeline (and singularity image) for viral calling of metagenomic data

>__Conclusions__: Parallelizing a proven viral calling pipeline and systematizing with snakemake makes for faster, reproducible metagenomic research....<br>
>__Authors__: Daniel Morgan; Haobin Yao; Joshua Ho

<space>\
<space>

Clone current version & run [see build file](https://cloud.sylabs.io/builder/6456f006a9850b03b2ca613e) and [here](https://github.com/dcolinmorgan/viral_snake/blob/master/tryBuild.def)
--------------------------------------------------
```bash
    ml singularity
    singularity pull --arch amd64 library://dcolinmorgan/vc/viral_calling:v0.1
    ml snakemake
    cd  <working dir>
    snakemake --profile profile/ --jobs 8
```

This analysis was performed on the [BIO-ML](https://www.broadinstitute.org/infectious-disease-and-microbiome/broad-institute-openbiome-microbiome-library) metagenomics, time-series dataset, from the [2019 Nature Medicine paper (ref. below)](https://sci-hub.se/10.1038/s41591-019-0559-3). Data can be found in serveral places, and collecting and coallating was messy, so I've included files to download the data in this repo.
```Poyet, M., Groussin, M., Gibbons, S.M. et al. A library of human gut bacterial isolates paired with longitudinal multiomics data enables mechanistic microbiome research. Nat Med 25, 1442–1452 (2019).```

*   [Metadata](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0559-3/MediaObjects/41591_2019_559_MOESM3_ESM.xlsx)
*   [Raw Fastq Data](https://www.ebi.ac.uk/ena/browser/view/PRJNA544527), in a proper format for [ftp download here](https://github.com/dcolinmorgan/grph/blob/main/ftp_PRJNA544527.txt) via this [helper file](https://github.com/dcolinmorgan/grph/blob/main/PRJNA544527-meta_inf.txt)
```bash
    ml parallel
    function ftpall {
    wget ftp://$1
    }
    export -f ftpall
    parallel -j 20 ftpall :::: ftp_PRJNA544527.txt
```

IMPORTANT: singularity image gives user installed versions of required packages, also attained in any pyenv, as follows:
*   conda: bioconda pandas megahit numpy prodigal bowtie bbmap hmmer
*   conda3.6: dvf python=3.6 numpy theano=1.0.3 keras=2.2.4 scikit-learn Biopython h5py
*   pip: MetaPhlan
*   github: DeepVirFinder, viralrecall
  
<space>\
<space>
  
raw data processing, viral contigs identification, and viral taxonomy annotation
---------------------
1. Assemble contigs with [megahit](https://github.com/voutcn/megahit)
    1. intput: raw paired end fastq.gz files from ftp link ({sample}_1.fastq.gz,{sample}_2.fastq.gz)
    2. output: NGS assembled contigs ({sample}.fa)
2. Identify bacterial contigs with [MetaPhlan4](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4) -- confirm data quality in comparison to BIO-ML publication, remove from downstream analysis
    1. intput: raw paired end fastq.gz files from ftp link ({sample}_1.fastq.gz,{sample}_2.fastq.gz)
    2. output: metaphlan table ({sample}.txt)
2. Identify viral contigs with [DeepVirFinder](https://github.com/jessieren/DeepVirFinder)
    1. intput: assembled contigs (({sample}.fa))
    2. output: viral contig predictions (final.contigs.fa_gt{params}bp_dvfpred.txt, where params is bp cutoff)
2. Predict viral genes with [Progical](https://github.com/hyattpd/Prodigal), run script [here](https://github.com/dcolinmorgan/viral_snake/blob/master/scripts/run_prodigal.py)
    1. intput: viral contig predictions (final.contigs.fa_gt{params}bp_dvfpred.txt, where params is bp cutoff)
    2. output: viral genes from predicted contigs (final_contigs.fna) 
2. ensure contigs are viral via blastp and filtering against [viral refseq database](https://support.nlm.nih.gov/knowledgebase/article/KA-03474/en-us)
    1. intput: viral genes from predicted contigs (final_contigs.fna) + refseq viral database
    2. output: (viral_hits.blast) & filter down to (blast_contigs.fa)
2. bowtie and samtools to produce feature table per sample, merge into single table
    1. intput: (blast_contigs.fa)
    2. output: (idxstats.txt) per sample & merge
3. Remove bacterial contigs and identify viral taxonomy with [viralrecall](https://github.com/faylward/viralrecall)
    1. intput: format then (cat_blast_contigs.faa)
    2. output: {sample}/VR_out and merge to get count table (contig_counts.tsv)
5. [plot metaplhlan abundance](https://github.com/dcolinmorgan/viral_snake/blob/master/scripts/plot_abundance.py)
   
Abundance, UMAP& Time-Series analysis
---------------------
requirements :
cudf cuml graphistry cu-cat seaborn

```bash
    pip install --extra-index-url=https://pypi.nvidia.com cuml-cu11 cudf-cu11 cugraph-cu11 pylibraft_cu11 raft_dask_cu11 dask_cudf_cu11 pylibcugraph_cu11 pylibraft_cu11
    pip install -U --force git+https://github.com/graphistry/pygraphistry.git@cudf
    pip install -U git+https://github.com/graphistry/cu-cat.git@DT3
    nvidia-smi
```
    
Following this [jupyter notebook]([https://github.com/dcolinmorgan/viral_snake/....](https://github.com/dcolinmorgan/grph/blob/main/accelearting_metagenomic_demo.ipynb)) for species abundance, umap, time-series analysis
    
Among other things, these checks are performed herewithin:

>Workflow figure from Nature Comms paper [previous manuscript](https://www.nature.com/articles/s41467-023-37975-y#Sec17)
>--------------------------------------------------


>![Figure 1. Pipeline: raw data processing, viral contigs identification, and viral taxonomy annotation](https://github.com/dcolinmorgan/viral_snake/blob/master/viral_snake_partial_pipe.png)\
> __Figure 1. The workflow of the raw data processing, viral contigs identification, and viral taxonomy annotation. The main workflow of this study was composed of three parts: raw data preprocessing, viral contigs identification, and viral taxonomy annotation. Viral contigs identification involves viral contigs identification and 2 rounds bacterial genome removal. The viral taxonomy
annotation also includes taxonomy annotation and NCLDV contigs verification.
```Wang, L., Yao, H., Morgan, D.C. et al. Altered human gut virome in patients undergoing antibiotics therapy for Helicobacter pylori. Nat Commun 14, 2196 (2023). https://doi.org/10.1038/s41467-023-37975-y```
