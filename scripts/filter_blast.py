import pandas as pd
import glob,os
import subprocess

reads=pd.read_csv(snakemake.input[0],names=['name','seq'])
hits=pd.read_csv('viral_hits.blast',sep='\t',names=['name','tmp','identity','tmp1','tmp2','tmp3','tmp4','tmp5','tmp8','tmp12','tmp13','tmp11'])
hits=hits[['name','identity']]
AA=set(reads.name.str.split(' ').str[0]) & set('>'+hits.name.str.split('_').str[0]+'_'+hits.name.str.split('_').str[1])
reads[reads['name'].str.split(' ').str[0].isin(AA)].to_csv('snakemake.output[0]',header=False,index=False)