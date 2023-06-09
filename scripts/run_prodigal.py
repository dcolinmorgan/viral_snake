import pandas as pd
import glob,os
import subprocess

A=pd.read_csv(snakemake.input[0],sep='\t')
vir_pred=A[A.pvalue<.05].name.str.split(' ').str[0]
subprocess.Popen("paste -d ',' - - < final.contigs.fa >final_contigs.fa", shell=True).communicate()
cont=pd.read_csv(snakemake.output[0],sep=',',header=None)
contA=cont[0].str.split(' ').str[0].str.split('>').str[1]
AA=set(contA) & set(vir_pred)
cont['contA']=contA
contB=cont[cont.contA.isin(vir_pred)]
del contB['contA']
contB.to_csv(snakemake.output[0],sep='\n',header=False,index=False)
subprocess.Popen("prodigal -i final_contigs.fna  -o prod.genes -a prod.prots", shell=True).communicate()
