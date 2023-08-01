
import pandas as pd
import numpy as np
# import chart_studio.plotly as cspy
import plotly.graph_objects as go
import plotly.express as px
import seaborn as se
# from google.colab import files
# print("hello world")
import glob,os,sys,re
import matplotlib.pyplot as plt
import seaborn as sns, numpy as np, pandas as pd
from numpy import inf


data=pd.read_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt',sep='\t',usecols=range(1,540))
data.fillna(0,inplace=True)
# data
# !wget https://raw.githubusercontent.com/dcolinmorgan/grph/main/PRJNA544527-meta_inf.txt

meta=pd.read_csv('~/snakemake_viral_calling/PRJNA544527-meta_inf.txt',sep='\t',header=None)
meta=meta.rename(columns={3:'id',5:'pat'})

data.index=data.species
data=data.T

mm=pd.merge(data,meta[['id','pat']],left_index=True,right_on='id')
# mm=mm.T

# mm['contig']=data.contig
# mm=mm.T
# mm.species=data.T.index
mm['id']=mm['pat'].str.split('-').str[0]
mm['time']=mm['pat'].str.split('_').str[0].str.split('-').str[1]

# !wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0559-3/MediaObjects/41591_2019_559_MOESM3_ESM.xlsx --no-check-certificate
metaa=pd.read_excel('~/snakemake_viral_calling/41591_2019_559_MOESM3_ESM.xlsx',sheet_name='SupTable2',skiprows=3)
metaa=metaa[['Donor','Age','Sex','BMI']]


Full_table=pd.merge(mm,metaa,left_on='id',right_on='Donor')
# Full_table=Full_table.drop(columns=[3,	'pat',	'id'])
# Full_table.Label=pd.to_dateLabel(Full_table.Label,unit='d')
# Full_table.Label=Full_table.Label.values.astype('dateLabel64[M]')

data2=Full_table.melt(id_vars=['time','Donor','Age','Sex','BMI','id','pat'])#,'species','contig'])
# data2=data2.rename(columns={'variable':'species'})
data2=data2.sort_values(by=['Donor','time','value'])
# data2=data2[data2.value>0]
# data2.to_csv('PRJNA544527_mpa4_annot_table.txt',sep='\t')

# final df stored here also
# !wget https://raw.githubusercontent.com/dcolinmorgan/grph/main/PRJNA544527_mpa4_annot_table.txt
# data2=pd.read_csv('PRJNA544527_mpa4_annot_table.txt',sep='\t',index_col=0)

data2=data2[data2.value>0]
data2=data2.reset_index(drop = True)
# data2=data2.groupby(['Donor','Label','value'],group_keys=False).apply(lambda x : x.sort_values(by = 'value', ascending = False).head(10))
data2=data2.drop_duplicates()

data2["Label"] = (
    data2.groupby("Donor")
    .apply(lambda x: x.groupby("time", sort=False).ngroup() + 1)
    .values
)


cc=pd.unique(data2[data2.Label<5].Donor)
# len(pd.unique(data2.Donor))
data2=data2.loc[ data2.Donor.isin(cc), : ]
data2=data2[data2.Label<5]
data2=data2[data2['rank']<10.0]
data2.to_csv('~/PRJNA544527_viral10top.txt',sep='\t')