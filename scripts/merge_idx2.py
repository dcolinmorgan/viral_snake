import glob,os,sys,re
import pandas as pd

# # ## make proper annotations
# # ## based on /groups/cgsd/dcmorgan/virusdb/viral.1.protein.faa
map1=pd.read_csv('~/MAP1.txt',names=['YP','species'],index_col=None,sep='\t')

os.chdir('/home/dcmorgan')

DS=['epilepsy/assemble/*/idxstats.txt',
    '/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/idxstats.txt',
    '/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt',
    '/groups/cgsd/dcmorgan/HT/*/idxstats.txt'
   ]


# R=glob.glob('epilepsy/assemble/*/idxstats.txt')
R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/idxstats.txt')
# R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt')

for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=0,header=None)
    BB=pd.DataFrame(AA[2]+AA[3])
    jeff=j.split('/')[2] # [7] for other projects
    BB.rename(columns={0:jeff},inplace=True)

    idxstat=BB

    idxstat=idxstat.rename_axis('contig').reset_index()#drop=True)

    # idxstat['yp_count']=idxstat.groupby('contig').count()
    # aa=idxstat.sort_values(by='yp_count').drop_duplicates(['contig'],keep='last')
    # print(i,j)

# for i in /groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/viral_hits.blast:
    RR=j.split('idx')[0]
    map2=pd.read_csv(RR+'viral_hits.blast',index_col=None,sep='\t',header=None)
    map2=map2[map2[10]<0.05]
    map2=map2[[0,1,10]]
    map2.columns=['contig','YP','evalue']

# map2=pd.read_csv('~/crc_map3.txt',names=['contig','YP'],index_col=None,sep='\t')
    MAP2=pd.merge(map1,map2)
    MAP2['contig1']=MAP2['contig']
    MAP2['contig']= MAP2.contig.str.split('_').str[0]+'_'+MAP2.contig.str.split('_').str[1]



    idx_out=pd.merge(idxstat,MAP2)
    SS=RR.split('assemble')[1].strip('/')

    jj=idx_out.groupby(['contig','species']).mean()
    jj['count']=idx_out.groupby(['contig','species']).count()['YP'].values
    jj.reset_index(inplace=True)
    jj.evalue=jj.evalue.replace(0,1)
    idx=jj.groupby(['contig'])['count'].transform(max) == jj['count']
    jj=jj[idx]
    jj.evalue=jj.evalue.replace(0,1)
    idx=jj.groupby(['contig'])['evalue'].transform(min) == jj['evalue']
    kk=jj[idx]

    # jj['sum']=jj.drop(columns=['count']).sum(axis=1)
    # kk=jj[~jj.contig.duplicated(keep='first')]
    # kk.drop(columns=['count','evalue']).to_csv("~/epilepsy/annot/"+SS+"_idx_table_M.txt",header=True,index=False,sep='\t')
    
    kk.drop(columns=['count','evalue']).to_csv("/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot/"+SS+"_idx_table_M.txt",header=True,index=False,sep='\t')
    # kk.drop(columns=['count','evalue']).to_csv("/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/"+SS+"_idx_table_M.txt",header=True,index=False,sep='\t')
    del kk, jj, idx
    
    jj=idx_out.groupby(['contig','species']).sum()
    jj['count']=idx_out.groupby(['contig','species']).count()['YP'].values
    jj.reset_index(inplace=True)
    jj.evalue=jj.evalue.replace(0,1)
    idx=jj.groupby(['contig'])['count'].transform(max) == jj['count']
    jj=jj[idx]
    jj.evalue=jj.evalue.replace(0,1)
    idx=jj.groupby(['contig'])['evalue'].transform(min) == jj['evalue']
    kk=jj[idx]
    
    # kk.drop(columns=['count','evalue']).to_csv("~/epilepsy/annot/"+SS+"_idx_table_S.txt",header=True,index=False,sep='\t')
    
    kk.drop(columns=['count','evalue']).to_csv("/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot/"+SS+"_idx_table_S.txt",header=True,index=False,sep='\t')
    # kk.drop(columns=['count','evalue']).to_csv("/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/"+SS+"_idx_table_S.txt",header=True,index=False,sep='\t')

# R=glob.glob('epilepsy/annot/*idx_table_M.txt')
R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot/*_idx_table_M.txt')
# R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/*idx_table_M.txt')

for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=None)
    SS=j.split('annot')[1].split('_')[0].strip('/')
    AA.columns=['contig','species',SS]
    if i!=0:
        BB=AA.merge(BB,on=['contig','species'],how='outer')#,left_on=['contig','species'])
    else:
        BB=AA
# BB.fillna(0).to_csv('epilepsy/annot_contigs_M.txt',sep='\t')
BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_contigs_M.txt',sep='\t')
# BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_contigs_M.txt',sep='\t')

# R=glob.glob('epilepsy/annot/*idx_table_S.txt')
R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot/*_idx_table_M.txt')
# R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/*idx_table_M.txt')
for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=None)
    SS=j.split('annot')[1].split('_')[0].strip('/')
    AA.columns=['contig','species',SS]
    if i!=0:
        BB=AA.merge(BB,on=['contig','species'],how='outer')#,left_on=['contig','species'])
    else:
        BB=AA
# BB.fillna(0).to_csv('epilepsy/annot_contigs_S.txt',sep='\t')
BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_contigs_M.txt',sep='\t')
# BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_contigs_M.txt',sep='\t')


# data=pd.read_csv('epilepsy/annot_contigs.txt',sep='\t')
# # data=pd.read_csv('contig_counts.tsv',sep='\t')
# data=data.melt(['species','contig']).replace({'-DNA':' '}, regex=True)
# data.variable.replace(' ','',regex=True,inplace=True)
# label=pd.read_csv('epilepsy/assemble/label.txt',sep='\t')
# data=pd.merge(data,label,left_on='variable',right_on='id')
# data.rename(columns={'value':'abundance'},inplace=True)

# data1=data[data.label=='Non-refractory']
# sub_gen=data1.groupby('species').count().sort_values('contig',ascending=False).reset_index().species[1:20]
# data1=data1[data1['species'].isin(sub_gen)]

# data2=data[data.label=='Refractory']
# sub_gen=data2.groupby('species').count().sort_values('contig',ascending=False).reset_index().species[1:20]
# data2=data2[data2['species'].isin(sub_gen)]

# plt.figure(figsize=(12,8))
# ax = sns.histplot(data1, x='variable', hue='species', weights='abundance',
#              multiple='stack', palette='tab20c', shrink=0.8)
# ax.set_ylabel('log_abundance')
# ax.set_xlabel('non-refractory')
# ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

# # Fix the legend so it's not on top of the bars.
# legend = ax.get_legend()
# legend.set_bbox_to_anchor((1, 1))
# plt.savefig('non-refractory.png')

# plt.figure(figsize=(12,8))
# ax = sns.histplot(data2, x='variable', hue='species', weights='abundance',
#              multiple='stack', palette='tab20c', shrink=0.8)
# ax.set_ylabel('log_abundance')
# ax.set_xlabel('refractory')
# ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

# # Fix the legend so it's not on top of the bars.
# legend = ax.get_legend()
# legend.set_bbox_to_anchor((1, 1))
# plt.savefig('refractory.png')