import glob,os,sys,re
import pandas as pd

# ## make proper annotations
# ## based on /groups/cgsd/dcmorgan/virusdb/viral.1.protein.faa
map1=pd.read_csv('~/map.txt',names=['YP','species','orf'],index_col=None,sep='[')

map1.YP=map1.YP.str.split('.1').str[0]
map1.YP=map1.YP.str.replace('>','')
map1.species=map1.species.str.replace(']','')
map1.YP=map1.YP+'.1'
map1.to_csv('~/MAP1.txt',sep='\t',index=False)

## based on cut -f1,2 /groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/viral_hits.blast > ~/map_CRC.txt
# sort  map_CRC.txt | uniq -u > crc_map3.txt

## based on cut -f1,2 /groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/viral_hits.blast > ~/map2.txt
# sort  map2.txt | uniq -u > crc_map2.txt
map2=pd.read_csv('~/crc_map3.txt',names=['contig','YP'],index_col=None,sep='\t')
MAP2=pd.merge(map1,map2)
MAP2['contig1']=MAP2['contig']
MAP2['contig']= MAP2.contig.str.split("_")[0]+"_"+MAP2.contig.str.split("_")[1]
MAP2.to_csv('~/MAP2_crc.txt',sep='\t')
# MAP2=pd.read_csv('~/MAP2_crc.txt',sep='\t',index_col=None)#,nrows=10000)
# MAP2.contig=MAP2.contig.str.split('_').str[0]+'_'+MAP2.contig.str.split('_').str[1]
del MAP2

# R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/idxstats.txt')
R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt')

for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=0,header=None)
    BB=pd.DataFrame(AA[2]+AA[3])
    jeff=j.split('/')[7]#jeff=os.path.basename(j).split('.')[0]
    BB.rename(columns={0:jeff},inplace=True)
    if i==0:
        idxstat=pd.DataFrame(BB)
    else:
        idxstat=pd.merge(idxstat,BB,how='outer',left_index=True, right_index=True)

# idxstat.to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/merge_idxstats.txt',sep='\t',header=True,index=True)
idxstat.to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/merge_idxstats.txt',sep='\t',header=True,index=True)

# idxstat=pd.read_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/merge_idxstats.txt',header=0,index_col=0,sep='\t')

idxstat=idxstat.rename_axis('contig').reset_index()#drop=True)

idxstat.['yp_count']=idxstat.groupby('contig').count()
aa=idxstat.sort_values(by='yp_count').drop_duplicates(['contig'],keep='last')
# idxstat['sum']=idxstat.sum(axis=1,numeric_only = True)
# aa=idxstat.sort_values(by='sum').drop_duplicates(['contig'],keep='last')

# idx_out=pd.merge(MAP2,aa)#,left_on='contig',right_index=True)


def preprocess(x):
    idx_out=pd.merge(aa,MAP2)#, left_on = "Colname1", right_on = "Colname2")
    # if not os.path.exists('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt') :
    if not os.path.exists('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_idx_table2.txt') :

        # idx_out.to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt',sep='\t',index=True, header=True)
        idx_out.to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_idx_table2.txt',sep='\t',index=True, header=True)

    else:
        # idx_out.to_csv("/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt",mode="a",header=False,index=True,sep='\t')
        idx_out.to_csv("/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_idx_table2.txt",mode="a",header=False,index=True,sep='\t')


# reader = pd.read_csv("/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt", chunksize=1000) # chunksize depends with you colsize
MAP2=pd.read_csv('~/MAP2_crc.txt',sep='\t',index_col=None,nrows=1000000)
MAP2.contig=MAP2.contig.str.split('_').str[0]+'_'+MAP2.contig.str.split('_').str[1]
MAP2=MAP2.drop(columns=['ff','YP','Unnamed: 0'])

[preprocess(r) for r in MAP2]

# idx_out.to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/new_idx_table.txt',sep='\t',index=True, header=True)

# idx_out=pd.read_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/new_idx_table.txt',sep='\t',low_memory=False)


# dd=pd.read_csv('~/data/annot_contigs.txt',sep='\t',index_col=0)
jj=dd.groupby(['contig','species']).sum()
jj['count']=dd.groupby(['contig','species']).count()['50131']
jj.reset_index(inplace=True)
idx=jj.groupby(['contig'])['count'].transform(max) == jj['count']
jj=jj[idx]

jj['sum']=jj.drop(columns=['count']).sum(axis=1)
kk=jj[~jj.contig.duplicated(keep='first')]
kk.drop(columns=['count','sum']).to_csv('annot_tbl_fin.txt',sep='\t')


# idx_out['sum']=idx_out.sum(axis=1,numeric_only = True)
# aa=idx_out.sort_values(by='sum').drop_duplicates(['contig'],keep='last')
# idx_out.drop(columns='sum').to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt',sep='\t',index=True, header=True)

zip -R assemble '*/*.txt' '*/*.prots' '*/*.genes' '*/*.blast' '*.tsv'