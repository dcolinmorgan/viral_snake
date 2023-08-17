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

# R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/idxstats.txt')
R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt')

for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=0,header=None)
    BB=pd.DataFrame(AA[2]+AA[3])
    jeff=j.split('/')[7]#jeff=os.path.basename(j).split('.')[0]
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
    map2=map2[[0,1]]
    map2.columns=['contig','YP']

# map2=pd.read_csv('~/crc_map3.txt',names=['contig','YP'],index_col=None,sep='\t')
    MAP2=pd.merge(map1,map2)
    MAP2['contig1']=MAP2['contig']
    MAP2['contig']= MAP2.contig.str.split('_').str[0]+'_'+MAP2.contig.str.split('_').str[1]



    idx_out=pd.merge(aa,MAP2)#, left_on = "Colname1", right_on = "Colname2")
    SS=RR.split('assemble')[1].strip('/')
    # idx_out


    jj=idx_out.groupby(['contig','species']).sum()
    jj['count']=idx_out.groupby(['contig','species']).count()['YP'].values
    jj.reset_index(inplace=True)
    idx=jj.groupby(['contig'])['count'].transform(max) == jj['count']
    jj=jj[idx]

    jj['sum']=jj.drop(columns=['count']).sum(axis=1)
    kk=jj[~jj.contig.duplicated(keep='first')]
    kk.drop(columns=['count','sum','yp_count']).to_csv("/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/"+SS+"_idx_table2.txt",header=True,index=False,sep='\t')


    
    R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/*idx_table2.txt')
    for i,j in enumerate(R):
        AA=pd.read_csv(j,sep='\t',index_col=None)
        SS=j.split('annot')[1].split('_')[0].strip('/')
        AA.columns=['contig','species',SS]
        if i!=0:
            BB=AA.merge(BB,on=['contig','species'],how='outer')#,left_on=['contig','species'])
        else:
            BB=AA
    BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_contigs.txt',sep='\t')