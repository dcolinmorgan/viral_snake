import glob,os,sys,re
import pandas as pd

## make proper annotations
## based on /groups/cgsd/dcmorgan/virusdb/viral.1.protein.faa
map1=pd.read_csv('~/map.txt',names=['YP','species','ff'],index_col=None,sep='[')

map1.YP=map1.YP.str.split('.1').str[0]
map1.YP=map1.YP.str.replace('>','')
map1.species=map1.species.str.replace(']','')
map1.YP=map1.YP+'.1'
map1.drop(columns='ff').to_csv('~/MAP1.txt',sep='\t',index=False)

## based on cut -f1,2 /groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/viral_hits.blast > ~/map2.txt
map2=pd.read_csv('~/map2.txt',names=['contig','YP'],index_col=None,sep='\t')
MAP2=pd.merge(map1,map2)
# MAP2.to_csv('~/MAP2.txt',sep='\t')
# MAP2=pd.read_csv('~/MAP2.txt',sep='\t',index_col=None,nrows=10000)
MAP2.contig=MAP2.contig.str.split('_').str[0]+'_'+MAP2.contig.str.split('_').str[1]


R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*.idxstats.txt')
for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=0,header=None)
    BB=pd.DataFrame(AA[2]+AA[3])
    jeff=os.path.basename(j).split('.')[0]
    BB.rename(columns={0:jeff},inplace=True)
    if i==0:
        CC=pd.DataFrame(BB)
    else:
        CC=pd.merge(CC,BB,how='outer',left_index=True, right_index=True)
        
# CC.to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/merge_idxstats.txt',sep='\t',header=True,index=True)
# idxstat=pd.read_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/merge_idxstats.txt',header=0,index_col=0,sep='\t')
idxstat=CC
idx_out=pd.merge(MAP2,idxstat,left_on='contig',right_index=True)
idx_out.to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/new_idx_table.txt',sep='\t',index=True, header=True)

