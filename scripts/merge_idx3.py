import glob,os,sys,re
import pandas as pd

os.chdir('/home/dcmorgan')

DS=['epilepsy/assemble/*/idxstats.txt',
    '/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/idxstats.txt',
    '/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt',
    '/groups/cgsd/dcmorgan/HT/*/idxstats.txt'
   ]


map1=pd.read_csv('~/long.species',sep='\t',header=None)
map1[0]=map1[0].str.split('#').str[0]
# map1.groupby('contig').count()
map1.rename(columns={0:'contig',1:'species'},inplace=True)
# map1.to_csv('~/contig2spec.txt',sep='\t')


# R=glob.glob('epilepsy/assemble/*/idxstats.txt')
# R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/*/idxstats.txt')
R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt')

for i,j in enumerate(R):
    kk=j.split('/')[7]
    # AA=pd.read_csv('epilepsy/assemble/'+kk+'/idxstats.txt',sep='\t',header=None,index_col=0)
    # AA=pd.read_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/assemble/'+kk+'/idxstats.txt',sep='\t',header=None,index_col=0)
    AA=pd.read_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/'+kk+'/idxstats.txt',sep='\t',header=None,index_col=0)
    BB=pd.DataFrame(AA[2]*1000*1000000/(AA[1]*1000000))
    BB.rename_axis('contig',inplace=True)
    BB.rename(columns={0:'50045'},inplace=True)
    CC=pd.merge(map1,BB.reset_index())
    # CC.to_csv("epilepsy/annot/"+kk+"_idx_table_M4.txt",header=True,index=False,sep='\t')
    # CC.to_csv("/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot/"+kk+"_idx_table_M4.txt",header=True,index=False,sep='\t')
    CC.to_csv("/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/"+kk+"_idx_table_M4.txt",header=True,index=False,sep='\t')

# R=glob.glob('epilepsy/annot/*idx_table_M4.txt')
# R=glob.glob('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot/*idx_table_M4.txt')
R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot/*idx_table_M4.txt')

for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=None)
    SS=j.split('annot')[1].split('_')[0].strip('/')
    AA.columns=['contig','species',SS]
    if i!=0:
        BB=AA.merge(BB,on=['contig','species'],how='outer')#,left_on=['contig','species'])
    else:
        BB=AA
# BB.fillna(0).to_csv('epilepsy/annot_contigs_M4.txt',sep='\t')
# BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_contigs_M4.txt',sep='\t')
BB.fillna(0).to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/annot_contigs_M4.txt',sep='\t')
