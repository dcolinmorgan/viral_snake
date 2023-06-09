import glob,os,sys,re
import pandas as pd
R=glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*.idxstats.txt')
for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=0,header=None)
    BB=pd.DataFrame(AA[2]+AA[3])
    jeff=os.path.basename(j).split('.')[0]
    BB.rename(columns={'0':jeff})
    if i==0:
        CC=pd.DataFrame(BB)
    else:
        CC=pd.merge(CC,BB,how='outer',left_index=True, right_index=True)
CC.to_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/merge_idxstats.txt',sep='\t',header=True,index=True)