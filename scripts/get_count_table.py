# import subprocess
# subprocess.Popen("for subdir in *; do cp {snakemake.input[0]} {snakemake.output[0]}; done;", shell=True).communicate()
# subprocess.Popen("find . -size 0 -delete", shell=True).communicate()

# def get_count_table(wildcards):
#     import sys,glob
#     flag = 0
#     ids = []
#     length = []
#     count = []
#     fname = [{snakemake.output[1]}]
#     R=glob.glob('assemble/{snakemake.output[0]}.idxstats.txt')
#     for f in R:
#         fname.append(f)
#         co = []
#         for n,line in enumerate(open(f)):
#             #print line
#             if(line[:1] == '*'):
#                 continue
#             spl = line.strip().split('\t')
#             # print spl[0]
#             if(flag == 0):
#                 ids.append(spl[0])
#                 length.append(spl[1])
#             su = int(spl[2])+int(spl[3])
#             co.append(str(su))
#         count.append(co)
#         co = []
#         flag = 1
#     names = '\t'.join(fname)
#     print('\t'.join(['contig','length',names]))

#     for i in range(0,len(ids)):
#         tco = []
#         for j in range(0,len(count)):
#             tco.append(count[j][i])
#         result =[ids[i],length[i],'\t'.join(tco)]
#         print('\t'.join(result))

import pandas as pd
R=glob.glob(snakemake.input[0])
for i,j in enumerate(R):
    AA=pd.read_csv(j,sep='\t',index_col=0,header=None)
    BB=pd.DataFrame(AA[2]+AA[3])
    jeff=os.path.basename(j).split('.')[0]
    BB.rename(columns={'0':jeff})
    if i==0:
        CC=pd.DataFrame(BB)
    else:
        CC=pd.merge(CC,BB,how='outer',left_index=True, right_index=True)
CC.to_csv(snakemake.output[0].str.strip('#').str[0],sep='\t')


# %%bash
# sed -n '/>/p' /groups/cgsd/dcmorgan/virusdb/viral.1.protein.faa > map.txt
# cut -f1,2 */viral_hits.blast > ~/map2.txt

map1.YP=map1.YP.str.split('.1').str[0]
map1.YP=map1.YP.str.replace('>','')
map1.species=map1.species.str.replace(']','')
map1.drop(columns='ff').to_csv('~/MAP1.txt',sep='\t',index=False)

map2=pd.read_csv('~/map2.txt',names=['contig','YP'],index_col=None,sep='\t')
MAPA=pd.merge(MAP1,map2)

idxstat=pd.read_csv('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/merge_idxstats.txt',header=0,index_col=0,sep='\t')

idx_out=pd.merge(MAPA,idxstat,left_on='contig',right_index=True)

idx_out.to_csv(new_idx_table.txt,sep='\t',index=True, header=True)