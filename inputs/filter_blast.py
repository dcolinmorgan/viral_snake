import pandas as pd
import glob,os
import subprocess


# rc = call("./sleep.sh")
os.chdir('/home/dcmorgan/')

for i in glob.glob('assemble/unfinished/*'):
    os.chdir(i)
    # try:
    reads=pd.read_csv('final_contigs.fa',names=['name','seq'])

    hits=pd.read_csv('viral_hits.blast',sep='\t',names=['name','tmp','identity','tmp1','tmp2','tmp3','tmp4','tmp5','tmp8','tmp12','tmp13','tmp11'])
    hits=hits[['name','identity']]

    AA=set(reads.name.str.split(' ').str[0]) & set('>'+hits.name.str.split('_').str[0]+'_'+hits.name.str.split('_').str[1])

    reads[reads['name'].str.split(' ').str[0].isin(AA)].to_csv('blast_contigs.fa',header=False,index=False)

    # except:
    #     print("didn't run: "+i)
    os.chdir('/home/dcmorgan/')