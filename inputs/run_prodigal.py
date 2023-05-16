import pandas as pd
import glob,os
import subprocess


# rc = call("./sleep.sh")
os.chdir('/home/dcmorgan/')

for i in glob.glob('assemble/*'):
    os.chdir(i)
    try:
        GMOREA0142=pd.read_csv('final.contigs.fa_gt500bp_dvfpred.txt',sep='\t')
        # print('A')
        vir_pred=GMOREA0142[GMOREA0142.pvalue<.05].name.str.split(' ').str[0]#.to_csv('filtered_contig.txt',sep='\t')
        # print('AB')
        # bashCommand = "paste -d ',' - - < final.contigs.fa >final_contigs.fa"
        # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        # output, error = process.communicate()
        subprocess.Popen("paste -d ',' - - < final.contigs.fa >final_contigs.fa", shell=True).communicate()

        # call("paste -d ',' - - < final.contigs.fa >final_contigs.fa")
        # print('B')
        cont=pd.read_csv('final_contigs.fa',sep=',',header=None)
        # print('C')
        contA=cont[0].str.split(' ').str[0].str.split('>').str[1]
        AA=set(contA) & set(vir_pred)
        # print('D')
        cont['contA']=contA
        contB=cont[cont.contA.isin(vir_pred)]

        del contB['contA']
        contB.to_csv('final_contigs.fna',sep='\n',header=False,index=False)
        # print('E')
        # bashCommand = "prodigal -i final_contigs.fna  -o prod.genes -a prod.prots"
        # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        # call("prodigal -i final_contigs.fna  -o prod.genes -a prod.prots")
        subprocess.Popen("prodigal -i final_contigs.fna  -o prod.genes -a prod.prots", shell=True).communicate()
        print("ran: "+i)
    
    except:
        print("didn't run: "+i)
    os.chdir('/home/dcmorgan/')