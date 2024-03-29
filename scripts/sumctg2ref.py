import sys, glob, argparse,os
from functools import cmp_to_key
from joblib import Parallel,delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()
# parser = argparse.ArgumentParser()
# parser.add_argument('acc2tax',type=str)
# parser.add_argument('ctg_num_orf',type=str)
# parser.add_argument('blasthits',type=str)
# args = parser.parse_args()

#python sumctg2ref.py acc2tax ctg2num_orf blastpres
#python snakemake_viral_calling/scripts/sumctg2ref.py snakemake_viral_calling/acc2species /groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/50310/idxstats.txt /groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/50310/blast_contigs.fa2
acc2family = {}
os.chdir('/home/dcmorgan')

def readtax(filename):
    global acc2family
    global id2name
    taxfile = open(filename,'r')
    for line in taxfile:
        line = line.strip()
        info = line.split('\t')
        acc = info[0].split('.')[0]
    
        name = info[2]
        acc2family[acc] =name
        

def readnorf(filename):
    ctg2norf = {}
    for line in open(filename,'r'):
        line = line.strip()
        info = line.split('\t')
        ctg2norf[info[0]] = int(info[2])
    return ctg2norf
def compare(anno1,anno2):
    if anno1[1] < anno2[1]:
        return -1
    elif anno1[1] == anno2[1]:
        if anno1[2]<anno2[2]:
            return -1
        elif anno1[2] == anno2[2]:
            return 0
        else:
            return 1
    else:
        return 1

def checktax(ctg,assignments,ctg2norf):
    if len(assignments) == 0:
        return 'not'
    norf = ctg2norf[ctg]
    tax2num = {}
    tax2score = {}

    num = 0
    for tax,score in assignments:
        num += 1
        if tax in tax2num.keys():
            tax2num[tax] += 1
            tax2score[tax] += score
        else:
            tax2num[tax] = 1
            tax2score[tax] = score

    #tax2numlist = list(tax2num.items())
    tax2anno = []
    for tax,num in tax2num.items():
        tax2anno.append([tax,num,tax2score[tax]])
    tax2anno = sorted(tax2anno,key=cmp_to_key(compare),reverse=True)

    if tax2anno[0][1]/norf >= 0.5:

        return tax2anno[0][0]
    else:
        return 'not'

def sumtab(filename):
    dirn=os.path.dirname(filename)
    readtax('snakemake_viral_calling/acc2species')
    ctg2norf = readnorf(dirn+'/idxstats.txt')
    global acc2family
    
    tabfile = open(dirn+'/blast_contigs.fa2','r')
    ctg2tax = {}
    for line in tabfile:
        line = line.strip()
        info = line.split('\t')
        orfinfo = info[0]
        ctginfo = orfinfo.split('_')
        ctg = ctginfo[0]+'_'+ctginfo[1]
        ctgid = ctg
        if not ctgid in ctg2norf.keys():
            continue
        acc = info[1].split('.')[0]
        score = float(info[2])     
        if not acc in acc2family.keys():
            continue
        assignment = [acc2family[acc],score]
        if ctgid in ctg2tax.keys():
            ctg2tax[ctgid].append(assignment)
        else:
            ctg2tax[ctgid] = [assignment]

    for ctg,assignments in ctg2tax.items():
        family = checktax(ctg,assignments,ctg2norf)
        if family != 'not':
            with open(dirn+"/ctg_spc.txt", "a") as f:
                print(ctg+'\t'+family,file=f)
                f.close()
# readtax(args.acc2tax)
# ctg2norf = readnorf(args.ctg_num_orf)
# sumtab(args.blasthits,ctg2norf)

# for i in glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt'):
    # dirn=os.path.dirname(i)
    # readtax('snakemake_viral_calling/acc2species')
    # ctg2norf = readnorf(dirn+'/idxstats.txt')
    # sumtab(dirn+'/blast_contigs.fa2',ctg2norf,dirn)
    
Parallel(n_jobs=num_cores)(delayed(sumtab)(i) for i in glob.glob('/groups/cgsd/dcmorgan/CRC_data/no_duplicated/assemble/*/idxstats.txt'))

# # for i in np.arange(1,len(data)):
# #     puc_cal(data,i,genes,str(tissue),str(omic))
# # # with open('data/Pipeline_consolidate_220301/gcn/'+tissue+'_'+omic+'cdf_ALL.pkl', 'wb') as f:
# # #     pickle.dump(cdf, f)
# # # with open('data/Pipeline_consolidate_220301/gcn/'+tissue+'_'+omic+'pucN_ALL.pkl', 'wb') as f:
# # #     pickle.dump(pucN, f)

    