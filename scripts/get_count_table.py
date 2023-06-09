import subprocess
subprocess.Popen("for subdir in *; do cp {snakemake.input[0]} {snakemake.output[0]}; done;", shell=True).communicate()
subprocess.Popen("find . -size 0 -delete", shell=True).communicate()

def get_count_table(wildcards):
    import sys,glob
    flag = 0
    ids = []
    length = []
    count = []
    fname = [{snakemake.output[1]}]
    R=glob.glob('assemble/{snakemake.output[0]}.idxstats.txt')
    for f in R:
        fname.append(f)
        co = []
        for n,line in enumerate(open(f)):
            #print line
            if(line[:1] == '*'):
                continue
            spl = line.strip().split('\t')
            # print spl[0]
            if(flag == 0):
                ids.append(spl[0])
                length.append(spl[1])
            su = int(spl[2])+int(spl[3])
            co.append(str(su))
        count.append(co)
        co = []
        flag = 1
    names = '\t'.join(fname)
    print('\t'.join(['contig','length',names]))

    for i in range(0,len(ids)):
        tco = []
        for j in range(0,len(count)):
            tco.append(count[j][i])
        result =[ids[i],length[i],'\t'.join(tco)]
        print('\t'.join(result))

