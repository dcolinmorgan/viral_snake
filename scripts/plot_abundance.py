import seaborn as sns, numpy as np, pandas as pd
from numpy import inf

def plot_violin_and_box():
    data=pd.read_csv('assemble/contig_counts.tsv',sep='\t')
    data=data.drop(columns='length').melt(['contig']).replace({'-DNA.idxstats.txt':' '}, regex=True)

    plt.figure(figsize=(8,12))
    data['log_value']=np.log10(data.value)
    data['log_value'][np.isneginf(data['log_value'])] = 0
    data['log_abundance']=data['log_value'].fillna(0)
    tmp=sns.violinplot(data=data, y="variable", x="log_abundance",scale="count")

    plt.figure(figsize=(8,12))
    tmp=sns.boxplot(data=data, y="variable", x="log_abundance",hue='label',dodge=False)

    data['over_fifty'] = np.where( data.age > 50, 1, 0)
    plt.figure(figsize=(8,12))
    tmp=sns.boxplot(data=data, y="variable", x="log_abundance",hue='over_fifty',dodge=False)
    tmp.savefig('box_plot.png')

def abundance_plot():
    data=pd.read_csv('assemble/contig_counts.tsv',sep='\t')
    data=data.drop(columns='length').melt(['contig']).replace({'-DNA.idxstats.txt':' '}, regex=True)
    data.variable.replace(' ','',regex=True,inplace=True)
    label=pd.read_csv('assemble/label.txt',sep='\t')
    data=pd.merge(data,label,left_on='variable',right_on='id')
    species=pd.read_csv('assemble/dvf_species.txt',sep='\t',names=['contig','species'])
    data=pd.merge(data,species,how='left')
    data=data.fillna('unknown')
    data=data[data.value!=0]

    plt.figure(figsize=(8,12))
    data['log_abundance']=np.log10(data.value)
    data['log_abundance'][np.isneginf(data['log_abundance'])] = 0
    data['log_abundance']=data['log_abundance'].fillna(0)
    data['genus']=data.species.str.split(' ').str[0]

    data1=data[data.label=='Non-refractory']
    sub_gen=data1.groupby('genus').count().sort_values('contig',ascending=False).reset_index().genus[1:20]
    data1=data1[data1['genus'].isin(sub_gen)]

    plt.figure(figsize=(12,8))
    ax = sns.histplot(data1, x='variable', hue='genus', weights='log_abundance',
                multiple='stack', palette='tab20c', shrink=0.8)
    ax.set_ylabel('log_abundance')
    ax.set_xlabel('non-refractory')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

    # Fix the legend so it's not on top of the bars.
    legend = ax.get_legend()
    legend.set_bbox_to_anchor((1, 1))

    data2=data[data.label=='Refractory']
    sub_gen=data2.groupby('genus').count().sort_values('contig',ascending=False).reset_index().genus[1:20]
    data2=data2[data2['genus'].isin(sub_gen)]

    plt.figure(figsize=(12,8))
    ax = sns.histplot(data2, x='variable', hue='genus', weights='log_abundance',
                multiple='stack', palette='tab20c', shrink=0.8)
    ax.set_ylabel('log_abundance')
    ax.set_xlabel('refractory')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

    # Fix the legend so it's not on top of the bars.
    legend = ax.get_legend()
    legend.set_bbox_to_anchor((1, 1))
    plt.savefig('abundance_plot.png')

def plot_index():
    import ecopy
    data=pd.read_csv('assemble/contig_counts.tsv',sep='\t',index_col=0)
    data=data.drop(columns='length')#.replace({'-DNA.idxstats.txt':' '}, regex=True)
    data.columns=data.columns.str.replace('-DNA.idxstats.txt','',regex=True)
    data=data.T
    label=pd.read_csv('assemble/label.txt',sep='\t')
    data=pd.merge(data,label,left_index=True,right_on='id')
    data.index=data.label
    data=data.drop(columns={'id','sex','age','label'})
    A=ecopy.diversity(data, method='shannon')
    C=ecopy.diversity(data, method='simpson',breakNA=True)
    A=pd.DataFrame(A)
    A.index=data.index
    A.rename(columns={0:'Shannon_Index'},inplace=True)

    from statannot import add_stat_annotation

    ax=sns.boxplot(data=A,y=A.Shannon_Index,x=A.index)
    add_stat_annotation(ax, data=A,y=A.Shannon_Index,x=A.index,
                        box_pairs=[("Refractory", "Non-refractory")],# ("Thur", "Sat"), ("Fri", "Sun")],
                        test='t-test_ind', text_format='star', loc='outside', verbose=2)


    C=pd.DataFrame(C)
    C.index=data.index
    C.rename(columns={0:'Simpson_Index'},inplace=True)

    from statannot import add_stat_annotation
    C=C.replace(inf,0)
    ax=sns.boxplot(data=C,y=C.Simpson_Index,x=C.index)
    add_stat_annotation(ax, data=C,y=C.Simpson_Index,x=C.index,
                        box_pairs=[("Refractory", "Non-refractory")],# ("Thur", "Sat"), ("Fri", "Sun")],
                        test='t-test_ind', text_format='star', loc='outside', verbose=2)
    ax.savefig('eco_plot.png')