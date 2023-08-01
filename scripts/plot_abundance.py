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
    
def vir_top10():
    
    data=pd.read_csv('/groups/cgsd/dcmorgan/PRJNA544527/V2_results/annot_idx_table.txt',sep='\t',usecols=range(1,540))
    data.fillna(0,inplace=True)
    # data
    # !wget https://raw.githubusercontent.com/dcolinmorgan/grph/main/PRJNA544527-meta_inf.txt

    meta=pd.read_csv('PRJNA544527-meta_inf.txt',sep='\t',header=None)
    meta=meta.rename(columns={3:'id',5:'pat'})

    data.index=data.species
    data=data.T

    mm=pd.merge(data,meta[['id','pat']],left_index=True,right_on='id')
    # mm=mm.T

    # mm['contig']=data.contig
    # mm=mm.T
    # mm.species=data.T.index
    mm['id']=mm['pat'].str.split('-').str[0]
    mm['time']=mm['pat'].str.split('_').str[0].str.split('-').str[1]

    # !wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0559-3/MediaObjects/41591_2019_559_MOESM3_ESM.xlsx --no-check-certificate
    metaa=pd.read_excel('41591_2019_559_MOESM3_ESM.xlsx',sheet_name='SupTable2',skiprows=3)
    metaa=metaa[['Donor','Age','Sex','BMI']]


    Full_table=pd.merge(mm,metaa,left_on='id',right_on='Donor')
    # Full_table=Full_table.drop(columns=[3,	'pat',	'id'])
    # Full_table.Label=pd.to_dateLabel(Full_table.Label,unit='d')
    # Full_table.Label=Full_table.Label.values.astype('dateLabel64[M]')

    data2=Full_table.melt(id_vars=['time','Donor','Age','Sex','BMI','id','pat'])#,'species','contig'])
    # data2=data2.rename(columns={'variable':'species'})
    data2=data2.sort_values(by=['Donor','time','value'])
    # data2=data2[data2.value>0]
    # data2.to_csv('PRJNA544527_mpa4_annot_table.txt',sep='\t')

    # final df stored here also
    # !wget https://raw.githubusercontent.com/dcolinmorgan/grph/main/PRJNA544527_mpa4_annot_table.txt
    # data2=pd.read_csv('PRJNA544527_mpa4_annot_table.txt',sep='\t',index_col=0)

    data2=data2[data2.value>0]
    data2=data2.reset_index(drop = True)
    # data2=data2.groupby(['Donor','Label','value'],group_keys=False).apply(lambda x : x.sort_values(by = 'value', ascending = False).head(10))
    data2=data2.drop_duplicates()

    data2["Label"] = (
        data2.groupby("Donor")
        .apply(lambda x: x.groupby("time", sort=False).ngroup() + 1)
        .values
    )


    cc=pd.unique(data2[data2.Label<5].Donor)
    # len(pd.unique(data2.Donor))
    data2=data2.loc[ data2.Donor.isin(cc), : ]
    data2=data2[data2.Label<5]
    data2=data2[data2['rank']<10.0]
    data2.to_csv('PRJNA544527_viral10top.txt',sep='\t')
    
def plot_vir_top10():
    import random
    plt.figure(figsize=(18,8))

    chars = '0123456789ABCDEF'
    palette_colors=['#'+''.join(random.sample(chars,6)) for i in range(len(pd.unique(data2.variable)))]
    # palette_colors = sns.color_palette('tab20c')
    palette_dict = {genus: color for genus, color in zip(pd.unique(data2.variable), palette_colors)}
    data2["rank"] = data2.groupby("Donor")["value"].rank(method="dense", ascending=False)
    data2=data2[data2['rank']<10.0]

    ax = sns.histplot(data2[data2['Sex']=='Male'], x='Donor', hue='variable', weights='value',
                 multiple='stack', palette=palette_dict, shrink=0.8)
    ax.set_ylabel('abundance')
    # ax.set_xlabel('non-refractory')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    ax.set_title('Male')

    plt.figure(figsize=(12,8))

    ax = sns.histplot(data2[data2['Sex']=='Female'], x='Donor', hue='variable', weights='value',
                 multiple='stack', palette=palette_dict, shrink=0.8)
    ax.set_ylabel('abundance')
    # ax.set_xlabel('non-refractory')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    ax.set_title('Female')
