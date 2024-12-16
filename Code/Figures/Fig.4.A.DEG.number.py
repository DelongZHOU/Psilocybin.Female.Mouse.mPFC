def deg_number_barplot(csv='..//..//Data//Fig.6.A.Libra.count.csv',out='..//..//Figures//Fig.6.A.DEG.barplot'):
    '''
    From a matrix of celltype X up/down X count
    produce a two way bar plot
    with up = upregulated number of genes
    down = down regulated genes, count should be positive
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore")
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    plt.rcParams['figure.figsize'] = [3,2]
    df=pd.read_csv(csv,header=0)
    up=df[df['State']=='up']
    dn=df[df['State']=='down']
    dn['Neg_Count']=-dn['Count']
    fig=plt.figure()
    ax=plt.subplot(111)
    ax.bar(up['Cell_Type'],up['Count'],color='#E6308A')
    ax.bar(dn['Cell_Type'],dn['Neg_Count'],color='gray')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(bottom=False)
    plt.ylim([-40,60])
    plt.ylabel('Number of DEGs')
    plt.yticks([-40,-20,0,20,40,60],[40,20,0,20,40,60])
    plt.xticks(rotation=45)
    plt.subplots_adjust(left=0.1, right=1, top=0.9, bottom=0.3)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()