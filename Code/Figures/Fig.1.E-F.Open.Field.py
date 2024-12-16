def open_field_zscore(file='..//..//Data//Fig.1.E-F.Open.Field.csv',sex='F',out=None):
    import seaborn as sns, matplotlib.pyplot as plt, pandas as pd
    import numpy as np
    if out==None:
        panel={'F':'E','M':'F'}[sex]
        out='..//..//Figures//Fig.1.{}.Open.Field.{}'.format(panel,sex)
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    plt.rcParams['figure.figsize'] = [1.5,2.5]
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    data=pd.read_csv(file)
    df=data[data['Sex']==sex]
    sns.barplot(x='Treatment', y='Z_score', data=df,
                hue='Treatment',legend=False, #to cope with palette
                errorbar='se',palette=['grey','#E6308A'])
    sns.swarmplot(x='Treatment', y='Z_score',data=df,
                    color='0', alpha=1)
    #set y axis range
    plt.ylim([-4,4])
    #set y axis tick frequency
    plt.yticks(np.arange(-4,5, 2.0))
    plt.xticks(rotation=45)
    #set title
    title={'F':'Female','M':'Male'}[sex]
    plt.title(title,y=1.2)
    #set axis titles
    plt.ylabel('Z_score(Time in Center)')
    plt.xlabel('Treatment',fontsize=12)
    plt.subplots_adjust(left=0.4,top=0.8,bottom=0.4)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()