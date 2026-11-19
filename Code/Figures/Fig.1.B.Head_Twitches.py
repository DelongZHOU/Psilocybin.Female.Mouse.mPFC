def head_twitch(file='..//..//Data//Fig.1.B.headtwitches.csv',sex='F',out=None):
    import seaborn as sns, matplotlib.pyplot as plt, pandas as pd
    import numpy as np
    if out==None:
        out='..//..//Figures//Fig.1.B.headtwitches'
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    plt.rcParams['figure.figsize'] = [2.5,2]
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    data=pd.read_csv(file)
    df=data[data['Sex']==sex]
    sns.barplot(x='Dose', y='HTS', data=df,
                hue='Dose',legend=False, #to cope with palette
                errorbar='se',palette=['grey','#F18DBF','#EC60A6','#E6308A','#9F1359'])
    sns.swarmplot(x='Dose', y='HTS',data=df,
                    color='0', alpha=1)
    #set x axis frequency and labels
    plt.xticks([0,1,2,3,4],
               ['Saline','0.25','0.5','1.0','2.0'],
               rotation=45)
    #set y axis range
    plt.ylim([0,40])
    #set y axis tick frequency
    plt.yticks([0,10,20,30,40])
    #set title
    title={'F':'Female','M':'Male'}[sex]
    plt.title(title,y=1)
    #set axis titles
    plt.ylabel('Total Head Twitches')
    plt.xlabel('Dose (mg/kg)',fontsize=12)
    plt.subplots_adjust(left=0.2,top=0.85,bottom=0.4)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()