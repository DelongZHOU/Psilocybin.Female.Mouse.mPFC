def neuronchat_bar(file='..//..//Data//Fig.3.E.neuronchat.summary.csv',out='..//..//Figures//Fig.S4.D.neuronchat.bar'):
    import seaborn as sns, matplotlib.pyplot as plt, pandas as pd
    import numpy as np
    df = pd.read_csv(file)
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    plt.rcParams['figure.figsize'] = [2,2.5]
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    df=pd.read_csv(file)
    sns.barplot(x='Treatment', y='Strength', data=df,
                hue='Treatment',legend=False, #to cope with palette
                errorbar='se',palette=['grey','#E6308A'])
    plt.ylim([0,15])
    plt.xticks(rotation=45)
    plt.ylabel('Interaction Strength')
    plt.xlabel('Treatment',fontsize=12)
    plt.subplots_adjust(left=0.3,top=0.9,bottom=0.4)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()