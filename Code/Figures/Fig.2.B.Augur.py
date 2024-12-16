def augur_barplot(file='..//..//Data//Fig.2.B.Augur.csv',out='..//..//Figures//Fig.2.B.Augur'):
    '''
    The Augur result was manually re-ordered to separate neurons and glias. An empty line was inserted between neurons and glias to produce a space between tew two.
    '''
    import pandas as pd
    import seaborn
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    plt.rcParams['figure.figsize'] = [4,2.5]
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    data=pd.read_csv(file,header=0,index_col=0)
    g=seaborn.barplot(data,x='Cell Type',y='Augur',
                        color='white',edgecolor="black")
    g.xaxis.set_tick_params(rotation=45)
    g.tick_params(bottom=False)
    plt.ylim(0.5, 0.65)
    loc = plticker.MultipleLocator(base=0.05)
    g.yaxis.set_major_locator(loc)
    plt.subplots_adjust(left=0.22,bottom=0.35)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',format='pdf',dpi=600)
    plt.close()