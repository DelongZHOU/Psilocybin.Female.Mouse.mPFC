def gse_ancestor_heatmap(fi='..//..//Data//Fig.6.B.Selected.Ancestors.csv',out='..//..//Figures//Fig.6.B.GSE.Ancestor.Heatmap',na='NE'): #just heatmap
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable
    from math import floor
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap
    sns.set(color_codes=True)
    matplotlib.rcParams['text.color'] = 'b'
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    df=pd.read_csv(fi,sep=",",index_col=0)
    df.columns=['Nervous System\nDevelopment','Ribosome','Mitochondrion']
    df=df.transpose()
    fig, ax = plt.subplots(figsize=(5,2))
    color=LinearSegmentedColormap.from_list(name=' ',colors=['grey','white','#E6308A'])
    #data=df.ix[:,:].dropna() #drop any rows with NaN
    ax=sns.heatmap(df,
            cmap=color,linewidths=0.1,
            yticklabels=True,
            vmin=-10,vmax=15,center=0,
            cbar_kws={'ticks':[-10,-5,0,5,10,15],'shrink': 0.6}
            )
    ax.collections[0].colorbar.set_ticklabels(['10 (↑SAL)','5','0','5','10','15 (↑PSI)'])
    ax.set_xticklabels(ax.get_xticklabels(),fontsize=12,rotation=45)
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=12)
    plt.subplots_adjust(left=0.3,bottom=0.35,right=0.9)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()