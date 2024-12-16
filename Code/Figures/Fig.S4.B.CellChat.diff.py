def cellchat_diff_heatmap(fi_1='..//..//Data//Fig.S4.B.1.CellChat.SAL.csv',fi_2='..//..//Data//Fig.S4.B.2.CellChat.PSI.csv',out='..//..//Figures//Fig.S4.B.CellChat.Diff.Heatmap',na='NE'):
    '''
    Take two csv of CellChat interaction numbers
    Both csv should be the same dimention
    Create the difference matrix by substracting element wise
    plot heatmap
    '''
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable
    from math import floor
    import pandas as pd
    import seaborn as sns
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    vmin=-80
    vmax=20
    #set color
    matplotlib.rcParams['figure.figsize']=3.5,3
    matplotlib.rcParams['text.color'] = 'b'
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    sns.set(color_codes=True)
    df1=pd.read_csv(fi_1,index_col=0,na_values=na).fillna(0)
    df2=pd.read_csv(fi_2,index_col=0,na_values=na).fillna(0)
    data=df2-df1
    plt.figure(figsize=(3.5,3))
    color=LinearSegmentedColormap.from_list(name=' ',colors=['grey','white','#E6308A'])
    ax=sns.heatmap(data,
                   linewidths=0.5,cmap=color,
                   yticklabels=True,
                   center=0,vmin=vmin,vmax=vmax,
                   cbar_kws={'shrink': 0.6}
                   )
    ax.set(xlabel='Target', ylabel='Source')
    ax.set_title('PSI - SAL')
    make_axes_area_auto_adjustable(ax)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.subplots_adjust(left=0,right=1,bottom=0,top=1)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()
    