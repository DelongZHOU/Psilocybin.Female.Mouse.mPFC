def cell_count(file='..//..//Data//Fig.S1.E.cell.count.csv',out='..//..//Figures//Fig.S1.E.cell.count'):
    import seaborn as sns, matplotlib.pyplot as plt, pandas as pd
    colors=['#c6dbef','#9ecae1','#3182bd','#08519c','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['patch.edgecolor'] = 'none'
    plt.rcParams["figure.figsize"] = (3,3)
    data=pd.read_csv(file,index_col=0)
    celltypes=data.columns
    #calculate the total cell count for each sample
    data['total']=data.sum(axis=1)
    #calculate the ratio of each cell type in every sample 
    #in order to plot normalized stacked bar plot
    for celltype in celltypes:
        data[celltype]=data[celltype]/data['total']
    #remove the total entry
    data=data.drop('total',axis=1)
    f = plt.figure()
    data.plot.bar(stacked=True,ax=f.gca(),color=colors)
    #set legend outside of the plot area
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),
                fontsize=8,
                handlelength=1,
                handleheight=0.5,
                handletextpad=0.5,
                labelspacing=0.3,
                borderpad=0.3,
                borderaxespad=0.3,
                )
    plt.xticks(rotation=45)
    #increase margin for long sample name and outside legend
    plt.subplots_adjust(top=0.9,bottom=0.35,right=0.75,left=0.15)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()