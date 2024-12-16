def plot_deg_violin_strip(file='..//..//Data//Fig.S7.Mito.genes.csv',out='..//..//Figures//Figure S7'):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    data=pd.read_csv(file)
    long_data = data.melt(id_vars=['Cell_ID', 'Cell_Type', 'Condition'], 
                      var_name='Gene', 
                      value_name='Expression')
    #manual encode the highlight for significance
    genes = long_data['Gene'].unique()
    annotates={'Baalc':[4,'**'],'Bag4':[4,'**'],'Idh3a':[4,'***'],'Epha4':[0,'****'],'Galc':[1,'*'],'Ppp1cb':[1,'**'],'Slc25a53':[1,'**'],'S100a1':[2,'****'],'Micu3':[3,'**'],'Snca':[3,'*'],'Pi4k2a':[8,'***'],'Ppp3ca':[9,'****'],'Sqstm1':[9,'****']}
    fig, axes = plt.subplots(7, 2, figsize=(8,10), sharex=False, sharey=False)
    axes = axes.ravel()
    for i, gene in enumerate(genes):
        ax = axes[i]
        gene_data = long_data[long_data['Gene'] == gene]
        sns.violinplot(data=gene_data, x='Cell_Type', y='Expression',
                       hue='Condition', split=True, inner=None, 
                       ax=ax, 
                       palette={'SAL': 'gray', 'PSI': '#E6308A'}, 
                       hue_order=['SAL', 'PSI'],
                       order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA',' ','Astro','Endo','Micro','Oligo','OPC'],
                       linewidth=0)
        sns.stripplot(data=gene_data, x='Cell_Type', y='Expression', 
                      hue='Condition', dodge=True,
                      size=2, jitter=True, 
                      ax=ax, 
                      palette='dark:black',
                      hue_order=['SAL', 'PSI'],
                      order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA',' ','Astro','Endo','Micro','Oligo','OPC'],
                      )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(f'{gene}', color='#386cb0')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(0, 5)
        if gene in annotates:
            x,star=annotates[gene]
            ax.annotate(star, xy=(x, 4), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
        ax.get_legend().remove()
        ax.tick_params(axis='x', rotation=45)
    fig.supxlabel('Cell Types')
    fig.supylabel('Expression')
    legend_elements = [Line2D([0], [0], color='gray', lw=4, label='SAL'),
                       Line2D([0], [0], color='#E6308A', lw=4, label='PSI')
                       ]
    fig.legend(handles=legend_elements, loc='center right', fontsize=12, title='Conditions')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()