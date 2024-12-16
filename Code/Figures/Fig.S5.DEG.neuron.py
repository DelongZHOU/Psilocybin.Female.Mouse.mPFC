def plot_deg_violin_strip(file='..//..//Data//Fig.S5.neuron.csv',out='..//..//Figures//Figure S5'):
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
    annotates={'Kctd13':[4,'**'],'Myo6':[4,'**'],'Nr4a3':[0,'****'],'Nrn1':[2,'**'],'Shank3':[2,'**'],'Slitrk3':[4,'**']}
    genes = long_data['Gene'].unique()
    fig, axes = plt.subplots(4, 2, figsize=(8,10), sharex=False, sharey=False)
    axes = axes.ravel()
    for i, gene in enumerate(genes):
        ax = axes[i]
        gene_data = long_data[long_data['Gene'] == gene]
        sns.violinplot(data=gene_data, x='Cell_Type', y='Expression',
                       hue='Condition', split=True, inner=None, 
                       ax=ax, 
                       palette={'SAL': 'gray', 'PSI': '#E6308A'}, 
                       hue_order=['SAL', 'PSI'],
                       order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA'],
                       linewidth=0)
        sns.stripplot(data=gene_data, x='Cell_Type', y='Expression', 
                      hue='Condition', dodge=True,
                      size=2, jitter=True, 
                      ax=ax, 
                      palette='dark:black',
                      hue_order=['SAL', 'PSI'],
                      order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA'],
                      )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(f'{gene}', color='#FFA500')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(0, 6)
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
        if gene=='Arc': #mannual annotate since 2 cell types
            ax.annotate('**', xy=(0, 5), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
            ax.annotate('***', xy=(3, 5), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
        elif gene=='Gria1':
            ax.annotate('*', xy=(0, 5), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
            ax.annotate('#', xy=(1, 5), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
        elif gene in annotates:
            x,star=annotates[gene]
            ax.annotate(star, xy=(x, 5), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
        ax.get_legend().remove()
    fig.supxlabel('Neuron Subtypes')
    fig.supylabel('Expression')
    legend_elements = [Line2D([0], [0], color='gray', lw=4, label='SAL'),
                       Line2D([0], [0], color='#E6308A', lw=4, label='PSI')
                       ]
    fig.legend(handles=legend_elements, loc='center right', fontsize=12, title='Conditions')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()