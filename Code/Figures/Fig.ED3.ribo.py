def plot_deg_violin_strip(file='..//..//Data//Fig.ED3.Ribosome.genes.csv',out='..//..//Figures//Fig.ED3.ribo'):
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
    annotates={'Hsp90b1':[9,'*'],'Nsa2':[8,'*'],'Ubr1':[1,'***']}
    genes = long_data['Gene'].unique()
    fig, axes = plt.subplots(2, 2, figsize=(8,10), sharex=False, sharey=False)
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.5, top=0.9, wspace=0.2, hspace=0.4)
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
        ax.set_title(f'{gene}', color='#4ea72e')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(0, 5)
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
        if gene in annotates:
            x,star=annotates[gene]
            ax.annotate(star, xy=(x, 5), xycoords='data',
                ha='center', va='center', fontsize=8, color='k')
        ax.get_legend().remove()
    ax.tick_params(axis='x', rotation=45)
    fig.supxlabel('Cell Types',y=0.5)
    fig.supylabel('Expression',y=0.8)
    legend_elements = [Line2D([0], [0], color='gray', lw=4, label='SAL'),
                       Line2D([0], [0], color='#E6308A', lw=4, label='PSI')
                       ]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=12, title='Conditions',bbox_to_anchor=(1, 0.82))
    plt.tight_layout(rect=[0, 0.5, 0.85, 1])

    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()