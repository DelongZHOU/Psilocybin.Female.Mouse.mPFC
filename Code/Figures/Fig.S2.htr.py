def plot_htr(file1='..//..//Data//Fig.S2.htr2a.csv',file2='..//..//Data//Fig.S2.htr.csv',out='..//..//Figures//Figure S2-2'):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '12'
    htr2a=pd.read_csv(file1)
    htr=pd.read_csv(file2)
    genes=['Htr2c','Htr1a','Htr3a','Htr1b','Htr7','Htr5a']
    fig = plt.figure(figsize=(8,10))
    gs=fig.add_gridspec(4,2)
    ax=fig.add_subplot(gs[0,:])
    #use first row to plot Htr2a in all cell types
    sns.violinplot(data=htr2a, x='Cell_Type', y='Expression',inner=None,
                   ax=ax, color='b',
                   order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA',' ','Astro','Endo','Micro','Oligo','OPC'],
                   linewidth=0)
    sns.stripplot(data=htr2a, x='Cell_Type', y='Expression',
                  size=2, jitter=True,
                  ax=ax,color='k',
                  order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA',' ','Astro','Endo','Micro','Oligo','OPC'],
                  )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', rotation=45)
    ax.set_ylim(0, 3)
    ax.set_xlabel('Cell Types')
    ax.set_title('Htr2a')
    #add panel letter
    ax.text(x=-1.5,y=3.5,s='A',fontsize=18)
    #on the left 3 panels plot excitatory receptors 2c, 3a, 7
    #on the right 3 panels plot inhibitory receptors 1a, 1b, 5a
    for i, gene in enumerate(genes):
        ax = fig.add_subplot(gs[i//2+1,i%2])
        gene_data = htr[htr['Gene'] == gene]
        sns.violinplot(data=gene_data, x='Cell_Type', y='Expression',inner=None,
                       ax=ax,  color='b',
                       order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA'],
                       linewidth=0)
        sns.stripplot(data=gene_data, x='Cell_Type', y='Expression',
                      size=2, jitter=True,
                      ax=ax, color='k',
                      order=['L2-6 IT','L5 PT','L5/6 NP','L6 CT','GABA'],
                      )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(f'{gene}')
        ax.set_xlabel('')
        ax.set_ylim(0, 5)
        #add panel letter
        ax.text(x=-1.5,y=5.5,s='BECFDG'[i],fontsize=18,fontweight='bold')
        #add y axis title for the middle panels
        ax.set_ylabel('Expression')
        #add x ticks
        ax.tick_params(axis='x', rotation=45)
        if i>3:
            ax.set_xlabel('Neuron Subtypes')
    fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()