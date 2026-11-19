def deg_volcano_plot(csv='..//..//Data//Fig.2.F.1.Libra.csv',highlight='..//..//Data//Fig.2.F.2.Libra.highlight.csv',yaxis=True,p=1.3,l2fc=0.1):
    '''
    csv: a csv table of Different Expression Analysis result from Libra; contains DE result for all cell types / subtypes
    prefix: output png file
    highlight: a text file containing list of genes to high light for each cell type
    yaxis: whether to plot the y axis or not
    p: cutoff of neg. log_10 of adjusted p val. 1.3 for p-val<=0.05
    l2fc: cutoff of log2 fold change
    '''
    import math
    import pandas as pd
    import plotly.graph_objects as go
    import plotly.express as px
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from adjustText import adjust_text
    import warnings
    warnings.filterwarnings('ignore')
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '10'
    plt.rcParams['figure.figsize'] = [2,2]
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['text.color'] = 'k'
    #by default use input csv name as prefix
    if yaxis:
        out_folder='..//..//Figures//Fig.2.F.volcano'
    else:
        out_folder='..//..//Figures//Fig.2.F.volcano.no_y'
    os.makedirs(out_folder, exist_ok=True)
    prefix='Fig.2.F.DEG.volcano'
    #load data
    hl_df=pd.read_csv(highlight,sep=',')
    def row_to_color(row):
        if row.significant=='Yes':
            return {'Nervous System Development':'#d01c8b','Mitochondrion':'#386cb0','Ribosome':'#4ea72e'}[row.function]
        else:
            return 'k'
    hl_df['color']=hl_df.apply(row_to_color,axis=1)
    df_full=pd.read_csv(csv,sep=',',header=0)
    df_full['neg_log_p']=df_full['p_val_adj'].apply(lambda x: -math.log(max(x,10**-10),10))
    #extract cell types
    cell_types=df_full['cell_type'].unique()
    #make volcano plot for each cell type
    for cell_type in cell_types:
        #slice highlight table
        hl_df_ct=hl_df[hl_df['celltype']==cell_type]
        colors=dict(zip(hl_df_ct['gene'],hl_df_ct['color']))
        #prepare output save path
        out_png='{}.{}.png'.format(prefix,cell_type).replace('/','_').replace(' ','_')
        out_pdf='{}.{}.pdf'.format(prefix,cell_type).replace('/','_').replace(' ','_')
        #subset dataframe to the cell type
        df=df_full[df_full['cell_type']==cell_type]
        plt.scatter(df.avg_logFC, df.neg_log_p,s=1,color='black')
        ##plot down regulated genes with grey
        deg_neg=df.query('(avg_logFC<=-{0}) & neg_log_p>={1}'.format(l2fc,p))
        plt.scatter(deg_neg.avg_logFC, deg_neg.neg_log_p, color="grey",s=4)
        ##plot up regulated genes with pink
        deg_pos=df.query('(avg_logFC>={0}) & neg_log_p>={1}'.format(l2fc,p)).copy()
        plt.scatter(deg_pos.avg_logFC, deg_pos.neg_log_p, color="#E6308A",s=4)
        #highlight Gene id from highlight file
        df_name=df[df['gene'].isin(hl_df_ct['gene'])]
        # #reset highlight genes
        # texts=[]
        # for index, row in df_name.iterrows():
            # texts.append(plt.text(row.avg_logFC, row.neg_log_p, row.gene,c=colors[row.gene]))
        plt.xlim(-1,1)
        plt.xticks(np.arange(-1,1.5, 0.5))
        #plt.ylabel('-Log10(adjusted p value)')
        plt.ylim(0,10.5)
        plt.yticks(np.arange(0,10.5, 2))
        plt.title(cell_type)
        if not yaxis:
            plt.gca().get_yaxis().set_visible(False)
            plt.gca().spines['left'].set_visible(False)
        # adjust_text(texts)
        plt.subplots_adjust(left=0.3,bottom=0.25)
        plt.savefig(os.path.join(out_folder,out_png),dpi=600)
        plt.savefig(os.path.join(out_folder,out_pdf),dpi=600)
        plt.close()