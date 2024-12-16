def volcano_plot_l56np(csv='..//..//Data//Fig.3.C.L56NP.marker.csv',out='..//..//Figures//Fig.3.C.L56NP.markers',n=10,p=1.3,l2fc=1):
    '''
    csv: a csv table of FindMarkers output
    out: output prefix
    n: number of top up/down regulated genes to highlight, default is 10
    p: cutoff of neg. log_10 of adjusted p val. 1.3 for p-val<=0.05
    l2fc: cutoff of log2 fold change, default is 1
    '''
    import math
    import pandas as pd
    import plotly.graph_objects as go
    import plotly.express as px
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    from adjustText import adjust_text
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = '10'
    plt.rcParams['patch.edgecolor'] = 'none'
    plt.rcParams['figure.figsize'] = [3,3]
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False
    #load data
    df=pd.read_csv(csv,sep=',',index_col=0,header=0)
    df['gene']=df.index
    #calculate the negtive log 10 p value; pval < 1e-300 including 0 are set at 1e-300
    df['neg_log_p']=df['p_val_adj'].apply(lambda x: -math.log(max(x,10**-300),10))
    plt.scatter(df.avg_log2FC, df.neg_log_p, alpha=0.5,s=8,color='grey')
    #plot down regulated genes with blue
    neg_markers=df.query('(avg_log2FC<=-{0}) & neg_log_p>={1}'.format(l2fc,p))
    plt.scatter(neg_markers.avg_log2FC, neg_markers.neg_log_p, color="white",s=8,edgecolor="red")
    #plot up regulated genes with red
    pos_markers=df.query('(avg_log2FC>={0}) & neg_log_p>={1}'.format(l2fc,p)).copy()
    plt.scatter(pos_markers.avg_log2FC, pos_markers.neg_log_p, color="red",s=8)
    #highlight top n genes with Gene id
    pos_markers.sort_values('avg_log2FC',ascending=False,inplace=True)
    n_markers=min(n,len(pos_markers.gene.tolist())) # reduce n if positive markers are fewer than n
    texts=[]
    for i in range(n_markers):
        x=pos_markers.avg_log2FC.tolist()[i]
        y=pos_markers.neg_log_p.tolist()[i]
        gene=pos_markers.gene.tolist()[i]
        if gene=='Htr2c': #highlight Htr2c with red
            texts.append(plt.text(x,y,gene,c='r',fontsize=10))
        else:
            texts.append(plt.text(x,y,gene,c='k',fontsize=10))
    #label axis and title
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(adjusted p value)')
    plt.xlim([-6,6])
    plt.ylim([0,320])
    loc = plticker.MultipleLocator(base=2)
    plt.gca().xaxis.set_major_locator(loc)
    adjust_text(texts, only_move={'points':'y', 'texts':'y'})
    plt.subplots_adjust(left=0.2,bottom=0.2)
    plt.savefig(out+'.png',dpi=600)
    plt.savefig(out+'.pdf',dpi=600)
    plt.close()