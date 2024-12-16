def merge_glmm(folder,out):
    '''
    Given a folder containing Libra output for each cell type
    Merge the tables to produce one table containing the result of all cell types
    '''
    import pandas as pd
    import os
    files=os.listdir(folder)
    files=map(lambda x: os.path.join(folder,x), files)
    dfs=map(lambda x:pd.read_csv(x), files)
    merged=pd.concat(dfs,ignore_index=True)
    merged.to_csv(out,index=False)

def filter_l2fc(fi,fo,l2fc=0.1,padj=0.05):
    '''
    Given fi a (merged) Libra output
    Keep only rows with abs(L2FC) over l2fc treshold, and adjusted pval under padj threshold
    Save output to file fo
    '''
    import pandas as pd
    df=pd.read_csv(fi,index_col=0,header=0)
    #filter base on adjusted p
    df_clean=df[df['p_val_adj']<=padj]
    #filter base on abs(l2fc)
    df['abs_avg_logFC']=abs(df['avg_logFC'])
    df_clean=df[df['abs_avg_logFC']>=l2fc]
    #remove temp abs column
    df_clean.drop(columns='abs_avg_logFC',inplace=True)
    df_clean.to_csv(fo,index=False)