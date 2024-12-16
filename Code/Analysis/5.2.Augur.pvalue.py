def compute_pval(augur,raw='raw',permute='permute',out='.'):
    '''
    Given augur file
    Assume it's under raw folder
    And a permute auc metric file under permute folder with ending .permute.csv
    add pval and adj pval to the raw table and save under out
    '''
    import pandas as pd
    from scipy import stats
    df_raw=pd.read_csv(os.path.join(raw,augur),index_col=0)
    df_permute=pd.read_csv(os.path.join(permute,augur[:-4]+'.permute.csv'))
    def celltype_compare(celltype,auc,df_permute):
        '''
        Given a celltype, its auc under AUGUR without permute, and the permute
        Given ther permute table from permute run
        find how many permuted auc is larger than the real auc and the number of iterations
        calculate the p val with pseudo count
        '''
        #slice by celltype
        df_permute_celltype=df_permute[df_permute['cell_type']==celltype]
        #slice by metric == auc
        df_permute_celltype=df_permute_celltype[df_permute_celltype['metric']=='roc_auc']
        #take average of folds for each iteration
        permute_auc=df_permute_celltype.groupby('subsample_idx')['estimate'].mean()
        #count numbers of auc larger than the real auc
        m=permute_auc[permute_auc>auc].count()
        #count the total number of iterations
        n=permute_auc.count()
        #return the permute pval with pseudocount
        return (m+1)/(n+1)
    df_raw['pval']=df_raw.apply(lambda row: celltype_compare(row['cell_type'],row['auc'],df_permute), axis=1)
    df_raw['qval']=stats.false_discovery_control(df_raw['pval'])
    df_raw.to_csv(os.path.join(out,augur),index=False)