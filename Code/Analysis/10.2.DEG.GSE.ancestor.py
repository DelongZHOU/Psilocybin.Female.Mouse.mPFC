def count_gse_ancestor(gse,head=50,out):
    '''
    Given a gseresult table
        with columns including ID (GO ID), Celltype and Database (GO=<flavor>)
    count how many times each ancestors appear in the ancestors list of terms of each cell type
    output a table
    '''
    import pandas as pd
    import re
    import requests, sys, json
    gse=pd.read_csv(gse)
    celltypes=gse['Celltype'].unique()
    def gse_to_ancestors(gse,celltype):
        '''
        Given a GSE table
        Given a celltype
        First slice the GSE table by celltype then remove reactome
        For each GO, retrieve its ancestor from the d_ancestors dictionary
        count each ancestor in a dictionary
        return the dictionary as a pandas DF
        '''
        print(celltype)
        #slice by celltype and remove reactome
        df_ct=gse.query("Celltype=='{}' and Database!='reactome'".format(celltype)).copy()
        df_ct.sort_values('p.adjust',inplace=True)
        go_terms=df_ct.head(head)['ID'].unique()
        d={}
        for go in go_terms:
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates".format(go)
            r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            responseBody = r.text
            ancestors=re.findall('GO:[0-9]*',responseBody)
            for ancestor in ancestors:
                if ancestor!=go:
                    try:
                        d[ancestor]+=1
                    except KeyError:
                        d[ancestor]=1
        df_ancestors=pd.DataFrame.from_dict(d,orient='index',columns=[celltype])
        return df_ancestors
    def get_go_name(go):
        '''
        Given a GO term, get the description of the term
        '''
        import requests, sys, json
        requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/terms/"+go
        r = requests.get(requestURL, headers={ "Accept" : "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        responseBody = r.text
        return json.loads(responseBody)['results'][0]['name']
    dfs=map(lambda x: gse_to_ancestors(gse,x), celltypes)
    df_merge=pd.concat(dfs,axis=1)
    df_merge.fillna(0)
    df_merge['GO Name']=df_merge2.apply(lambda row: get_name(row.name), axis=1)
    df_merge.to_csv(out)