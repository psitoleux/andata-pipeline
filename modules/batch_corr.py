from bbknn import bbknn
import scanpy as sc


def bbknn(adata : sc.AnnData, batch_keys : list,n_neighbors : int ) -> sc.AnnData:
    
    
    return adata

def bbknn_reg(adata : sc.AnnData, batch_keys : list,n_neighbors : int ) -> sc.AnnData:
    
    
    return adata 



"""

Park et al (2020)'s BBKNN-reg implementation utils

their actual call 
bdata, X_explained = regress_iter(adata,['method','source','donor','Cycle_score'],['anno_predict_B'],'donor_method',scale=False)

"""


def regress_batch_v2(adata,batch_key,confounder_key):
    '''batch regression tool
    batch_key=list of observation categories to be regressed out
    confounder_key=list of observation categories to be kept
    returns ndata with corrected X'''

    from sklearn.linear_model import Ridge
    
    dummy = pd.get_dummies(adata.obs[batch_key+confounder_key],drop_first=False)
    X_exp = adata.X # scaled data
    if scipy.sparse.issparse(X_exp):
        X_exp = X_exp.todense()
    LR = Ridge(fit_intercept=False,alpha=1.0)
    LR.fit(dummy,X_exp)

    if len(batch_key)>1:
        batch_index = np.logical_or.reduce(np.vstack([dummy.columns.str.startswith(x) for x in batch_key]))
    else:
        batch_index = np.vstack([dummy.columns.str.startswith(x) for x in batch_key])[0]
    
    dm = np.array(dummy)[:,batch_index]
    X_explained = dm.dot(LR.coef_[:,batch_index].T)
    X_remain = X_exp - X_explained
    ndata = sc.AnnData(X_remain)
    ndata.obs = adata.obs
    ndata.var = adata.var
    return ndata, X_explained

def regress_iter(adata,batch_key,confounder_key,bbknn_key,scale=True, approx = True,n_pcs=50):
    if scale == True:
        sc.pp.scale(adata,max_value=10)
    ndata, X_explained = regress_batch_v2(adata,batch_key=batch_key,confounder_key=confounder_key)
    sc.pp.pca(ndata)
    bbknn(ndata, batch_key = bbknn_key,n_pcs=n_pcs, approx=approx)
    return ndata, X_explained
