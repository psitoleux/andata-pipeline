from abmodule import Module
import scanpy as sc


def scdbl(adata : sc.AnnData) -> sc.AnnData:
    
    data_mat = adata.X.T

    %%R -i data_mat -o doublet_score -o doublet_class

    set.seed(123)
    sce = scDblFinder(
        SingleCellExperiment(
            list(counts=data_mat),
        ) 
    )
    doublet_score = sce$scDblFinder.score
    doublet_class = sce$scDblFinder.class

    adata.obs["scDblFinder_score"] = doublet_score
    adata.obs["scDblFinder_class"] = doublet_class
    adata.obs.scDblFinder_class.value_counts()
    
    return adata