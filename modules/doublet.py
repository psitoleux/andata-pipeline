from abmodule import Module
import scanpy as sc


def scdbl(adata: sc.AnnData, seed : int = 123) -> sc.AnnData:
    # Import R libraries
    scDblFinder = importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")

    # Prepare the data matrix
    data_mat = adata.X.T

    # Transfer data to the R environment
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["seed"] = seed

    # Run scDblFinder in R
    ro.r('''
        library(scDblFinder)
        library(SingleCellExperiment)

        set.seed(seed)
        sce <- scDblFinder(
            SingleCellExperiment(
                list(counts = data_mat)
            )
        )
        doublet_score <- sce$scDblFinder.score
        doublet_class <- sce$scDblFinder.class
    ''')

    # Retrieve the results
    doublet_score = np.array(ro.globalenv["doublet_score"])
    doublet_class = np.array(ro.globalenv["doublet_class"])

    # Store the results in the AnnData object
    adata.obs["scDblFinder_score"] = doublet_score
    adata.obs["scDblFinder_class"] = doublet_class
    print(adata.obs["scDblFinder_class"].value_counts())
    
    
    
 
    return adata

