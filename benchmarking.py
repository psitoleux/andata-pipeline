from pipeline import Pipeline



def benchmark_normalization(adata, methods, config):
    
    for norm in methods:

        config["normalization"] = norm
        pipeline = Pipeline(config)
        pipeline.run(adata)
        
    # TODO do the relevant comparisons