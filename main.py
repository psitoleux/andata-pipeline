from pipeline import Pipeline
from parser import get_config


if __name__ == "__main__":
    
    config = get_config()
    pipeline = Pipeline(config)
    pipeline.preprocess()
    #pipeline.directions_analysis()
    #pipeline.directions_analysis('donor')
    #pipeline.grid_direction_analysis(key = 'donor')
    #pipeline.visualize()
    pipeline.analysis(PCA=True, )
