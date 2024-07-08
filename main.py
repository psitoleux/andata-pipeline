from pipeline import Pipeline
from parser import get_config


if __name__ == "__main__":
    
    config = get_config()
    pipeline = Pipeline(config)
    pipeline.preprocess()
    pipeline.visualize()
    pipeline.analysis()
