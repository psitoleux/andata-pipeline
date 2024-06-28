from argparse import ArgumentParser, Namespace
from configparser import ConfigParser


def get_parser() -> ArgumentParser:
    parser = ArgumentParser()
    
    parser.add_argument("--input", type=str, help="input")
    parser.add_argument("--output", type=str, help="output")
    
    parser.add_argument("--normalization", type=str, help="normalization", default="log1p")
    parser.add_argument("--feature_selection", type=str, help="feature_selection", default="hvg")
    parser.add_argument("--ambient", type=str, help="ambient", default=None)
    parser.add_argument("--doublets", type=str, help="doublets", default="scdbl")
    parser.add_argument("--dim_reduction", type=str, help="dim_reduction", default="pca")
    parser.add_argument("--batch_corr", type=str, help="batch_corr", default=None)
    parser.add_argument("--visualization", type=str, help="visualization", default="umap")
     
    
    return parser

def parse_args() -> Namespace:

    return get_parser().parse_args()

def get_config() -> dict:
    
    return vars(parse_args())