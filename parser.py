import argparse
from argparse import ArgumentParser, Namespace
from ast import literal_eval
from configparser import ConfigParser



def get_parser() -> ArgumentParser:
    parser = ArgumentParser()
    
    parser.add_argument("--input", type=str, help="input"
                        , default="HTA07.A01.v02.entire_data_raw_count.h5ad")  
    parser.add_argument("--output", type=str, help="output", default="./out/")
    parser.add_argument("--save", action=argparse.BooleanOptionalAction)
    
    parser.add_argument("--seed", type=int, help="seed", default=1998)
    
    parser.add_argument("--outlier_keys", type=literal_eval, help="outlier_keys", default="['mt']")
    parser.add_argument("--qc", action=argparse.BooleanOptionalAction)
    
    
    
     
    parser.add_argument("--normalization", type=str, help="normalization", default="log1p")
    parser.add_argument("--feature_selection", type=str, help="feature_selection", default="hvg")
    parser.add_argument("--ambient", type=str, help="ambient", default=None)
    parser.add_argument("--doublets", type=str, help="doublets", default=None)
    
    parser.add_argument("--dim_reduction", type=str, help="dim_reduction", default="pca")
    parser.add_argument("--pca_n_comps", type=int, help="pca_n_comps", default=50)
    
    parser.add_argument("--batch_corr", type=str, help="batch_corr", default=None)
    parser.add_argument("--batch_key", type=str, help="batch_key", default="method")
    
    parser.add_argument("--visualization", type=str, help="visualization", default="umap")

    
    return parser

def parse_args() -> Namespace:

    return get_parser().parse_args()

def get_config() -> dict:
    
    return vars(parse_args())