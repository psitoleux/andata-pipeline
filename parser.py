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
    
    parser.add_argument("--shuffle", action=argparse.BooleanOptionalAction)
    parser.add_argument("--obs_mask", type=literal_eval, help = "dict observations to keep", default="{ }")
     
    parser.add_argument("--normalization", type=str, help="normalization", default="log1p")
    parser.add_argument("--pearson_overdis", type = float, help="pearson_overdispersion parameter", default=100)
    
    parser.add_argument("--feature_selection", type=str, help="feature_selection", default="hvg")
    parser.add_argument("--n_top_genes", type=int, help="n_top_genes", default=2000)
    
    parser.add_argument("--ambient", type=str, help="ambient", default=None)
    
    parser.add_argument("--doublets", type=str, help="doublets", default=None)
    parser.add_argument("--filter_dbl", action=argparse.BooleanOptionalAction)
    
    parser.add_argument("--dim_red", type=str, help="dim_reduction", default="pca")
    parser.add_argument("--pca_n_comps", type=int, help="pca_n_comps", default=50)
    
    parser.add_argument("--batch_corr", type=str, help="batch_corr", default=None)
    parser.add_argument("--batch_key", type=str, help="batch_key", default="method")
    
    parser.add_argument("--viz", type=str, help="visualization", default="umap")

    
    return parser

def parse_args() -> Namespace:

    return get_parser().parse_args()

def get_config() -> dict:
    
    return vars(parse_args())