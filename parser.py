from argparse import ArgumentParser, Namespace
from configparser import ConfigParser


def get_parser() -> ArgumentParser:
    parser = ArgumentParser()
    
    parser.add_argument("--input", type=str, help="input")
    
    parser.add_argument
    
    return parser

def parse_args() -> Namespace:

    return get_parser().parse_args()

def get_config() -> dict:
    
    return vars(parse_args())