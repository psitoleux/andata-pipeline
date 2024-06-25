from argparse import ArgumentParser


def get_config() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument("--input", type=str, help="input")
    return parser