import os
import argparse
from dotenv import load_dotenv
parser = argparse.ArgumentParser( description='scRNA Pipelines')
parser.add_argument( '--infile', type=str, help='input  file')
parser.add_argument( '--outdir', type=str, help='output file')
args = parser.parse_args()

# load env
load_dotenv()
cell_ranger = os.getenv('CellRanger')
