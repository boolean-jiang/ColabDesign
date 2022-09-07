import numpy as np
import pandas as pd
import sys, os
import argparse
import math, random
from IPython.display import HTML
from colabdesign import mk_afdesign_model, clear_mem

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="pdb filename or filepath", type=str)
    parser.add_argument("-c", "--chainid", help="pdb chain ID of target(s)", type=str)
    parser.add_argument("-l", "--length", help="length of binder", type=int)
    parser.add_argument("-r", "--residues", help="hotspot residues on target chain to optimize binding", type=str)
    parser.add_argument("-d", "--design", help="design protocol", type=int, choices=[2,3])
    parser.add_argument("-s", "--steps", help="comma-separated list for # design steps", type=str)

    args = parser.parse_args()

    model = mk_afdesign_model(protocol="binder")
    model.prep_inputs(pdb_filename=args.pdb, chain=args.chainid, binder_len=args.length, hotspot=args.residues)
    
    steps = [int(step) for step in args.steps.strip().split(',')]

    if args.design == 2:
        model.design_2stage(steps[0], steps[1])
    else:
        model.design_3stage(steps[0], steps[1], steps[2])
    
    print("-"*30 + '\n')
    print("Best design:", model.get_seqs())
    print("Last design:", model.get_seqs(get_best=False))

if __name__ == "__main__":
    main()
