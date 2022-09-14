import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import sys, os
import argparse
import math, random
from colabdesign import mk_afdesign_model, clear_mem

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="pdb filename or filepath", type=str)
    parser.add_argument("-c", "--chainid", help="pdb chain ID(s) of target(s)", type=str)
    parser.add_argument("-l", "--length", help="length of binder", type=int)
    parser.add_argument("-m", "--use_multimer", help="Use AF-multimer", type=int, default=0)
    parser.add_argument("-r", "--hotspot_residues", help="hotspot residues on target chain to optimize binding", type=str)
    parser.add_argument("-g", "--semigreedy", help="use semigreedy design protocol", type=int, default=0)
    parser.add_argument("--logits", help="number of logit optimization steps", type=int, default=100)
    parser.add_argument("--soft", help="number of softmax(logits) optimization steps", type=int, default=100)
    parser.add_argument("--hard", help="number of one_hot(logits) optimization steps", type=int, default=10)
    parser.add_argument("--pdb_out_name", help="filename of design PDB structure", type=str)
    #parser.add_argument("--fas_out_name", help="filename of design FASTA file", type=str) 

    args = parser.parse_args()

    #using setup dictionary format to specify arguments, a la binder_hallucination example notebook
    setup = {"pdb_filename": args.pdb, "chain": args.chainid, "binder_len": args.length, "hotspot": args.hotspot_residues,}

    model = mk_afdesign_model(protocol="binder", use_multimer = args.use_multimer > 0)
    model.prep_inputs(**setup)

    #default hyperparams from binder hallucination example notebook - accept CLI arguments for these later 
    weights = {"con":0.0, "pae":0.0, "plddt":0.1, "i_pae":0.0, "i_con":1.0}
    opt = {
        "con":{"seqsep":9, "cutoff":14.0, "num":1, "binary":False},
        "i_con":{"cutoff":21.6875, "num":1, "num_pos":float("inf"), "binary":False}
    }

    model.restart(opt=opt, weights=weights)

    print("Script actually startingggg")

    #Add ability to disable residues, set seq, etc. with CLI arguments later 

    #accept CLI arguments for num_recycles later 
    if args.semigreedy > 0:
        model.design_semigreedy(args.semigreedy, num_recycles=0)
    else:
        model.design_3stage(args.logits, args.soft, args.hard, num_recycles=int(args.use_multimer>0))
    
    model.save_pdb(f"{args.pdb_out_name}.pdb")

    #For now, just return PDB file - later, figure out how to make both FASTA and PDB
    #with open(f"{args.fas_out_name}.fasta", "w") as outfas: outfas.write(model.get_seqs())
    #return f"{args.pdb_out_name}.pdb", ####f"{args.fas_out_name}.fasta", 

if __name__ == "__main__":
    main()
