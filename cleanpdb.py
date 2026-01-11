#!/usr/bin/env python3
import os
import argparse
import sys
import openmmawsem
import helperFunctions.myFunctions


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__author__ = 'Wei Lu'

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create\
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein(1r69 for example): \
                        python3 ~/OPENAWSEM_LOCATION/mm_create_project.py 1r69 or you can specify the target pdb file you want\
                        to simulation. for example do: python3 ~/OPENAWSEM_LOCATION/mm_create_project.py YOUR_PDB_LOCATION/1r69.pdb")
parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False, help="Generate fragment memories")
parser.add_argument("--extended", action="store_true", default=False, help="Start from extended structure")
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--hybrid", action="store_true", default=False)
parser.add_argument("--verbose", action="store_true", default=False)
parser.add_argument("--predict_ssweight_from_fasta", action="store_true", default=False)
parser.add_argument("--keepIds", action="store_true", default=False, help="Set to True if you want to preserve the chain and residue index. default will rename chains from 'A', and index from 1")
parser.add_argument("--keepLigands", action="store_true", default=True)
args = parser.parse_args()

# Print if in debug
if args.debug:
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# Log the command to a file
with open('create_project_commandline_args.txt', 'w') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

# if you provide the pdb then we will use it for the project(move to folder named original_pdbs to prevent from overwritten), otherwise we download it online.
if args.protein[-4:] == '.pdb':
    if not os.path.exists(args.protein):
        print("ERROR: the pdb you specified is not exist")
        exit()
    name = os.path.basename(args.protein)[:-4]
    pdb = os.path.basename(args.protein)
    do("mkdir -p original_pdbs")
    do(f"cp {args.protein} original_pdbs/")
print (name)
removeHeterogens = False if args.keepLigands is True else True
chain = args.chain

if not os.path.exists(f"crystal_structure.pdb"):
    helperFunctions.myFunctions.cleanPdb([name], chain=chain, toFolder="cleaned_pdbs", verbose=args.verbose, keepIds=True, removeHeterogens=removeHeterogens)
    do(f"cp cleaned_pdbs/{pdb} crystal_structure.pdb")


# If the chain is not specified then select all the chains
if chain == "-1":
    chain = helperFunctions.myFunctions.getAllChains("crystal_structure.pdb")
    print("Chains info read from crystal_structure.pdb, chains to simulate: ", chain)

# for compute Q
input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb("crystal_structure.pdb", chain,keepIds=args.keepIds, removeHeterogens=True)
openmmawsem.ensure_atom_order(input_pdb_filename)


# get fasta, pdb, seq file ready
chain = helperFunctions.myFunctions.getAllChains("crystal_structure-cleaned.pdb")
openmmawsem.getSeqFromCleanPdb(input_pdb_filename, chains=chain, writeFastaFile=True)
do(f"cp crystal_structure.fasta {name}.fasta")

input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb(pdb, chain, keepIds=args.keepIds, removeHeterogens=removeHeterogens)
#input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb("crystal_structure.pdb", chain,keepIds=args.keepIds, removeHeterogens=removeHeterogens)
openmmawsem.ensure_atom_order(input_pdb_filename)
if args.keepLigands:
    cleaned_pdb_filename = f"{name}-cleaned.pdb"
    input_pdb_filename = f"{name}-openmmawsem.pdb"
    do(f"grep 'ATOM' {input_pdb_filename} > tmp.pdb")
    do(f"grep 'HETATM' {cleaned_pdb_filename} >> tmp.pdb")
    do(f"mv tmp.pdb {input_pdb_filename}")
    do(f"sed -i \"/TER/d\" crystal_structure-openmmawsem.pdb")
    do(f"sed -i \"/END/d\" crystal_structure-openmmawsem.pdb")
    do(f"grep 'HETATM' {cleaned_pdb_filename} >> crystal_structure-openmmawsem.pdb")

