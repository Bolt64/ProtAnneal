#!/usr/bin/env python

import pdb_parser as parser
import get_atom_distances as gad
import residue_mutator as mut
from config import *


def mutate_to_alanines(protein_pdb):
    protein_lines = list(parser.parse_file(protein_pdb))
    res_lines = list(parser.parse_file(residue_path.format("ALA")))
    buried_residues = gad.get_buried_residues_pdb(protein_pdb)
    current_pdb = protein_lines
    for i in buried_residues:
        current_pdb = mut.mutate(current_pdb, i, res_lines)
    return current_pdb

if __name__=="__main__":
    from sys import argv
    orig_protein = argv[-1]
    new_protein = mutate_to_alanines(orig_protein)
    for l in new_protein:
        print(parser.pdb_format(l))
