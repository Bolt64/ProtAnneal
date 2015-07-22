#!/usr/bin/env python

import pdb_parser as parser
import get_atom_distances as gad
from annealing import get_current_residue

def get_sequence(pdb_file):
    buried_residues = gad.get_buried_residues_pdb(pdb_file)
    buried_residues.sort()
    return [get_current_residue(pdb_file, i) for i in buried_residues]

def similarity(seq1, seq2):
    assert len(seq1)==len(seq2), "Sequences not of same length"
    matches = sum(1.0 for i in range(len(seq1)) if seq1[i]==seq2[i])
    return matches/len(seq1)

if __name__=="__main__":
    from sys import argv
    print similarity(get_sequence(argv[-1]), get_sequence(argv[-2]))
