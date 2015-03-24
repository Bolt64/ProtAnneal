#!/usr/bin/env python

import get_atom_distances as gad
import itertools as it
import pdb_parser as parser

lj_potential = lambda distance, rm: (rm/distance)**12 - 2*(rm/distance)**6

def get_energy(pdb_file, rm_file):
    energy=0
    rm_values=gad.read_labeled_pairs(rm_file)
    for atom1,atom2 in it.combinations(parser.parse_file(pdb_file), 2):
        key=tuple(sorted(((atom1['resName'], atom1['name']),(atom2['resName'], atom2['name']))))
        if key in rm_values:
            energy+=lj_potential(gad.get_distance(atom1, atom2), rm_values[key])
    return energy
