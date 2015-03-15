#!/usr/bin/env python

import pdb_parser as parser
import itertools as it
import math
import os
import get_buried_residues as gbr

atoms_to_look = (" C", " H", " S")

def get_distance(atom1, atom2):
    return math.sqrt(sum((atom1[i]-atom2[i])**2 for i in 'xyz'))

def separate_by_res_no(pdb_file):
    separated_residues = {}
    for atom in parser.parse_file(pdb_file):
        if not atom['resSeq'] in separated_residues:
            separated_residues[atom['resSeq']] = []
        separated_residues[atom['resSeq']].append(atom)
    return separated_residues

def get_pairwise_distances(res1, res2):
    for atom1 in res1:
        if atom1['element'] in atoms_to_look:
            for atom2 in res2:
                if atom1['element'] in atoms_to_look:
                    yield get_distance(atom1, atom2)

def make_pairs(separated_residues):
    for res1, res2 in it.combinations(separated_residues, 2):
        yield res1, res2

def get_distances(pdb_file):
    separated_residues = separate_by_res_no(pdb_file)
    for res1, res2 in make_pairs(separated_residues):
        for distance in get_pairwise_distances(separated_residues[res1], separated_residues[res2]):
            yield distance

def get_distances_from_core(pdb_file):
    os.popen("/home/bolt/protein_lab/third-party/binaries/naccess2.1.1/naccess "+pdb_file)
    temp_output=open("temp.pdb", "w")
    filename=pdb_file.split("/")[-1][:-3]
    for l in gbr.get_inner_residues(pdb_file, gbr.get_buried_residues(filename+"rsa")):
        temp_output.write(l+"\n")
    temp_output.close()
    for distance in get_distances("temp.pdb"):
        print(distance)
    os.popen("rm temp.pdb")
    os.popen("rm {0}asa".format(filename))
    os.popen("rm {0}rsa".format(filename))
    os.popen("rm {0}log".format(filename))

if __name__ == "__main__":
    from sys import argv
    get_distances_from_core(argv[-1])
