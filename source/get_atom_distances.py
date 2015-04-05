#!/usr/bin/env python

"""
This defines the functions to get inter atomic
distances in the protein.
"""

import pdb_parser as parser
import itertools as it
import math
import os
from config import *

# Some global constants. I'm so sorry. :(
#ATOMS_TO_LOOK = (" C",)# " H", " S")
#THRESHOLD = 10 # Angstroms
#TO_CONSIDER = (' CD1',)# " CD2")
#RESIDUE = "LEU"
#naccess_location = "/home/bolt/protein_lab/third-party/binaries/naccess2.1.1/naccess"
#reduce_location = "/home/bolt/protein_lab/third-party/binaries/reduce"


def get_distance(atom1, atom2):
    """
    Returns distance between atom1 and atom2.
    """
    return math.sqrt(sum((atom1[i]-atom2[i])**2 for i in 'xyz'))

def get_buried_residues(rsa_file, threshold=30):
    """
    Takes an rsa file and returns the residues which less
    than threshold accessible.
    """
    for line in open(rsa_file):
        if line.startswith("RES"):
            temp=line.strip().split()
            res_no=int(temp[3])
            rel_accessibility=float((temp[7]))
            if rel_accessibility<=threshold:
                yield res_no

def get_buried_residues_pdb(pdb_file, threshold=30):
    """
    Takes a pdb file and returns a list of residue numbers
    which are buried more than the threshold
    """
    filename = pdb_file.split("/")[-1][:-4]
    os.popen("{0} {1}".format(naccess_location, pdb_file))
    result = list(get_buried_residues(filename+".rsa", threshold))
    os.popen("rm {}.rsa".format(filename))
    os.popen("rm {}.asa".format(filename))
    os.popen("rm {}.log".format(filename))
    return result

def get_inner_residues(prot_pdb, list_of_res):
    """
    Takes in a prot_pdb file and a list of residues
    and returns a pdb file with only the selected residues.
    """
    res_nos={i for i in list_of_res}
    for atom in parser.parse_file(prot_pdb):
        if atom['resSeq'] in res_nos:
            yield parser.pdb_format(atom)

def filter_by_residue(residue, carbon_atoms):
    """
    Higher order function that takes in a residue and carbon atoms
    and returns a function that only evaluates to true if the atom
    satisfies the criteria.
    """
    def predicate(line):
        if line['resName']==residue and line['name'] in carbon_atoms:
            return True
        else:
            return False
    return predicate

def get_distances(pdb_file):
    """
    Gets the pairwise distances of all the valid atoms in the protein.
    """
    predicate=filter_by_residue(RESIDUE, TO_CONSIDER)
    for atom1,atom2 in it.combinations(parser.filter_line(parser.parse_file(pdb_file), predicate), 2):
        yield get_distance(atom1, atom2)

def filter_distance(list_of_distances, predicate):
    """
    Takes a list of distances and filter through
    according to predicate.
    """
    for distance in list_of_distances:
        if predicate(distance):
            yield distance

def distance_threshold(threshold):
    """
    A higher order function that returns a predicate which
    is only true if the distance is less than threshold.
    """
    def predicate(distance):
        if distance<=threshold:
            return True
        else:
            return False
    return predicate

def get_labeled_pairs(pdb_file):
    labeled_pairs={}
    for atom1, atom2 in it.combinations(parser.parse_file(pdb_file), 2):
        label=tuple(sorted(((atom1['resName'], atom1['name']),(atom2['resName'], atom2['name']))))
        labeled_pairs[label]=get_distance(atom1, atom2)
    return labeled_pairs

def write_labeled_pairs(labeled_pairs, output_file):
    with open(output_file, "w") as output:
        for pair in labeled_pairs:
            res1,res2=pair
            out_string="{0} {1} {2} {3} {4}\n".format(res1[0], res1[1], res2[0], res2[1], labeled_pairs[pair])
            output.write(out_string)

def read_labeled_pairs(pair_file):
    labeled_pairs={}
    for line in open(pair_file):
        a=line[0:3]
        b=line[4:8]
        c=line[9:12]
        d=line[13:17]
        e=float((line.strip().split()[-1]))
        labeled_pairs[((a,b),(c,d))]=e
    return labeled_pairs

def get_distances_from_core(pdb_file):
    """
    This function does everything, reduces the protein, gets the core
     and then the distances. Very ugly though. :P
    """
    os.popen("{0} {1} > temp1.pdb".format(reduce_location, pdb_file))
    os.popen("{0} temp1.pdb".format(naccess_location))
    temp_output=open("temp2.pdb", "w")
    for l in get_inner_residues("temp1.pdb", get_buried_residues("temp1.rsa")):
        temp_output.write(l+"\n")
    temp_output.close()
    for distance in filter_distance(get_distances("temp2.pdb"), distance_threshold(THRESHOLD)):
        print(distance)
    os.popen("rm temp1.pdb")
    os.popen("rm temp2.pdb")
    os.popen("rm temp1.rsa")
    os.popen("rm temp1.asa")
    os.popen("rm temp1.log")

if __name__ == "__main__":
    from sys import argv
    get_distances_from_core(argv[-1])

"""
Discarded functions

def get_subunit(subunit):
    def predicate(line):
        if line['chainID']==subunit:
            return True
        else:
            return False
    return predicate

def separate_by_res_no(pdb_file, subunit):
    separated_residues = {}
    for atom in parser.filter_line(parser.parse_file(pdb_file), get_subunit(subunit)):
        if atom['resName'] == "LEU" and atom['name'] in TO_CONSIDER:
            if not (atom['chainID'],atom['resSeq']) in separated_residues:
                separated_residues[(atom['chainID'], atom['resSeq'])] = []
            separated_residues[(atom['chainID'], atom['resSeq'])].append(atom)
    return separated_residues

def get_pairwise_distances(res1, res2):
    for atom1 in res1:
        if atom1['element'] in ATOMS_TO_LOOK:
            for atom2 in res2:
                if atom1['element'] in ATOMS_TO_LOOK:
                    yield get_distance(atom1, atom2)

def make_pairs(separated_residues):
    for res1, res2 in it.combinations(separated_residues, 2):
        yield res1, res2

"""
