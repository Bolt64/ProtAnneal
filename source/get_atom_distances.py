#!/usr/bin/env python

import pdb_parser as parser
import itertools as it
import math
import os

atoms_to_look = (" C", " H", " S")

def get_distance(atom1, atom2):
    return math.sqrt(sum((atom1[i]-atom2[i])**2 for i in 'xyz'))

def get_subunit(subunit):
    def predicate(line):
        if line['chainID']==subunit:
            return True
        else:
            return False
    return predicate

def get_buried_residues(rsa_file, threshold=30):
    for line in open(rsa_file):
        if line.startswith("RES"):
            temp=line.strip().split()
            res_no=int(temp[3])
            rel_accessibility=float((temp[7]))
            if rel_accessibility<=threshold:
                yield res_no

def get_inner_residues(prot_pdb, list_of_res):
    res_nos={i for i in list_of_res}
    for atom in parser.parse_file(prot_pdb):
        if atom['resSeq'] in res_nos:
            yield parser.pdb_format(atom)

def separate_by_res_no(pdb_file, subunit):
    separated_residues = {}
    for atom in parser.filter_line(parser.parse_file(pdb_file), get_subunit(subunit)):
        if not (atom['chainID'],atom['resSeq']) in separated_residues:
            separated_residues[(atom['chainID'], atom['resSeq'])] = []
        separated_residues[(atom['chainID'], atom['resSeq'])].append(atom)
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

def get_distances(pdb_file, subunit):
    separated_residues = separate_by_res_no(pdb_file, subunit)
    for res1, res2 in make_pairs(separated_residues):
        for distance in get_pairwise_distances(separated_residues[res1], separated_residues[res2]):
            yield distance

def get_distances_from_core(pdb_file, subunit):
    os.popen("~/protein_lab/third-party/binaries/reduce {0} > temp1.pdb".format(pdb_file))
    os.popen("/home/bolt/protein_lab/third-party/binaries/naccess2.1.1/naccess temp1.pdb")
    temp_output=open("temp2.pdb", "w")
    for l in get_inner_residues("temp1.pdb", get_buried_residues("temp1.rsa")):
        temp_output.write(l+"\n")
    temp_output.close()
    for distance in get_distances("temp2.pdb", subunit):
        print(distance)
    os.popen("rm temp1.pdb")
    os.popen("rm temp2.pdb")
    os.popen("rm temp1.rsa")
    os.popen("rm temp1.asa")
    os.popen("rm temp1.log")

if __name__ == "__main__":
    from sys import argv
    get_distances_from_core(argv[-2], argv[-1])
