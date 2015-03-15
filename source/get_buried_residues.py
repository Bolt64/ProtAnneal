#!/usr/bin/env python

import pdb_parser as parser

def get_buried_residues(rsa_file, threshold=30):
    for line in open(rsa_file):
        if line.startswith("RES"):
            temp=line.strip().split()
            res_no=int(temp[3])
            rel_accessibility=float((temp[7]))
            if rel_accessibility<=threshold:
                yield res_no

def generate_residue_numbers(rsa_file, threshold=30):
    with open(rsa_file[:-3]+"rnos", "a") as output:
        for num in get_buried_residues(rsa_file, threshold):
            output.write(str(num)+'\n')

def get_inner_residues(prot_pdb, list_of_res):
    res_nos={i for i in list_of_res}
    for atom in parser.parse_file(prot_pdb):
        if atom['resSeq'] in res_nos:
            yield parser.pdb_format(atom)

if __name__=="__main__":
    from sys import argv
    rsa_file=argv[-1]
    prot_pdb=argv[-2]
    for l in get_inner_residues(prot_pdb, get_buried_residues(rsa_file)):
        print(l)
