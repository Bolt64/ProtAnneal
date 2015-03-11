#!/usr/bin/env python

import numpy as np
import kabsch as kb
import pdb_parser as parser

# This script may not be 100% correct. It has not been tested
# rigourously enough. Unless tested further, do not base any scientific
# data on the results of this script. You have been warned.

def get_backbone(pdb_file):
    """
    Gets the N, CA, C backbone of the pdb file
    """
    backbone=[]
    for line in parser.parse_file(pdb_file):
        if line['ATOM']=="ATOM  ":
            if line['name'] in (' N  ', ' CA ', ' C  '):
                backbone.append([line['x'],line['y'],line['z']])
    return np.matrix(backbone)

def mutate(protein_pdb, res_no, res_pdb):
    """
    It replaces the residue on protein at res_no
    position with the new residue, res_pdb.
    This function may be dodgy. DO NOT USE without further testing.
    You have been warned.
    """
    backbone_original=[]
    for line in parser.parse_file(protein_pdb):
        if line['ATOM']=="ATOM  " and line['resSeq']==res_no:
            if line['name'] in (' N  ', ' CA ', ' C  '):
                backbone_original.append([line['x'], line['y'], line['z']])
    backbone_replacement=get_backbone(res_pdb).tolist()
# Should the arguments be switched around?
    rotation=kb.kabsch(np.matrix(backbone_original), np.matrix(backbone_replacement))
    centroid_original=kb.get_centroid(np.matrix(backbone_original))
    centroid_replacement=kb.get_centroid(np.matrix(backbone_replacement))
    inverted_centroid=[-1*i for i in centroid_replacement]
    new_protein=[]
    amino_acid_inserted=False
    for l in parser.parse_file(protein_pdb):
        if l['resSeq']!=res_no:
            new_protein.append(l)
        if l['resSeq']==res_no and not amino_acid_inserted:
            for j in parser.parse_file(res_pdb):
                location=np.matrix([[j['x'],j['y'], j['z']]])
                new_loc=kb.translate((kb.translate(location, inverted_centroid))*rotation, centroid_original)
                x,y,z=new_loc.tolist()[0]
                j['x']=x
                j['y']=y
                j['z']=z
                j['resSeq']=res_no
                new_protein.append(j)
            amino_acid_inserted=True
    return new_protein

if __name__=="__main__":
    from sys import argv
    protein_pdb=argv[-3]
    res_no=int(argv[-2])
    res_pdb=argv[-1]
    for l in mutate(protein_pdb, res_no, res_pdb):
        print(parser.pdb_format(l))
