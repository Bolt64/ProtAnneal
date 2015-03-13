#!/usr/bin/env python

import numpy as np
import kabsch as kb
import pdb_parser as parser
import math

# This script may not be 100% correct. It has not been tested
# rigourously enough. It's probably correct, but don't hold me
# to my word. The residue mutator function works fine, but the other
# sidechain rotating function still needs to be tested.

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

def mutate_to_file(protein_pdb, res_no, res_pdb, output_pdb):
    new_protein=mutate(protein_pdb, res_no, res_pdb)
    with open(output_pdb, "a") as output:
        for line in new_protein:
            output.write(parser.pdb_format(line))
            output.write("\n")

def get_rotation_matrix(angle, axis):
    """
    angle: angle in degrees
    axis: numpy.matrix
    """
    angle = -(angle/180.0)*math.pi
    matrix = []
    axis = axis/(np.linalg.norm(axis))
    axis = axis.tolist()
    x = axis[0][0]
    y = axis[0][1]
    z = axis[0][2]
    cos = math.cos(angle)
    sin = math.sin(angle)
    matrix.append([
                  cos+((x**2)*(1-cos)),
                  x*y*(1-cos)-(z*sin),
                  x*z*(1-cos)+(y*sin)
                  ])
    matrix.append([
                  y*x*(1-cos)+(z*sin),
                  cos+(y**2)*(1-cos),
                  y*z*(1-cos)-(x*sin)
                  ])
    matrix.append([
                  z*x*(1-cos)-(y*sin),
                  z*y*(1-cos)+(x*sin),
                  cos+(z**2)*(1-cos)
                  ])
    return np.matrix(matrix)

def rotate(set_of_points, pivot, rotation_matrix):
    """
    set_of_points: numpy.matrix
    pivot: tuple(x, y, z)
    rotation_matrix: numpy.matrix
    """
    return kb.translate(kb.translate(set_of_points, [-i for i in pivot])*rotation_matrix, pivot)

def remove_hydrogens(list_of_lines):
    """
    Removes hydrogen from the pdb file.
    To add back the hydrogens, run the reduce program on the file.
    """
    return (line for line in list_of_lines if line['element']!=" H")

def rotate_by_chi(chis, list_of_lines):
    """
    chis: tuple of angles
    list_of_lines: A list of lines that correspond to a residue pdb
    This function takes a list of lines, all corresponding to one amino
    acid, and rotates it by the chis supplied in the function
    """
    atoms=[]
    for line in remove_hydrogens(list_of_lines):
        if line['name']==' C  ':
            C_line=line
        else:
            atoms.append(line)
    for index,chi in enumerate(chis):
        #print(atoms[index]['name'], atoms[index+1]['name'], atoms[index+2]['name'], atoms[index+3]['name'])
        axis=np.matrix([atoms[index+2]['x']-atoms[index+1]['x'], atoms[index+2]['y']-atoms[index+1]['y'], atoms[index+2]['z']-atoms[index+1]['z']])
        pivot=[atoms[2+index][i] for i in 'xyz']
        set_of_points=np.matrix([[atom['x'], atom['y'], atom['z']] for atom in atoms[index+3:]])
        rotation_matrix=get_rotation_matrix(chi, axis)
        for offset,new_location in enumerate(rotate(set_of_points, pivot, rotation_matrix)):
            x,y,z=new_location.tolist()[0]
            atoms[index+3+offset]['x']=x
            atoms[index+3+offset]['y']=y
            atoms[index+3+offset]['z']=z
    return atoms[0:2]+[C_line]+atoms[2:]


if __name__=="__main__":
    from sys import argv
    protein_pdb=argv[-3]
    res_no=int(argv[-2])
    res_pdb=argv[-1]
    for l in mutate(protein_pdb, res_no, res_pdb):
        print(parser.pdb_format(l))
