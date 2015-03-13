#!/usr/bin/env python

import numpy as np
import kabsch as kb
import pdb_parser as parser
import math

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
    res_pdb=argv[1]
    chis =[int(i) for i in argv[2:]]
    for l in rotate_by_chi(chis, parser.parse_file(res_pdb)):
        print(parser.pdb_format(l))
