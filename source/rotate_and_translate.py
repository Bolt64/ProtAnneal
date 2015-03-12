#!/usr/bin/env python

import numpy as np
import kabsch as kb
import pdb_parser as parser
import math

def rotation_matrix(angle, axis):
    """
    angle: float in degress
    axis: numpy.matrix
    """
	angle = -(angle/180)*math.pi
	matrix = []
	axis = axis/(np.linalg.norm(axis))
	axis = axis.tolist()
	x = axis[0][0]
	y = axis[0][1]
	z = axis[0][2]
	cos = math.cos(angle)
	sin = math.sin(angle)
	matrix.append([cos+((x**2)*(1-cos)), (x*y*(1-cos))-(z*sin), (x*z*(1-cos))+(y*sin)])
	matrix.append([(y*x*(1-cos))+(z*sin), cos+((y**2)*(1-cos)), (y*z*(1-cos))-(x*sin)])
	matrix.append([(z*x*(1-cos))-(y*sin), (z*y*(1-cos))+(x*sin), cos+((z**2)*(1-cos))])
	return np.matrix(matrix)

def rotate(set_of_points, pivot, rotation_matrix):
    """
    set_of_points: numpy.matrix
    pivot: tuple(x, y, z)
    rotation_matrix: numpy.matrix
    """
    return kb.translate(kb.translate(set_of_points, [-i for i in pivot])*rotation_matrix, pivot)

def remove_hydrogens(list_of_lines):
    return (line for line in list_of_lines if line['element']!=" H")

def rotate_by_chi(chis, res_pdb):
    atoms=[]
    for line in remove_hydrogens(parser_parse_file(res_pdb)):
        if l['name']==' C  ':
            C_line=l
        else:
            atoms.append(l)
    for index,chi in enumerate(chis):
        axis=np.matrix([atoms[index+2]['x']-atoms[index+1]['x'], atoms[index+2]['y']-atoms[index+1]['y'], atoms[index+2]['z']-atoms[index+1]['z']])
        pivot=[atoms[2][i] for i in 'xyz']
        set_of_points=np.matrix([[atom['x'], atom['y'], atom['z']] for atom in atoms[i+3:]])
        rotation_matrix
