#!/usr/bin/env python

import random
from math import sin,cos
import numpy as np
import pdb_parser as parser
import kabsch as kb


def rotate_x(theta):
    return np.matrix([
            [1, 0, 0],
            [0, cos(theta), -sin(theta)],
            [0, sin(theta), cos(theta)]
            ])


def rotate_y(theta):
    return np.matrix([
            [cos(theta), 0, sin(theta)],
            [0, 1, 0],
            [-sin(theta), 0, cos(theta)]
            ])

def rotate_z(theta):
    return np.matrix([
            [cos(theta), -sin(theta), 0],
            [sin(theta), cos(theta), 0],
            [0, 0, 1]
            ])

def randomize_protein(protein_pdb):
    rotation=rotate_x(random.randint(0,100))*rotate_y(random.randint(0,100))*rotate_z(random.randint(0,100))
    translation=[random.randint(-100,100), random.randint(-100,100), random.randint(-100,100)]
    for l in parser.parse_file(protein_pdb):
        loc=np.matrix([[l['x'],l['y'],l['z']]])
        new_loc=kb.translate(loc*rotation, translation)
        x,y,z=new_loc.tolist()[0]
        l['x']=x
        l['y']=y
        l['z']=z
        yield l

if __name__=="__main__":
    from sys import argv
    for l in randomize_protein(argv[-1]):
        print(parser.pdb_format(l))
