#!/usr/bin/env python

"""
Check the ramachandran function for any errors.
"""

import residue_mutator as rm
import numpy as np

class TerminalResidueError(Exception):
    pass

def get_phi(n, backbone_dict):
    if n<=1:
        raise TerminalResidueError
    point1 = np.matrix(backbone_dict[n-1]['C'])
    point2 = np.matrix(backbone_dict[n]['N'])
    point3 = np.matrix(backbone_dict[n]['CA'])
    point4 = np.matrix(backbone_dict[n]['C'])
    return rm.get_angle((point1-point2), (point2-point3), (point4-point3))

def get_psi(n, backbone_dict):
    if n>=len(backbone_dict):
        raise TerminalResidueError
    point1 = np.matrix(backbone_dict[n]['N'])
    point2 = np.matrix(backbone_dict[n]['CA'])
    point3 = np.matrix(backbone_dict[n]['C'])
    point4 = np.matrix(backbone_dict[n+1]['N'])
    return rm.get_angle((point1-point2), (point2-point3), (point4-point3))

def ramachandran(n, backbone_dict):
    args = n, backbone_dict
    return get_phi(*args),get_psi(*args)
