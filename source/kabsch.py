#!/usr/bin/env python

import numpy as np
import math

"""
Kabsch algorithm and rmsd calculator based on the wiki article.
Kabsch and rmsd functions take in two numpy arrays. Please be type
conscious while passing parameters since strict error checking is not
implemented.
"""

class ImproperSizeError(Exception):
    pass

def sign(num):
    """
    Determines the sign of a number
    """
    if num >= 0:
        return 1
    else:
        return -1

def average(iterable):
    """
    Finds average of the members in the iterable. Wrote because of no
    statistics module in standard library
    """
    number = 0
    total = 0
    for i in iterable:
        number += 1
        total += i
    return total/float(number)

def get_array_i(array, i):
    """
    Built to get ith element of a numpy array
    """
    return array.tolist()[0][i]

def translate(set_of_points, translation_matrix):
    points = set_of_points.tolist()
    x, y, z = translation_matrix
    translated= [[i[0]+x, i[1]+y, i[2]+z] for i in points]
    return np.matrix(translated)

def rmsd(set1, set2):
    """
    set1: numpy.matrix
    set2:numpy.matrix
    Built to get rmsd. Speaks for itself.
    """
    if set1.shape != set2.shape:
        raise ImproperSizeError
    distance_squared_sum= []
    for i in range(len(set1)):
        x1,y1,z1 = set1[i].tolist()[0]
        x2,y2,z2 = set2[i].tolist()[0]
        distance_squared_sum.append((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    return math.sqrt(average(distance_squared_sum))

def get_centroid(set_of_points):
    """
    set_of_points: numpy.matrix
    """
    centroid_x = average((get_array_i(i, 0) for i in set_of_points))
    centroid_y = average((get_array_i(i, 1) for i in set_of_points))
    centroid_z = average((get_array_i(i, 2) for i in set_of_points))
    return centroid_x, centroid_y, centroid_z

def centre(set_of_points):
    """
    set_of_points: numpy.matrix
    """
    x,y,z = get_centroid(set_of_points)
    return np.matrix([[get_array_i(i, 0)-x, get_array_i(i, 1)-y, get_array_i(i, 2)-z] for i in set_of_points])

def centred_covariance(set1, set2):
    """
    set1: numpy.matrix
    set2:numpy.matrix
    """
    return set1.transpose()*set2

def optimal_rotation_matrix(A):
    """
    A: numpy.matrix
    """
    v,s,w_trans = np.linalg.svd(A)
    d = sign(np.linalg.det(w_trans.transpose()*v.transpose()))
    temp_matrix = np.matrix([[1,0,0],[0,1,0],[0,0,d]])
    return w_trans.transpose()*temp_matrix*v.transpose()

def kabsch(set1, set2):
    """
    set1: numpy.matrix
    set2:numpy.matrix
    """
    if set1.shape != set2.shape:
        print(set1)
        print(set2)
        raise ImproperSizeError
    centred_1 = centre(set1)
    centred_2 = centre(set2)
    A = centred_covariance(centred_1, centred_2)
    return optimal_rotation_matrix(A)

def optimal_rmsd(set1, set2):
    """
    set1: numpy.matrix
    set2:numpy.matrix
    """
    centred_1 = centre(set1)
    centred_2 = centre(set2)
    return rmsd(centred_1, centred_2*kabsch(set1, set2))
