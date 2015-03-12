#usr/bin/python!

import numpy
import sys
import math

def rotation_matrix(angle, axis):
	angle = -(angle/180)*math.pi
	matrix = []
	axis = axis/(numpy.linalg.norm(axis))
	axis = axis.tolist()
	#print(axis)
	x = axis[0][0]
	y = axis[0][1]
	z = axis[0][2]
	print(x, y, z)
	cos = math.cos(angle)
	sin = math.sin(angle)
	matrix.append([cos+((x**2)*(1-cos)), (x*y*(1-cos))-(z*sin), (x*z*(1-cos))+(y*sin)])
	matrix.append([(y*x*(1-cos))+(z*sin), cos+((y**2)*(1-cos)), (y*z*(1-cos))-(x*sin)])
	matrix.append([(z*x*(1-cos))-(y*sin), (z*y*(1-cos))+(x*sin), cos+((z**2)*(1-cos))])
	matrix = numpy.matrix(matrix)
	return matrix

if __name__=="__main__":
	a = float(sys.argv[1])
	b = numpy.matrix(sys.argv[2])
	c = rotation_matrix(a, b)
	print(c)
