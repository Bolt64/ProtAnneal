import os
import csv
import sys
import copy
import math
import time
import numpy
import string
import random
import shutil
import getopt
import subprocess
import matplotlib.pyplot as plt

##########
## MAIN ##
##########
def main():
	protein = open(sys.argv[1])
	protein_data = protein.readlines()
	z = 0
	for k in range(0, len(protein_data)):
		if (protein_data[k][0:6] == 'ATOM  '):
			b = int(protein_data[k][22:26])
			if (b>z):
				z = b
	c1_values=[]
	for i,_ in enumerate(protein_data):
		if protein_data[i][0:6] == 'ATOM  ':
			a = int(protein_data[i][22:26])
			if (a>z-2):
				break
			if (protein_data[i][12:16] == ' C  ' and int(protein_data[i][22:26]) == a):
				c1_values.append([float(protein_data[i][30:38]), float(protein_data[i][38:46]), float(protein_data[i][46:54])])
	n2_values=[]
	ca2_values=[]
	c2_values=[]
	for i,_ in enumerate(protein_data):
		if protein_data[i][0:6] == 'ATOM  ':
			a = int(protein_data[i][22:26])
			if (a>1 and a<z):
				if (protein_data[i][12:16] == ' N  ' and int(protein_data[i][22:26]) == a):
					n2_values.append([float(protein_data[i][30:38]), float(protein_data[i][38:46]), float(protein_data[i][46:54])])
				elif (protein_data[i][12:16] == ' CA ' and int(protein_data[i][22:26]) == a):
					ca2_values.append([float(protein_data[i][30:38]), float(protein_data[i][38:46]), float(protein_data[i][46:54])])
				elif (protein_data[i][12:16] == ' C  ' and int(protein_data[i][22:26]) == a):
					c2_values.append([float(protein_data[i][30:38]), float(protein_data[i][38:46]), float(protein_data[i][46:54])])
	n3_values=[]
	for i,_ in enumerate(protein_data):
		if protein_data[i][0:6] == 'ATOM  ':
			a = int(protein_data[i][22:26])
			if (a>2):
				if (protein_data[i][12:16] == ' N  ' and int(protein_data[i][22:26]) == a):
					n3_values.append([float(protein_data[i][30:38]), float(protein_data[i][38:46]), float(protein_data[i][46:54])])
	phi = []
	psi = []
	for i in range(0, z-2):
		vec_c1n2 = numpy.matrix([c1_values[i], n2_values[i]])
		vec_n2ca2 = numpy.matrix([n2_values[i], ca2_values[i]])
		vec_ca2c2 = numpy.matrix([ca2_values[i], c2_values[i]])
		vec_c2n3 = numpy.matrix([c2_values[i], n3_values[i]])
		temp = ramachandran(vec_c1n2, vec_n2ca2, vec_ca2c2, vec_c2n3)
		phi.append(temp[0])
		psi.append(temp[1])
	plt.title('Ramachandran Plot of '+str(sys.argv[1]))
	plt.xlabel('Phi ---->')
	plt.ylabel('Psi ---->')
	plt.plot(phi, psi, 'ro')
	plt.savefig('Rama-'+str(sys.argv[1])+'.jpeg')


###############
## FUNCTIONS ##
###############

#this function gives the angle betwen 2 vectors:
def dot_angle(vec_1, vec_2):
	a = (vec_1[1,0]-vec_1[0,0])*(vec_2[1,0]-vec_2[0,0]) + (vec_1[1,1]-vec_1[0,1])*(vec_2[1,1]-vec_2[0,1]) + (vec_1[1,2]-vec_1[0,2])*(vec_2[1,2]-vec_2[0,2])
	b = math.sqrt((vec_1[1,0]-vec_1[0,0])**2 + (vec_1[1,1]-vec_1[0,1])**2 + (vec_1[1,2]-vec_1[0,2])**2)
	c = math.sqrt((vec_2[1,0]-vec_2[0,0])**2 + (vec_2[1,1]-vec_2[0,1])**2 + (vec_2[1,2]-vec_2[0,2])**2)
	cos = a/(b*c)
	if (cos>1):
		cos = 1
	elif (cos<-1):
		cos = -1
	angle_rad = math.acos(cos)
	angle_deg = (angle_rad/math.pi)*180
	return angle_deg

#This function gives a vector that is the result of a cross product between 2 vectors:
def cross_vector(vec_1, vec_2):
	a1 = vec_1[1,0]-vec_1[0,0]
	b1 = vec_1[1,1]-vec_1[0,1]
	c1 = vec_1[1,2]-vec_1[0,2]
	a2 = vec_2[1,0]-vec_2[0,0]
	b2 = vec_2[1,1]-vec_2[0,1]
	c2 = vec_2[1,2]-vec_2[0,2]
	x = c2*b1 - b2*c1
	y = a2*c1 - c2*a1
	z = b2*a1 - a2*b1
	OP_vector = numpy.matrix([[0,0,0], [x, y, z]])
	return OP_vector    

#This function gives the angle between two planes, if three vectors are given as input, one in the first plane, the other in the second plane, and the third one common to both the planes
def dihedral(vec_1, vec_2, vec_3):
	a = cross_vector(vec_1, vec_2)
	b = cross_vector(vec_2, vec_3)
	angle = dot_angle(a, b)
	c = cross_vector(a, -b)
	d = dot_angle(c, vec_2)
	if (d < 90):
		angle = -angle
	return angle

def ramachandran(vec_1, vec_2, vec_3, vec_4):
	phi = dihedral(vec_1, vec_2, vec_3)
	psi = dihedral(vec_2, vec_3, vec_4)
	return [phi, psi]
	#print ("Residue no. = "+str(i)+"; phi = "+str(phi)+"; psi = "+str(psi))

if __name__ == '__main__':
    main()

