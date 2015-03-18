#!/usr/bin/python
#Praveen created this script to evaluate z-dope for a given model
import sys

#Ensuring whether the right number of arguments are being passed with the script
if len(sys.argv) != 2:
   print 'Usage: calculate_zdope.py <pdb file>'
   sys.exit(1)

#Ensuring that PDB file exists at the right place 
from os import path, access, R_OK

PDB_FILE = sys.argv[1]
if path.exists(PDB_FILE) and path.isfile(PDB_FILE) and access(PDB_FILE, R_OK):
   print "Following are the results for your PDB -", sys.argv[1]
else:
   print "Your PDB file does not exist. Please check the path of your PDB file."
   sys.exit(1)

#Exercising modeller library to get the z-dope score for the alignment

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = complete_pdb(env, sys.argv[1])

zscore = mdl.assess_normalized_dope()
