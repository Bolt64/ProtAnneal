#!/usr/bin/env python

"""
For the get_neighbour_state, the neighbour state is gotten by mutating one of the inner core residues to one of valine, leucine, isoleucine, alanine or methionine.
For starters, try just putting in the residues, and then try to use the rotamer library to align the sidechains for max benefit.

Step 1: Take a reduced pdb file with only one subunit.
Step 2: Get the core, residue numbers only, using the naccess tool.
Step 3: Randomly select one of the core residues, but only among the ones which
        belong to the allowed residues.
Step 4: Mutate it to one of the other residues.
Step 5: Use the rotamer library to align it in an allowed conformation.
Step 6: Return this new residue.
"""

from __future__ import print_function
import random as rd
import residue_mutator as mut
import get_atom_distances as gad
import pdb_parser as parser
import zdope_score as energy
import ramachandran_angles as ra
import dunbrack_parser as dp
from config import *

#allowed_residues = ["VAL", "LEU", "ILE", "MET", "ALA"]
#residue_path = "/home/bolt/protein_lab/database/residues/{0}_naiveH.pdb"
#dunbrack_file = "/home/bolt/protein_lab/other/dunbrack2.csv"

# Helper functions

def get_current_residue(protein_pdb, res_no):
    for line in parser.parse_file(protein_pdb):
        if line['resSeq']==res_no:
            return line['resName']

def memoize(func):
    cache = {}
    def inner(*args):
        if not args in cache:
            cache[args]=func(*args)
        return cache[args]
    return inner

def round_to_tens(num):
    return int(10.0*round(num/10.0))

def weighted_choice(choices):
    total = sum(w for w, c in choices)
    r = rd.uniform(0, total)
    upto = 0
    for w, c in choices:
        if upto + w > r:
            return c
        upto += w
    assert False, "Shouldn't get here"

# End of helper functions

def get_neighbour_state(protein_pdb, name, dunbrack_library, buried_residues):
    """
    Placeholder function
    This does not rotate the side chains. Only replaces the
    residues.
    """
    new_file=open(name, "a")
    residue_to_mutate=rd.choice(buried_residues)
    current_residue = get_current_residue(protein_pdb, residue_to_mutate)
    new_residue = rd.choice([i for i in allowed_residues if i!=current_residue])
    protein_lines = list(parser.parse_file(protein_pdb))
    residue_lines = list(parser.parse_file(residue_path.format(new_residue)))

    backbone_dictionary = mut.get_backbone_dictionary(mut.get_backbone(protein_lines))
    phi = round_to_tens(ra.get_phi(residue_to_mutate, backbone_dictionary))
    psi = round_to_tens(ra.get_psi(residue_to_mutate, backbone_dictionary))
    main_key = new_residue, phi, psi
    choices = dunbrack_library[main_key]
    weighted_choices = [(i, choices[i]) for i in choices]
    choice = weighted_choice(weighted_choices)
    new_lines = mut.rotate_to_chis(choice, residue_lines)

    new_protein = mut.mutate(protein_lines, residue_to_mutate, new_lines)
    for l in new_protein:
        line = parser.pdb_format(l)
        print(line, file=new_file)
    new_file.close()

    return name

@memoize
def score(state):
    return energy.get_dope_score(state)

def acceptor(old_state, new_state, temperature):
    sigma = rd.random()
    if sigma<=temperature:
        return new_state
    else:
        if score(new_state) <= score(old_state):
            return new_state
        else:
            return old_state

def linear_temperature_func(i, iterations):
    """
    Placeholder function
    """
    return 1-float(i)/iterations

def exp_temperature_func(i, iterations):
    return 2**(-10.0*i/iterations)

def sigmoid_temperature_func(i, iterations):
    return 1-1.0/(1+2.0**(-20*(float(i)/iterations-0.2)))

def anneal(start_state, iterations, temperature_function, prefix=""):
    scores=[]
    dunbrack_library = dp.parse_dunbrack(dunbrack_file)
    current_state=start_state
    buried_residues = gad.get_buried_residues_pdb(start_state)
    for i in range(iterations):
        scores.append(score(current_state))
        current_temperature = temperature_function(i, iterations)
        new_state = get_neighbour_state(current_state, prefix+str(i)+".pdb", dunbrack_library, buried_residues)
        current_state = acceptor(current_state, new_state, current_temperature)
        print(current_state)
    return current_state, scores
