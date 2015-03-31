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

import random as rd
import residue_mutator as mut
import get_atom_distances as gad
import pdb_parser as parser
from __future__ import print_function

allowed_residues = ["VAL", "LEU", "ILE", "MET", "ALA"]
residue_path = "~/protein_lab/database/residues/{0}_naiveH.pdb"

# Helper functions

def get_current_residue(protein_pdb, res_no):
    for line in parser.parse_file(protein_pdb):
        if line['resSeq']==res_no:
            return line['resName']

# End of helper functions

def get_neighbour_state(protein_pdb):
    """
    Placeholder function
    """
    residue_to_mutate=rd.choice(gad.get_buried_residues_pdb(protein_pdb))
    current_residue = get_current_residue(protein_pdb, residue_to_mutate)
    new_residue = rd.choice([i for i in allowed_choices if i!=current_residue])
    protein_lines = list(parser.parse_file(protein_pdb))
    residue_lines = list(parser.parse_file(residue_path.format(new_residue)))
    new_protein = mut.mutate(protein_lines, residue_to_mutate, residue_lines)


def score(state):
    """
    Placeholder function
    """
    pass

def acceptor(old_state, new_state, temperature):
    sigma = rd.random()
    if sigma<=temperature:
        return new_state
    else:
        if score(new_state) <= score(old_state):
            return new_state
        else:
            return old_state

def temperature_func(i, iterations):
    """
    Placeholder function
    """
    pass

def anneal(start_state, iterations, temperature_function):
    current_state=start_state
    for i in range(iterations):
        current_temperature = temperature_function(i, iterations)
        new_state = get_neighbour_state(current_state)
        current_state = acceptor(current_state, new_state, current_temperature)
    return current_state
