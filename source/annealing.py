#!/usr/bin/env python

import random as rd
import residue_mutator as mut


def get_neighbour_state(protein_pdb):
    """
    Placeholder function
    """

def score(state):
    """
    Placeholder function
    """
    pass

def acceptor(old_state, new_state, temperature):
    sigma = random.random()
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
