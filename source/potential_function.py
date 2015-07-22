#usr/bin/python

"""
This module uses a potential function on a distance file
to get total potential of a protein
"""

import matplotlib.pyplot as plt

# The currently used potential function
lj_potential = lambda distance, rm: (rm/distance)**12 - 2*(rm/distance)**6

def parse_distance_file(filename):
    for line in open(filename):
        yield float(line.strip())

def get_total_potential(distance_file, potential_function):
    """
    Gets the total potential given a distance file and
    a potential function.
    """
    return sum(potential_function(distance) for distance in parse_distance_file(distance_file))

def vary_across_rm(distance_file, rm_range, step_size, potential_function):
    """
    Varies rm over rm_range with step_size and yields
    potential.
    """
    lower,upper=rm_range
    current_rm=lower
    while current_rm < upper:
        current_potential = lambda distance : potential_function(distance, current_rm)
        yield current_rm, get_total_potential(distance_file, current_potential)
        current_rm += step_size

def get_min_rm(distance_file, rm_range, step_size, potential_function):
    """
    Gets the rm at which potential energy is minimized.
    """
    return min((i[1],i[0]) for i in vary_across_rm(distance_file, rm_range, step_size, potential_function))[1]


if __name__=="__main__":
    from sys import argv
    distance_file=argv[-3]
    lower=float(argv[-2])
    upper=float(argv[-1])
    step_size=0.1
    xy=[i for i in vary_across_rm(distance_file, (lower, upper), step_size, lj_potential)]
    x=[i[0] for i in xy]
    y=[i[1] for i in xy]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Total lennard-jones potential vs VdW radius")
    ax.set_xlabel("Van der Waals radius (Angstrom)")
    ax.set_ylabel("Total lennard-jones potential")
    ax.plot(x,y,'ro')
    plt.savefig("/tmp/output")
