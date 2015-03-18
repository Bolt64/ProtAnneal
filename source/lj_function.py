#usr/bin/python

import math
from sys import argv
import matplotlib.pyplot as plt

def lj_potential(distance, rm):
    return math.pow((rm/distance), 12) - 2*math.pow((rm/distance), 6)

def extract(filename):
    for line in open(filename):
        yield float(line.strip())

def get_total_potential(distance_file, rm, potential_function):
    return sum(potential_function(distance, rm) for distance in extract(distance_file))

def vary_across_rm(distance_file, rm_range, step_size, potential_function):
    lower,upper=rm_range
    current_rm=lower
    while current_rm < upper:
        yield current_rm, get_total_potential(distance_file, current_rm, potential_function)
        current_rm += step_size

def get_min_rm(distance_file, rm_range, step_size, potential_function):
    return min((i[1],i[0]) for i in vary_across_rm(distance_file, rm_range, step_size, potential_function))[1]


if __name__=="__main__":
    distance_file=argv[-3]
    lower=float(argv[-2])
    upper=float(argv[-1])
    step_size=0.1
    xy=[i for i in vary_across_rm(distance_file, (lower, upper), step_size, lj_potential)]
    x=[i[0] for i in xy]
    y=[i[1] for i in xy]
    plt.plot(x,y,'ro')
    plt.savefig("output")
