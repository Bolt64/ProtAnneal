#!/usr/bin/env python

"""
The main script that does it all
"""

from __future__ import print_function
import annealing as an
import os
import create_primer as primer
import matplotlib.pyplot as plt
from config import *


def main(original_pdb, iterations, temperature_function, result_dir):
    os.popen("mkdir {0}".format(output_directory))
    prefix = output_directory+"/"
    new_pdb, scores = an.anneal(original_pdb, iterations, temperature_function, prefix)
    os.popen("mkdir {0}".format(result_dir))
    os.popen("cp {0} {1}".format(original_pdb, result_dir+"/original.pdb"))
    os.popen("cp {0} {1}".format(new_pdb, result_dir+"/new.pdb"))
    os.popen("rm -r {0}".format(output_directory))
    x = list(range(iterations))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Energy of the protein vs iterations")
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Zdope score')
    ax.plot(x, scores)
    plt.savefig("{0}/result.png".format(result_dir))
    score_file = open("{0}/scores.txt".format(result_dir), "a")
    for s in scores:
        print(s, file=score_file)
    score_file.close()

if __name__=="__main__":
    from sys import argv
    try:
        orig_pdb = argv[1]
        iterations = int(argv[2])
        main(orig_pdb, iterations, an.sigmoid_temperature_func, result_directory)
    except (ValueError, IndexError):
        print("Invalid Syntax:")
        print("\tpython main.py <pdb_file> <number of iterations>")
