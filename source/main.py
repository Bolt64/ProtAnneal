#!/usr/bin/env python

"""
The main script that does it all
"""

import annealing as an
import os
import matplotlib.pyplot as plt
from config import *


def main(original_pdb, iterations, temperature_function, result_dir):
    os.popen("mkdir {0}".format(output_directory))
    prefix = output_directory+"/"
    new_pdb, scores = an.anneal(original_pdb, iterations, temperature_function, prefix)
    os.popen("mkdir {0}".format(result_dir))
    os.popen("cp {0} {1}".format(original_pdb, result_dir+"/original.pdb"))
    os.popen("cp {0} {1}".format(new_pdb, result_dir+"/new.pdb"))
    #os.popen("rm -r {0}".format(output_directory))
    x = list(range(iterations))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Energy of the protein vs iterations")
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Zdope score')
    ax.plot(x, scores)
    plt.savefig("{0}/result.png".format(result_dir))

