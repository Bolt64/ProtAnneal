#usr/bin/python

import math
import sys
import matplotlib.pyplot as plt

def lj_potential(dist, rm):
    lj = 0
    for i in range(0, len(dist)):
        lj += math.pow((rm/dist[i]), 12) - 2*math.pow((rm/dist[i]), 6)
    return lj

def extract(filename):
    distances = []
    for line in open(filename):
        distances.append(float(line.strip()))
    return distances

if __name__ == '__main__':
    pot = []
    rm = []
    d = extract(sys.argv[1])
    for i in range (0, 100):
        rm.append(i*0.1)
        v = lj_potential(d, i*0.5)
        pot.append(v)
    plt.title('Lennard-Johnes Plot of '+str(sys.argv[1]))
    plt.xlabel('Rm ---->')
    plt.ylabel('LJ-Potential ---->')
    plt.plot(pot, rm, 'ro')
    plt.savefig('LJ-potential-'+str(sys.argv[1])+'.jpeg')
