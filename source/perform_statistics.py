#!/usr/bin/env python

def parse_distance_file(distance_file):
    for line in open(distance_file):
        yield float(line.strip())

def put_in_bins(list_of_distances, bin_size):
    bins={}
    maximum=0
    for distance in list_of_distances:
        lower=int(distance/bin_size)
        upper=lower+1
        if not (lower, upper) in bins:
            bins[(lower, upper)]=0
        bins[(lower, upper)]+=1
        if lower>maximum:
            maximum=lower
    for lower in xrange(maximum+1):
        if not (lower, lower+1) in bins:
            bins[(lower, lower+1)]=0
    return {(i[0]*bin_size, i[1]*bin_size):bins[i]/float(i[1]*bin_size)**2 for i in bins}
