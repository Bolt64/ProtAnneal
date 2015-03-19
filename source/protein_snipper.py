#!/usr/bin/env python

"""
Snips off ends from either ends of a protein unit.
"""

import pdb_parser as parser

def snip(prot_pdb, threshold):
    """
    Takes in a prot_pdb file and threshold values
    and snips off the ends.
    """
    lower,upper=threshold
    for l in parser.parse_file(prot_pdb):
        if lower < l['resSeq'] < upper:
            yield l

if __name__ == "__main__":
    from sys import argv
    lower = int(argv[-2])
    upper = int(argv[-1])
    prot_pdb = argv[-3]
    for l in snip(prot_pdb, (lower, upper)):
        print(parser.pdb_format(l))
