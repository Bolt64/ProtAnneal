#!/usr/bin/env python

import pdb_parser as parser
from residue_mutator import subunit_predicate

def get_subunit(pdb_file, subunit):
    predicate=subunit_predicate(subunit)
    return parser.filter_line(parser.parse_file(pdb_file), predicate)

if __name__=="__main__":
    from sys import argv
    subunit=argv[-2]
    pdb_file=argv[-1]
    for l in get_subunit(pdb_file, subunit):
        print(parser.pdb_format(l))
