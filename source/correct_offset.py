#!/usr/bin/env python

import pdb_parser as parser

def correct_offset(list_of_lines):
    list_of_lines = list(list_of_lines)
    offset = list_of_lines[0]['resSeq']
    new_lines = []
    for line in list_of_lines:
        line['resSeq']-=(offset-1)
        new_lines.append(line)
    return new_lines

if __name__=="__main__":
    from sys import argv
    prot_pdb = argv[-1]
    for l in correct_offset(parser.parse_file(prot_pdb)):
        print(parser.pdb_format(l))
