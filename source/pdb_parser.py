#!/usr/bin/env python

"""
A parser for pdb files as well as a writer.
It also has a a nice filter function.
"""

def parse_file(pdb_file):
    """
    Parses the file and returns as dictionary
    with appropriate fields.
    """
    for line in open(pdb_file):
        temp = line.strip().ljust(80, " ")
        if temp[0:6] == "ATOM  ":
            line_output = {}
            line_output['ATOM'] = temp[0:6]
            line_output['serial'] = int(temp[6:11])
            line_output['name'] = temp[12:16]
            line_output['altLoc'] = temp[16]
            line_output['resName'] = temp[17:20]
            line_output['chainID'] = temp[21]
            line_output['resSeq'] = int(temp[22:26])
            line_output['iCode'] = temp[26]
            line_output['x'] = float(temp[30:38])
            line_output['y'] = float(temp[38:46])
            line_output['z'] = float(temp[46:54])
            line_output['occupancy'] = float(temp[54:60])
            line_output['tempFactor'] = float(temp[60:66])
            line_output['element'] = temp[76:78]
            line_output['charge'] = temp[78:80]
            yield line_output

def pdb_format(line):
    """
    Takes in a dictionary with the fields as defined
    by the last function and outputs a line.
    """
    out_string = ""
    out_string += line['ATOM']
    out_string += str(line['serial']).rjust(5)
    out_string += " "
    out_string += line['name']
    out_string += line['altLoc']
    out_string += line['resName']
    out_string += " "
    out_string += line['chainID']
    out_string += str(line['resSeq']).rjust(4)
    out_string += line['iCode']
    out_string += " "*3
    out_string += ("{0:.3f}".format(line['x'])).rjust(8)
    out_string += ("{0:.3f}".format(line['y'])).rjust(8)
    out_string += ("{0:.3f}".format(line['z'])).rjust(8)
    out_string += ("{0:.2f}".format(line['occupancy'])).rjust(6)
    out_string += ("{0:.2f}".format(line['tempFactor'])).rjust(6)
    out_string += " "*10
    out_string += line['element']
    out_string += line['charge']
    return out_string

def filter_line(line_generator, predicate):
    """
    Takes in an iterator that generates lines and
    a predicate function and filters through.
    """
    for line in line_generator:
        if predicate(line):
            yield line

if __name__ == "__main__":
    from sys import argv
    for l in parse_file(argv[-1]):
        print(pdb_format(l))
