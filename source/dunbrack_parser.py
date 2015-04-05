#!/usr/bin/env python

def parse_dunbrack(dunbrack_file):
    dunbrack_dictionary = {}
    for l in open(dunbrack_file):
        line = l.strip().split(",")
        main_key = line[0][0:3], int(line[0][3:7]), int(line[0][7:11])
        if not main_key in dunbrack_dictionary:
            dunbrack_dictionary[main_key] = {}
        prob = float(line[1])
        chis_raw = [float(i) for i in line[2:6]]
        chis =[i for i in chis_raw if i!=0]
        dunbrack_dictionary[main_key][prob]=chis
    return dunbrack_dictionary
