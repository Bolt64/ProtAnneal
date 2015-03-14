#!/usr/bin/env python

from os import system

def get_filenames(cluster_file):
    for line in open(cluster_file):
        if line.strip().endswith("*"):
            l=line.strip()
            index=l.index(">")
            yield l[index+1:index+5].lower()+".pdb"

def move(filename, source):
    system("mv "+filename+" "+source)

if __name__=="__main__":
    from sys import argv
    for l in get_filenames(argv[-2]):
        move(l, argv[-1])
