#!/usr/bin/env python

def parse_file(pdb_file):
    for line in open(pdb_file):
        temp=line.strip().ljust(80," ")
        line_output={}
        line_output['ATOM']=temp[0:6]
        line_output['serial']=int(temp[6:11])
        line_output['name']=temp[12:16]
        line_output['altLoc']=temp[16]
        line_output['resName']=temp[17:20]
        line_output['chainID']=temp[21]
        line_output['resSeq']=int(temp[22:26])
        line_output['iCode']=temp[26]
        line_output['x']=float(temp[30:38])
        line_output['y']=float(temp[38:46])
        line_output['z']=float(temp[46:54])
        line_output['occupancy']=float(temp[54:60])
        line_output['tempFactor']=float(temp[60:66])
        line_output['element']=temp[76:78]
        line_output['charge']=temp[78:80]
        yield line_output
