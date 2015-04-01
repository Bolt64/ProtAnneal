#!/usr/bin/env python

import os

Zdope_command = "/usr/local/bin/Zdope"

def get_dope_score(pdb_file):
    result = os.popen("{0} {1}".format(Zdope_command, pdb_file))
    result_lines = [l.strip() for l in result]
    score_line = result_lines[-2]
    return float(score_line.strip().split(":")[1])
