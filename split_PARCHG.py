#!/usr/bin/env python3
"""
Created on Fri Feb 24 11:12:50 2017

@author: Julija
"""

import argparse
import os, sys
import numpy as np
import math

parser = argparse.ArgumentParser(description="split PARCHG by spin")
parser.add_argument('PARCHG_path', help="path to PARCHG file")
parser.add_argument('-u', '--up_out_path', default="PARCHG_up.vasp", help="output path for up spin channel")
parser.add_argument('-d', '--down_out_path', default="PARCHG_dn.vasp", help="output path for down spin channel")

args = parser.parse_args()

file_path_PARCHG = args.PARCHG_path
file_path_up = args.up_out_path
file_path_dn = args.down_out_path

with open(file_path_PARCHG) as f:
    lines = f.readlines()
    atoms = [int(x) for x in lines[6].split()]
    n_header = sum(atoms) + 9
    block_lists = []
    add_block = True
    for line in lines[n_header:]:
        sline = line.split()
        if len(sline)==3:
            tmp_block = []
        elif len(sline)<3:
            if add_block:
                block_lists.append(tmp_block)
                add_block = False
        else:
            tmp_line = [float(x) for x in sline]
            tmp_block.append(tmp_line)
    block_lists.append(tmp_block)
    block1_float = block_lists[0]
    block2_float = block_lists[1]

up = []
dn = []

for en1, val1 in enumerate(block1_float):
    tmp_up = []
    tmp_dn = []
    for en2, val2 in enumerate(val1):
        tmp = (block1_float[en1][en2] + block2_float[en1][en2])/2
        tmp_up.append(tmp)
        tmp = (block1_float[en1][en2] - block2_float[en1][en2])/2
        tmp_dn.append(tmp)
    up.append(tmp_up)
    dn.append(tmp_dn)

with open(file_path_up, "w") as f1:
    f1.writelines(lines[0:n_header+1])
    for up_row in up:
        for up_item in up_row:
            f1.write("%s\t" % up_item)
        f1.write("\n")
    f1.close()

with open(file_path_dn, "w") as f2:
    f2.writelines(lines[0:n_header+1])
    for dn_row in dn:
        for dn_item in dn_row:
            f2.write("%s\t" % dn_item)
        f2.write("\n")
    f2.close()

