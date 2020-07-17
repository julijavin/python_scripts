#!/usr/bin/env python3
import argparse
import os, sys
import numpy as np

parser = argparse.ArgumentParser(description="find position of site located at the mean of specified atoms")
parser.add_argument('POSCAR', help="name of POSCAR-type file in current directory from which to retrieve atom positions")
parser.add_argument('atoms', nargs='*', type=int, help="coordinating atoms, specified by number in POSCAR")
parser.add_argument('-a', '--adjustments', help="file name with x y z adjustments for each atom per row, using spaces as delimiter")
parser.add_argument('-i', '--inside', default=True, help="moves coordinates into unit cell if true (default), else returns actual mean")
parser.add_argument('-v', '--vesta', help="file name that contains vesta terminal output for reading atoms and adjustments")

args = parser.parse_args()

if (args.adjustments and args.vesta) or (args.atoms and args.vesta):
    print("Reading atoms and adjustments from vesta output")

# read in x y z of user-specified atoms
with open(args.POSCAR, 'r') as f:
    lines = f.readlines()
    all_atoms = []
    for i,line in enumerate(lines[8:]):
        tmp = line.split()
        tmp = [float(x) for x in tmp[0:3]]
        all_atoms.append(tmp)

if args.atoms:
    atoms = args.atoms

# read in user-specified adjustments
if args.adjustments:
    with open(args.adjustments, 'r') as f:
        lines = f.readlines()
        adj = []
        for i,line in enumerate(lines):
            tmp = line.split()
            tmp = [float(x) for x in tmp]
            adj.append(tmp)
        adj = np.array(adj)

# read in atoms and adjustments from vesta print
if args.vesta:
    atoms = []
    adj = []
    with open(args.vesta, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            sline = line.split()
            if len(sline) > 0:
                if sline[0] == "Atom:":
                    atoms.append(int(sline[1]))
                    spline_1 = line.split('(')[1]
                    spline_2 = spline_1.split(')')[0]
                    spline = spline_2.split(',')
                    tmp = []
                    for j, val in enumerate(spline):
                        tmp.append(int(val))
                    adj.append(tmp)
    adj = np.array(adj)


# shifted atoms
atoms_shifted = []
for i,atom in enumerate(atoms):
    tmp = all_atoms[atom-1]
    if args.adjustments or args.vesta:
        atoms_shifted.append(tmp+adj[i])
    else:
        atoms_shifted.append(tmp)

atoms_shifted = np.array(atoms_shifted)

# find mean
site = np.mean(atoms_shifted, axis=0)

if args.inside:
    for i, s in enumerate(site):
        while s < 0:
            s = s + 1
        while s > 1:
            s = s - 1
        site[i] = s
    print("site inside u.c.:\n%17.16f  %17.16f  %17.16f" %(site[0], site[1], site[2]))
else:
    print("direct site:\n%17.16f  %17.16f  %17.16f" %(site[0], site[1], site[2]))
