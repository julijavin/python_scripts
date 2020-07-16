#!/usr/bin/env python3
"""
Created on Thu Jan 19 15:39:59 2017

@author: Julija

Plot band energies from VASP
"""

import argparse
import os, sys
import matplotlib.pyplot as pl
import matplotlib
import numpy as np

parser = argparse.ArgumentParser(description="plot DOS and band energies")
parser.add_argument('dir_path', help="directory that contains EIGENVAL and DOSCAR files")
parser.add_argument('-s', '--spin_pol', default=True, help="specifies that calculation is spin-polarized, default=True")
parser.add_argument('-f', '--shift_fermi', default=True, help="specifies to shift fermi level to zero, default=True")

args = parser.parse_args()

dir_path = args.dir_path
spin_pol = args.spin_pol
shift_fermi = args.shift_fermi

EIGENVAL_path = os.path.join(dir_path, 'EIGENVAL')
DOSCAR_path = os.path.join(dir_path, 'DOSCAR')

font_s = 12
font = {'family' : 'Arial',
        'size'   : font_s}
matplotlib.rc('font', **font)

with open(DOSCAR_path) as d:
    line1 = d.readline()
    line2 = d.readline()
    line3 = d.readline()
    line4 = d.readline()
    line5 = d.readline()
    info = d.readline().split()
    efermi = float(info[3])
    pts = int(float(info[2]))
    total_DOS_en = np.zeros((pts, 1))
    total_DOS_up = np.zeros((pts, 1))
    if spin_pol:
        total_DOS_dn = np.zeros((pts, 1))
    for ii in range(pts):
        sline = d.readline().split()
        total_DOS_en[ii] = float(sline[0])
        total_DOS_up[ii] = float(sline[1])
        if spin_pol:
            total_DOS_dn[ii] = float(sline[2])
    d.close()

if spin_pol == True:
    with open(EIGENVAL_path) as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        line4 = f.readline()
        comment = f.readline()
        unknown, npoints, nbands = [int(x) for x in f.readline().split()]
    
        blankline = f.readline()
    
        band_energies_1 = [[] for i in range(nbands)]
        band_energies_2 = [[] for i in range(nbands)]
    
        for i in range(npoints):
            x, y, z, weight = [float(x) for x in f.readline().split()]
    
            for j in range(nbands):
                fields = f.readline().split()
                id1, energy_1, energy_2 = int(float(fields[0])), float(fields[1]), float(fields[2])
                band_energies_1[id1-1].append(energy_1)
                band_energies_2[id1-1].append(energy_2)
            blankline = f.readline()
        f.close()
    
    rescale_up = npoints*0.5/np.amax(total_DOS_up)
    rescale_dn = npoints*0.5/np.amax(total_DOS_dn)
    rescale = max(rescale_up, rescale_dn)
    band_energies_1 = np.array(band_energies_1)
    band_energies_2 = np.array(band_energies_2)
    if shift_fermi:
        efermi_shift = efermi
    else:
        efermi_shift = 0
        
    band_energies_1 = band_energies_1 - efermi_shift
    band_energies_2 = band_energies_2 - efermi_shift
    efermi = efermi - efermi_shift
    total_DOS_en = total_DOS_en - efermi_shift
    
    fig, ax = pl.subplots(nrows=1, ncols=2, sharex=False, sharey=True, figsize=(10, 7))
    
    for i in range(nbands):
        ax[1].plot(range(npoints), np.array(band_energies_1[i]))
        ax[1].text(npoints-1, np.mean(np.array(band_energies_1[i])), i+1, fontdict=font) # up
        ax[1].plot(range(npoints), np.ones(len(range(npoints)))*efermi, '-', color = 'black', linewidth = 1.5)
        ax[1].text(1, efermi, 'E Fermi', fontdict=font)
        ax[1].plot(total_DOS_up*rescale, total_DOS_en, '-', color = 'black', linewidth = 1.5)
        ax[1].set_xticks([])  # no tick marks
        ax[1].set_xlabel('k-vector')
        #ax[1].set_title('spin up')
        
        ax[0].plot(-1*np.array(range(npoints)), np.array(band_energies_2[i]))
        ax[0].text(-npoints+1, np.mean(np.array(band_energies_2[i])), i+1, fontdict=font) # down
        ax[0].plot(-1*np.array(range(npoints)), np.ones(len(range(npoints)))*efermi, '-', color = 'black', linewidth = 1.5)
        #ax[0].text(-npoints+1, efermi, 'E Fermi', fontdict=font)
        ax[0].plot(-total_DOS_dn*rescale, total_DOS_en, '-', color = 'black', linewidth = 1.5)
        ax[0].set_xticks([])  # no tick marks
        ax[0].set_xlabel('k-vector')
        ax[0].set_ylabel('energy (eV)')
        #ax[0].set_title('spin down')
        
    pl.suptitle(dir_path)
    pl.tight_layout()
    pl.subplots_adjust(top=0.950)
        
elif spin_pol == False:
    fig = pl.figure()
    with open(EIGENVAL_path) as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        line4 = f.readline()
        comment = f.readline()
        unknown, npoints, nbands = [int(x) for x in f.readline().split()]
    
        blankline = f.readline()
    
        band_energies_1 = [[] for i in range(nbands)]
    
        for i in range(npoints):
            x, y, z, weight = [float(x) for x in f.readline().split()]
    
            for j in range(nbands):
                fields = f.readline().split()
                id1, energy_1 = int(float(fields[0])), float(fields[1])
                band_energies_1[id1-1].append(energy_1)
            blankline = f.readline()
        f.close()
    
    band_energies_1 = np.array(band_energies_1)
    if shift_fermi:
        efermi_shift = efermi
    else:
        efermi_shift = 0
    
    rescale = npoints/np.amax(total_DOS_up)
    
    if shift_fermi:
        band_energies_1 = band_energies_1 - efermi_shift
        efermi = efermi - efermi_shift
        total_DOS_en = total_DOS_en - efermi_shift
    ax1 = pl.subplot(111)
    for i in range(nbands):
        pl.plot(range(npoints), np.array(band_energies_1[i]))
        pl.text(npoints-1, np.mean(np.array(band_energies_1[i])), i+1, fontdict=font) # up
        pl.plot(total_DOS_up*rescale, total_DOS_en, '-', color = 'black', linewidth = 1.5)
    
    pl.plot(range(npoints), np.ones(len(range(npoints)))*efermi, '-', color = 'black', linewidth = 2)
    pl.text(npoints, efermi, 'E Fermi', fontdict=font)
    pl.xlabel('k-vector')
    pl.ylabel('Energy (eV)')
    pl.title(dir_path)

pl.show()
