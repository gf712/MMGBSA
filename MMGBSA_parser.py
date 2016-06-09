#!/usr/bin/env python
# Author: Gil Ferreira Hoben
# Script to process data collected from MMGBSA.py from Amber

import subprocess
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys


# Command line parser
parser = argparse.ArgumentParser(usage="""{} """.
                                 format(sys.argv[0]),
                                 epilog="""Script to extract information from 
                                 MMGBSA.py Amber 15 output.""")

parser.add_argument("-i", "--input", help="""Input directory containing
											output of MMGBSA.py.""", 
											type=str)
parser.add_argument("-o", "--output", help="""Output directory for all the 
											generated files.""", 
											type=str)
parser.add_argument("-fo", "--output_file", help="""Output file name.""", 
											type=str)
parser.add_argument("-pt", "--plot_title", help="""Plot title.""", 
											default='$\Delta$Total Energy - 50 frame rolling average',
											type=str)
parser.add_argument("-v", "--verbose", help="""Switch verbose on/off.
                    							Default is True.""", 
                    							default=True, type=bool)
parser.add_argument("-p", "--plot", help="""Switch plotting on/off.
                    							Default is True.""", 
                    							default=True, type=bool)
args = parser.parse_args()

def main(data_dir, output_dir, output_file, verbose, plot, plot_title):

	# first move to data directory and create an Analysis folder to dump all the data
	# mmgbsa_parser.sh is called to extract all the data
	os.chdir(glob(data_dir))
	subprocess.call(['bash', 'mmpbsa_parser.sh'])
	os.chdir(output_dir)

	# Get the data from each file created by mmgbsa_parser.sh to calculate delta total
	complex_gb_coul = np.loadtxt('_MMPBSA_complex_gb.coul')
	complex_gb_polar = np.loadtxt('_MMPBSA_complex_gb.polar')
	complex_gb_vdw = np.loadtxt('_MMPBSA_complex_gb.vdw')
	complex_gb_surf = np.loadtxt('_MMPBSA_complex_gb.surf')
	ligand_gb_coul = np.loadtxt('_MMPBSA_ligand_gb.coul')
	ligand_gb_polar = np.loadtxt('_MMPBSA_ligand_gb.polar')
	ligand_gb_vdw = np.loadtxt('_MMPBSA_ligand_gb.vdw')
	ligand_gb_surf = np.loadtxt('_MMPBSA_ligand_gb.surf')
	receptor_gb_coul = np.loadtxt('_MMPBSA_receptor_gb.coul')
	receptor_gb_polar = np.loadtxt('_MMPBSA_receptor_gb.polar')
	receptor_gb_vdw = np.loadtxt('_MMPBSA_receptor_gb.vdw')
	receptor_gb_surf = np.loadtxt('_MMPBSA_receptor_gb.surf')

	complex_total = complex_gb_coul + complex_gb_polar + complex_gb_vdw + complex_gb_surf
	receptor_total = receptor_gb_coul + receptor_gb_polar + receptor_gb_vdw + receptor_gb_surf
	ligand_total = ligand_gb_coul + ligand_gb_polar + ligand_gb_vdw + ligand_gb_surf

	if verbose:

		print 'Average complex Energy: %.2f kcal/mol' % complex_togtal.mean()
		print 'Average receptor Energy: %.2f kcal/mol' % receptor_total.mean()
		print 'Average ligand Energy: %.2f kcal/mol' % liand_total.mean()

	# calculate delta total
	delta_total = complex_total - receptor_total - ligand_total

	if verbose:

		print 'Average Delta Total: %.2f kcal/mol' % delta_total.mean()
		print 'Standard Deviation of Delta Total' % delta_total.std()


	if plot:
		plt.figure(figsize=(12,12))
		plt.title(plot_title, size = 22)
		plt.plot(pd.DataFrame(delta_total).rolling(window=50,center=False).mean()[0], 
			label = 'Average $\Delta$Total = %.2f kcal/mol\nStandard deviation = %.2f' % (delta_total.mean(), delta_total.std()))
		plt.ylabel('$\Delta$Total (kcal/mol)', size =15)
		plt.xlabel('Frame', size =15)
		plt.legend()

	np.savetxt(output_file+'_delta_total', delta_total)
	np.savetxt(output_file+'_complex_total', complex_total)
	np.savetxt(output_file+'_receptor_total', receptor_total)
	np.savetxt(output_file+'_ligand_total', ligand_total)

if __name__ == "__main__":
    main(args.input, args.output, args.output_file, args.verbose, args.plot, args.plot_title)  