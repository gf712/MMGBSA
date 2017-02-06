#!/usr/bin/env python
# Authors: Gil Ferreira Hoben, Juan Eiros
# Script to process data collected from MMGBSA.py from Amber

from __future__ import print_function
import subprocess
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import argparse
import sys
import seaborn as sns
sns.set_style('whitegrid')
matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['font.family'] = "sans-serif"


# Command line parser
def parse_args():
    parser = argparse.ArgumentParser(usage="""{} """.
                                     format(sys.argv[0]),
                                     epilog="""Script to extract information from
                                     MMGBSA.py Amber 15 output.""")

    parser.add_argument("-i", "--input_dir", help="""Input directory containing
                                                 output of MMGBSA.py.""",
                        type=str)
    parser.add_argument("-o", "--output_dir", help="""Output directory for all the
                                                generated files.""",
                        default='plots',
                        type=str)
    parser.add_argument("-fo", "--output_file", help="""Output file name.""",
                        default='plot',
                        type=str)
    parser.add_argument("-pt", "--plot_title", help="""Plot title.""",
                        default='$\Delta$Total Energy',
                        type=str)
    parser.add_argument("-ts", "--time_step", help="Time step (in ns) between frames",
                        default=0.02, type=float)
    parser.add_argument("-v", "--verbose", help="""Switch verbose on/off.
                                                    Default is True.""",
                        default=True, type=bool)
    parser.add_argument("-p", "--plot", help="""Switch plotting on/off.
                                                    Default is True.""",
                        default=True, type=bool)
    return parser.parse_args()


def main(args):

    # first move to data directory and create an Analysis folder to dump all the data
    # mmgbsa_parser.sh is called to extract all the data
    os.chdir(args.input_dir)
    subprocess.call(['bash', 'mmpbsa_parser.sh'])

    # Get the data from each file created by mmgbsa_parser.sh to calculate delta total
    # complex_gb_coul = np.loadtxt('_MMPBSA_complex_gb.coul')
    # complex_gb_polar = np.loadtxt('_MMPBSA_complex_gb.polar')
    # complex_gb_vdw = np.loadtxt('_MMPBSA_complex_gb.vdw')
    # complex_gb_surf = np.loadtxt('_MMPBSA_complex_gb.surf')
    # ligand_gb_coul = np.loadtxt('_MMPBSA_ligand_gb.coul')
    # ligand_gb_polar = np.loadtxt('_MMPBSA_ligand_gb.polar')
    # ligand_gb_vdw = np.loadtxt('_MMPBSA_ligand_gb.vdw')
    # ligand_gb_surf = np.loadtxt('_MMPBSA_ligand_gb.surf')
    # receptor_gb_coul = np.loadtxt('_MMPBSA_receptor_gb.coul')
    # receptor_gb_polar = np.loadtxt('_MMPBSA_receptor_gb.polar')
    # receptor_gb_vdw = np.loadtxt('_MMPBSA_receptor_gb.vdw')
    # receptor_gb_surf = np.loadtxt('_MMPBSA_receptor_gb.surf')

    # complex_total = complex_gb_coul + complex_gb_polar + complex_gb_vdw + complex_gb_surf
    # receptor_total = receptor_gb_coul + receptor_gb_polar + receptor_gb_vdw + receptor_gb_surf
    # ligand_total = ligand_gb_coul + ligand_gb_polar + ligand_gb_vdw + ligand_gb_surf
    complex_total = np.loadtxt('./Analysis/data._MMPBSA_complex_gb')
    receptor_total = np.loadtxt('./Analysis/data._MMPBSA_receptor_gb')
    ligand_total = np.loadtxt('./Analysis/data._MMPBSA_ligand_gb')
    delta_total = complex_total - receptor_total - ligand_total
    n_frames = delta_total.shape[0]
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    os.chdir(args.output_dir)

    print('Output Directory: ', os.getcwd())

    if args.verbose:

        print("""\nAverage complex Energy: %.2f kcal/mol\n
Average receptor Energy: %.2f kcal/mol\n
Average ligand Energy: %.2f kcal/mol\n
Average Delta Total: %.2f kcal/mol\n
Standard Deviation of Delta Total: %.2f\n
Loaded frames: %d\n""" %
              (complex_total.mean(), receptor_total.mean(), ligand_total.mean(),
               delta_total.mean(), delta_total.std(), n_frames))

    names = ['Complex Contribution', 'Receptor Contribution', 'Ligand Contribution', '$\Delta$ Total']

    if args.plot:
        plt.figure(figsize=(12, 12))
        plt.suptitle(args.plot_title, size=22)
        for i, data in enumerate([complex_total, receptor_total, ligand_total, delta_total]):
            # Create DataFrame, with Energy and Time columns
            df = pd.DataFrame({
                'Energy': pd.Series(data).rolling(window=max(1, int(n_frames / 100))).mean(),  # Moving avg
                'Time': pd.Series([x * args.time_step for x in range(0, n_frames)])
            })
            plt.subplot(2, 2, i + 1)
            plt.title(names[i], size=16)
            plt.plot(df['Time'], df['Energy'],
                     label='avg = %.2f kcal/mol\nstd = %.2f kcal/mol' % (data.mean(), data.std()))
            # Shared limits for Y axis for the complex and receptor plots
            if i < 2:
                y_min = min(complex_total.min(), receptor_total.min())
                y_max = max(complex_total.max(), receptor_total.max())
                # Make whitespace above and bottom to be 2.5% of plot
                # so diff between max and min values is 95%
                total_y_space = abs(y_min - y_max) / 0.95
                offsetY = total_y_space * 0.025
                plt.ylim(round(y_min - offsetY), round(y_max + offsetY))
            plt.ylabel('$\Delta$G (kcal/mol)', size=15)
            plt.xlabel('Time (ns)', size=15)
            plt.legend(prop={'size': 8})
            plt.tight_layout()
            plt.subplots_adjust(hspace=0.2, top=.9)
        plt.savefig((args.output_file + '.pdf'))

    np.savetxt(args.output_file + '_delta_total', delta_total)
    np.savetxt(args.output_file + '_complex_total', complex_total)
    np.savetxt(args.output_file + '_receptor_total', receptor_total)
    np.savetxt(args.output_file + '_ligand_total', ligand_total)


if __name__ == "__main__":
    args = parse_args()
    main(args)
