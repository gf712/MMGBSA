#!/usr/bin/env python
# Authors: Gil Ferreira Hoben, Juan Eiros
# Script to process data collected from MMGBSA.py from Amber

from __future__ import print_function
import subprocess
import os
from glob import glob
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import argparse
import sys
import seaborn as sns
from msmexplorer import plot_trace
from matplotlib.ticker import FuncFormatter


# Import MMPBSA API from local AMBERHOME installation
# append AMBERHOME/bin to sys.path
sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))
from MMPBSA_mods import API as MMPBSA_API

# Plotting settings
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

    parser.add_argument("-i", "--info_file", type=str, default="_MMPBSA_info",
                        help="""The information file that is printed after the
                        MMPBSA calculation. Default is _MMPBSA_info""")
    parser.add_argument("-o", "--output_dir", type=str, default='plots',
                        help='Output directory for all the generated files. Default is plots/')
    parser.add_argument("-fo", "--output_file", type=str, default='plot',
                        help='Output file name. Default is plot.')
    parser.add_argument("-pt", "--plot_title", type=str, default='$\Delta$Total Energy',
                        help='Plot title.')
    parser.add_argument("-ts", "--time_step", type=float, default=0.02,
                        help='Time step (in ns) between frames')
    parser.add_argument("-p", "--plot", type=bool, default=True,
                        help='Switch plotting on/off. Default is True.')
    return parser.parse_args()


def getDataFrames(info_file, calculation_type='gb'):
    '''
    Reads an MMPBSA info file and returns a dictionary with four keys
    'complex', 'receptor', 'ligand' and 'difference'
    each with an associated pd.DataFrame with all the corresponding data
    '''
    data = MMPBSA_API.load_mmpbsa_info(info_file)  # data is a dictionary
    data_dfs = OrderedDict.fromkeys(['complex', 'receptor', 'ligand', 'difference'])
    for key in data[calculation_type]:
        df = pd.DataFrame.from_dict(data[calculation_type][key])
        data_dfs[key] = df
    # Calculate the difference between complex, receptor and ligand contributions
    data_dfs['difference'] = data_dfs['complex'] - data_dfs['receptor'] - data_dfs['ligand']
    return data_dfs


def getFileName(info_file):
    '''
    Cleans the args.info_file by getting only the str representation of the file name
    Couldn't be bothered to write a regexp so:
        Starts from the end and stops whenever it finds a /
        new_str is built in reversed so that's why it's returned as [::-1] (to revert it again)
    '''
    new_str = ''
    for letter in reversed(info_file):
        if letter != '/':
            new_str += letter
        else:
            break
    return new_str[::-1]


def saveDataFrames(dict_of_dfs):
    """
    Saves the pd.DataFrames of each key in the dict_of_dfs as pickle objects
    """
    for key, val in dict_of_dfs.items():
        val.to_pickle(key + '.pkl')


def makeTimeSeriesPlots(dict_of_dfs, n_frames):
    names = ['Complex Contribution', 'Receptor Contribution',
             'Ligand Contribution', '$\Delta$ Total']
    plt.figure(figsize=(12, 12))
    plt.suptitle(args.plot_title, size=22)
    for i, df in enumerate(dict_of_dfs.values()):
        # Create DataFrame, with Energy and Time columns
        df_ene = pd.DataFrame({
            'Energy': df['TOTAL'],
            # Moving avg
            'Energy_avg': df['TOTAL'].rolling(window=max(1, int(n_frames / 100))).mean(),
            'Time': pd.Series([x * args.time_step for x in range(0, n_frames)])
        })
        plt.subplot(2, 2, i + 1)
        plt.title(names[i], size=16)
        # actual values
        plt.plot(df_ene['Time'], df_ene['Energy'],
                 alpha=0.2, color='#1f77b4')
        # moving avg
        plt.plot(df_ene['Time'], df_ene['Energy_avg'],
                 color='#1f77b4')
        # Shared limits for Y axis for the complex and receptor plots
        if i < 2:
            y_min = min(dict_of_dfs['complex']['TOTAL'].min(), dict_of_dfs['receptor']['TOTAL'].min())
            y_max = max(dict_of_dfs['complex']['TOTAL'].max(), dict_of_dfs['receptor']['TOTAL'].max())
            # Make whitespace above and bottom to be 2.5% of plot
            # so diff between max and min values is 95%
            total_y_space = abs(y_min - y_max) / 0.95
            offsetY = total_y_space * 0.025
            plt.ylim(round(y_min - offsetY), round(y_max + offsetY))
        plt.xlim(0, round(n_frames * args.time_step))
        plt.ylabel('$\Delta G (kcal \cdot mol^{-1}$)', size=15)
        plt.xlabel('Time (ns)', size=15)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.2, top=.9)
    plt.savefig((args.output_file + '.pdf'))
    plt.savefig((args.output_file + '.png'), dpi=300)


def to_ns(x, pos):
    return '{:.0f}'.format(x * args.time_step)


def plot_individual_energy(data, plot_trace_kwargs=None, save=True):
    """
    Individual plot of the Delta G time evolution
    :param data: np.array
    :param plot_trace_kwargs: kwargs for msmexplorer's plot_trace function
    :param save: Bool, save image or not
    :return ax, side_ax: The main and side axis produced by plot_trace
    """
    if plot_trace_kwargs is None:
        plot_trace_kwargs = {}
    # actual values
    ax, side_ax = plot_trace(data=data, window=1, alpha=0.2,
                             **plot_trace_kwargs)
    # mean and std lines in side_ax
    side_ax.axhline(y=data.mean())
    side_ax.axhline(y=data.mean() - data.std(), ls='--', lw=0.8)
    side_ax.axhline(y=data.mean() + data.std(), ls='--', lw=0.8)
    # annotate side_ax with actual values
    side_ax.annotate(f"\u03bc = {data.mean():.02f}", xy=(0.55, 0.90), xycoords='axes fraction')
    side_ax.annotate(f"\u03c3 = {data.std():.02f}", xy=(0.55, 0.85), xycoords='axes fraction')
    # moving avg in ax
    plot_trace(data=data, window=max(1, int(len(data) / 100)),
               ax=ax, alpha=1,
               **plot_trace_kwargs)
    ax.set(title=args.plot_title, xlabel='Time (ns)',
           ylabel='$\Delta G (kcal \cdot mol^{-1}$)')
    formatter = FuncFormatter(to_ns)
    ax.xaxis.set_major_formatter(formatter)
    if save:
        f = plt.gcf()
        f.savefig((args.output_file + '_individual.pdf'))
        f.savefig((args.output_file + '_individual.png'), dpi=300)
    return ax, side_ax


def main(args):

    # First move to data directory and create an the output folder and cd into it

    if os.path.abspath(os.getcwd()) != os.path.dirname(os.path.realpath(args.info_file)):
        os.chdir(os.path.dirname(os.path.realpath(args.info_file)))
    data_dfs = getDataFrames(getFileName(args.info_file))
    n_frames = data_dfs['complex'].shape[0]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    os.chdir(args.output_dir)

    print('Output Directory: ', os.getcwd())
    # complex
    print('\ncomplex energy:\t\t%.2f %s %.2f\n'
          % (data_dfs['complex']['TOTAL'].mean(), chr(177), data_dfs['complex']['TOTAL'].std()))
    # receptor
    print('receptor energy:\t%.2f %s %.2f\n'
          % (data_dfs['receptor']['TOTAL'].mean(), chr(177), data_dfs['receptor']['TOTAL'].std()))
    # ligand
    print('ligand energy:\t\t%.2f %s %.2f\n'
          % (data_dfs['ligand']['TOTAL'].mean(), chr(177), data_dfs['ligand']['TOTAL'].std()))
    # difference
    print('difference energy:\t%.2f %s %.2f\n'
          % (data_dfs['difference']['TOTAL'].mean(), chr(177), data_dfs['difference']['TOTAL'].std()))

    # Make a 2x2 plot with time series of Delta G binding
    if args.plot:
        makeTimeSeriesPlots(data_dfs, n_frames)

    saveDataFrames(data_dfs)

    # Make an individual plot of the Difference
    plot_individual_energy(data=np.array(data_dfs['difference']['TOTAL']),
                           plot_trace_kwargs={'color': '#1f77b4'})

if __name__ == "__main__":
    args = parse_args()
    main(args)
