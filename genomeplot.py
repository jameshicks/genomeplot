#!/usr/bin/python

import sys
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from math import ceil

try:
    import seaborn as sns
    seaborn = True
except ImportError:
    seaborn = False

if seaborn:
    sns.set(context='talk', style='nogrid')

# Parse arguments
parser = argparse.ArgumentParser(description='Create Manhattan style plots of results from'
                                 'genome-wide screens')
parser.add_argument('-f', '--file', help="Results file from screen", 
                    metavar='FILE', dest='file')
parser.add_argument('--ymax', type=float, help='Maximum value for the y axis',
                    dest='ymax', default=None)
parser.add_argument('-l','--lines', help="Make a line graph", dest='lines', action='store_true')
parser.add_argument('-p','--points', help="Make a manhattan plot", action='store_true', dest='points')
parser.add_argument('-t','--title', help="Set plot title")
args = parser.parse_args()


sigalpha = 3.0
sugalpha = -np.log10(0.05)

# Read data
gwas = pd.read_csv(args.file)

# Calculate positions on X axis
chroms = sorted(set(gwas.chr))
for c in chroms:
    cstart = min(gwas.ix[gwas.chr==c,'pos'])
    gwas.ix[gwas['chr'] == c,'cumpos'] = gwas.ix[gwas['chr'] == c,:]['pos'] - cstart
    lastcumpos = max(gwas.ix[gwas.chr == (c-1),'cumpos']) if c != 1 else 0        
    gwas.ix[gwas.chr==c,'cumpos'] = gwas.ix[gwas.chr==c,'pos'] + lastcumpos

# Set plot limits
xmin = gwas.cumpos.min()
xmax = gwas.cumpos.max()
maxstat = ceil(max(-np.log10(gwas['p']))) if not args.ymax else args.ymax

plt.xlim(xmin,xmax)
plt.ylim(0,maxstat)

# Draw lines separating chromosomes
chrombreaks = np.array([max(gwas.ix[gwas.chr == x,'cumpos']) for x in chroms])
plt.vlines(chrombreaks, 0, maxstat, color='gray', alpha=0.1)

# Plot the data!
for c in chroms:
    ss = gwas.ix[gwas.chr == c,:]
    if args.lines:
        plt.plot(ss.cumpos, -np.log10(ss.p), linewidth=0.75)
    if args.points:
        plt.plot(ss.cumpos, -np.log10(ss.p), '.')

# Set x-axis ticks
xticks=[gwas.ix[gwas.chr==x,'cumpos'].median() for x in chroms]
plt.xticks(xticks, chroms)
plt.tick_params(axis='x', which='major', labelsize=9)

# Draw significance lines
if sigalpha:
    plt.hlines(sigalpha, xmin, xmax, color='red')
if sugalpha:
    plt.hlines(sugalpha, xmin, xmax, color='blue')

# Write title
if args.title:
    plt.title(args.title)

# Label axes
plt.xlabel('Chromosome')
plt.ylabel('-log10(p)')

# Show plot
plt.show()
