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
    print 'Seaborn library not installed. Plots would be prettier with it'
    seaborn = False

if seaborn:
    sns.set(context='talk', style='nogrid')

# Parse arguments
parser = argparse.ArgumentParser(description='Create Manhattan style plots of results from'
                                 'genome-wide screens')
parser.add_argument('-f', '--file', help="Results file from screen", 
                    metavar='FILE', dest='file', required=True)
parser.add_argument('-s','--stat', help='Statistic to plot on Y axis',default='p')
parser.add_argument('--ymin', type=float, help='Minimum value for the y axis', default=0)
parser.add_argument('--ymax', type=float, help='Maximum value for the y axis',
                    dest='ymax', default=None)
parser.add_argument('-l','--lines', help="Make a line graph", dest='lines', action='store_true')
parser.add_argument('-p','--points', help="Make a manhattan plot", action='store_true', dest='points')
parser.add_argument('-t','--title', help="Set plot title")
parser.add_argument('--ylab', help='Y axis label', default='p')
parser.add_argument('--log10', help='Transform the Y axis values as -log10(y)', action='store_true')
parser.add_argument('--significant', help='Draw significance threshold at value', type=float)
parser.add_argument('--suggestive', help="Draw a 'suggestive' threshold at y value", type=float)
parser.add_argument('--bonferroni', help="Draw significance/suggestive lines at"
                    " Bonferroni corrected thresholds", action='store_true')
args = parser.parse_args()


if args.log10:
    transform = lambda y: -np.log10(y)
else:
    transform = lambda y: y

# Read data
print 'Reading data...',
gwas = pd.read_csv(args.file)
print 'Done'


# Calculate positions on X axis
print 'Calculating layout'

chroms = sorted(set(gwas.chr))
for c in chroms:
    cstart = gwas.ix[gwas.chr==c,'pos'].min()
    lastcumpos = max(gwas.ix[gwas.chr == (c-1),'cumpos']) + 1 if c != 1 else 0
    gwas.ix[gwas['chr'] == c,'cumpos'] = gwas.ix[gwas['chr'] == c,'pos'] - cstart + lastcumpos

# Set plot limits
xmin = gwas.cumpos.min()
xmax = gwas.cumpos.max()
maxstat = ceil(transform(gwas[args.stat]).max()) if not args.ymax else args.ymax

plt.xlim(xmin, xmax)
plt.ylim(args.ymin, maxstat)

# Draw lines separating chromosomes
chrombreaks = np.array([max(gwas.ix[gwas.chr == x,'cumpos']) for x in chroms])
plt.vlines(chrombreaks, args.ymin, maxstat, color='gray', alpha=0.1)

print 'Plotting data'
# Plot the data!
for c in chroms:
    ss = gwas.ix[gwas.chr == c,:]
    if args.lines:
        plt.plot(ss.cumpos, transform(ss[args.stat]), linewidth=0.75)
    if args.points:
        plt.plot(ss.cumpos, transform(ss[args.stat]), '.')

# Set x-axis ticks
xticks=[gwas.ix[gwas.chr==x,'cumpos'].median() for x in chroms]
plt.xticks(xticks, chroms)
plt.tick_params(axis='x', which='major', labelsize=8)

# Draw significance lines

if args.bonferroni:
    plt.hlines(transform(gwas.shape[0]), xmin, xmax, color='red')
    plt.hlines(transform(0.05), xmin, xmax, color='blue')
else:
    if args.significant:
        plt.hlines(transform(args.significant), xmin, xmax, color='red')
    if args.suggestive:
        plt.hlines(transform(args.suggestive), xmin, xmax, color='blue')

# Write title
if args.title:
    plt.title(args.title)

# Label axes
plt.xlabel('Chromosome')
ylab = '-log10({0})'.format(args.ylab) if args.log10 else args.ylab  
plt.ylabel(ylab)

# Show plot
plt.show()
