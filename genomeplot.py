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
parser = argparse.ArgumentParser(description='Create Manhattan style plots of '
                                 'results from genome-wide screens',
                                 prog='genomeplot')

parser.add_argument('-f', '--file', help="Results file from screen",
                    metavar='FILE', dest='file', required=True)
parser.add_argument('--explore', help=argparse.SUPPRESS, action='store_true',
                    default=False)

group = parser.add_argument_group('Plot style')
group.add_argument('-l', '--lines', help="Make a line graph", dest='lines',
                   action='store_true')
group.add_argument('-p', '--points', help="Make a manhattan plot",
                   action='store_true', dest='points')

group = parser.add_argument_group('Column options')
group.add_argument('-d', '--sep', default=',', help='Field delimiter',
                   dest='sep')
group.add_argument('--chr', default='chr',
                   help='Column label containing chromosome')
group.add_argument('--pos', default='pos',
                   help='Column label containing position')
group.add_argument('-s', '--stat',
                   help='Statistic to plot on Y axis', default='p')

group = parser.add_argument_group('Axis options')
group.add_argument('-t', '--title', help="Set plot title")
group.add_argument('--ylab', help='Y axis label', default=None)
group.add_argument('--log10', help='Transform the Y axis values as -log10(y)',
                   action='store_true')
group.add_argument('--ymin', type=float,
                   help='Minimum value for the y axis', default=0)
group.add_argument('--ymax', type=float, help='Maximum value for the y axis',
                   dest='ymax', default=None)

group = parser.add_argument_group('Significance lines')
group.add_argument('--significant', type=float,
                   help='Draw significance threshold at value')
group.add_argument('--suggestive', type=float,
                   help="Draw a 'suggestive' threshold at y value")
group.add_argument('--bonferroni', help="Draw significance/suggestive lines at"
                   " Bonferroni corrected thresholds", action='store_true')
group.add_argument('--kruglyak', help='Draws significance lines at the pvalue'
                   'thresholds proposed for linkage analysis by'
                   'Lander & Kruglyak (1995). Significant: 4.9e-5'
                   'Suggestive: 1.7e-3',
                   action='store_true')

args = parser.parse_args()


if args.log10:
    transform = lambda y: -np.log10(y)
else:
    transform = lambda y: y

if args.file.endswith('.gz'):
    comp = 'gzip'
elif args.file.endswith('.bz2'):
    comp = 'bz2'
else:
    comp = None

# Read data
try:
    print 'Reading data...',
    gwas = pd.read_csv(args.file, sep=args.sep, compression=comp)
    print 'Done'
except IOError:
    print 'Could not read file %s' % args.file

print 'Minimum statistic: %s' % gwas[args.stat].min()
print 'Maximum statistic: %s' % gwas[args.stat].max()
print

if args.explore:
    try:
        from IPython import embed
        embed()
    except ImportError:
        print 'ERROR: IPython not found!'
        exit(1)
# Calculate positions on X axis
print 'Calculating layout'

chroms = sorted(set(gwas[args.chr]))
for c in chroms:
    cstart = gwas.ix[gwas[args.chr] == c, args.pos].min()
    if c != 1:
        lastcumpos = max(gwas.ix[gwas[args.chr] == (c - 1), 'cumpos']) + 1
    else:
        lastcumpos = 0
    gwas.ix[gwas['chr'] == c, 'cumpos'] = \
        gwas.ix[gwas['chr'] == c, args.pos] - cstart + lastcumpos

# Set plot limits
xmin = gwas.cumpos.min()
xmax = gwas.cumpos.max()

if not args.ymax:
    maxstat = ceil(transform(gwas[args.stat]).max())
else:
    maxstat = args.ymax

plt.xlim(xmin, xmax)
plt.ylim(args.ymin, maxstat)

# Draw lines separating chromosomes
chrombreaks = np.array([max(gwas.ix[gwas[args.chr] == x, 'cumpos'])
                        for x in chroms])
plt.vlines(chrombreaks, args.ymin, maxstat, color='gray', alpha=0.1)

print 'Plotting data'
# Plot the data!
for c in chroms:
    ss = gwas.ix[gwas[args.chr] == c, :]
    if args.lines:
        plt.plot(ss.cumpos, transform(ss[args.stat]), linewidth=0.75)
    if args.points:
        plt.plot(ss.cumpos, transform(ss[args.stat]), '.')

# Set x-axis ticks
xticks = [gwas.ix[gwas[args.chr] == x, 'cumpos'].median() for x in chroms]
plt.xticks(xticks, chroms)
plt.tick_params(axis='x', which='major', labelsize=8)

# Draw significance lines

if args.bonferroni:
    plt.hlines(transform(gwas.shape[0]), xmin, xmax, color='red')
    plt.hlines(transform(0.05), xmin, xmax, color='blue')
elif args.kruglyak:
    plt.hlines(transform(4.9e-5), xmin, xmax, color='red')
    plt.hlines(transform(1.7e-3), xmin, xmax, color='blue')
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

if args.ylab:
    ylab = '-log10({0})'.format(args.ylab) if args.log10 else args.ylab
else:
    ylab = '-log10({0})'.format(args.stat) if args.log10 else args.stat
plt.ylabel(ylab)

# Show plot
plt.show()
