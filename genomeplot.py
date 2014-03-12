#!/usr/bin/python

import sys
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

sigalpha = 3.0
sugalpha = -np.log10(0.05)

gwas = pd.read_csv(sys.argv[1])

chroms = sorted(set(gwas.chr))
for c in chroms:
    if c == 1: 
        gwas.ix[gwas['chr'] == c,'cumpos'] = gwas.ix[gwas['chr'] == c,:]['pos']
    else:
        lastcumpos = max(gwas.ix[gwas.chr == (c-1),'cumpos'])
        gwas.ix[gwas.chr==c,'cumpos'] = gwas.ix[gwas.chr==c,'pos'] + lastcumpos

xmin = gwas.cumpos.min()
xmax = gwas.cumpos.max()
maxstat = ceil(max(-np.log10(gwas['p'])))

plt.xlim(xmin,xmax)
plt.ylim(0,maxstat)

chrombreaks = np.array([max(gwas.ix[gwas.chr == x,'cumpos']) for x in chroms])

plt.vlines(chrombreaks, 0, maxstat, color='gray', alpha=0.25)

for c in chroms:
    ss = gwas.ix[gwas.chr == c,:]
    plt.plot(ss.cumpos, -np.log10(ss.p))

xticks=[gwas.ix[gwas.chr==x,'cumpos'].min() for x in chroms]
plt.xticks(xticks, chroms)
if sigalpha:
    plt.hlines(sigalpha, xmin, xmax, color='red')
if sugalpha:
    plt.hlines(sugalpha, xmin, xmax, color='blue')

plt.xlabel('Chromosome')
plt.ylabel('-log10(p)')
plt.show()
