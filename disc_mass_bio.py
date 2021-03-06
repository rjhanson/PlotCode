#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

# Element indicies in the exoplanets.org database
PER = 0
MASS = 1
RMS = 2
K = 3
DATE = 4
MSTAR = 5
V = 6
ECC = 7
RSTAR = 8
TEFF = 9
A = 10
NOBS = 11

cols = ['PER','MSINI', 'RMS', 'K', 'DATE', 'MSTAR', 'V', 'ECC', 'RSTAR', 'TEFF', 'A','NOBS']

valueTable = np.empty(0)

names = []
with open('exoplanets.csv','r') as f:
    line = f.readline()
    l = line.split(',')
    index = [l.index(p) for p in cols]
    discovery_idx = l.index('PLANETDISCMETH')
    name_idx = l.index('NAME')
    for line in f:
        l = line.split(',')
        # Restrict results to RV detections only
        if l[discovery_idx] != 'RV': continue
        row = []
        for i in index:
            try:
                val = float(l[i])
            # If there was no value, default to zero
            except ValueError:
                val = 0
            
            row.append(val)
        names.append(l[name_idx].strip())
        valueTable = np.append(valueTable, np.array(row))

names = np.array(names)

valueTable = valueTable.reshape(len(valueTable)/len(cols), len(cols))

vt = valueTable


# Calculate the equilibrium temprature values
n1 = np.sqrt((vt[:,RSTAR]*0.00464913034) / (2*vt[:,A]))
n2 = vt[:,TEFF] / ((1-vt[:,ECC]**2)**(1.0/8.0))

Teq = n1*n2

# Calculate the AstroBiology Value
em = vt[:,MASS] * 317.828133
t1 = (1.0/vt[:,MSTAR])**(1.0/3.0)
t2 = np.exp(-1* (np.log10(em) / 0.2)**2.0 )
t3 = np.exp(-1* ((Teq - 273)/30.0)**2.0 )
t4 = (2.5**(12-vt[:,V]))**(1.0/2.0)


ABV = t1*t2*t3*t4

# Normalize the ABV values to between 0-75
jdx, = np.where(ABV == 0)
ABV[jdx] = 1e-300

ABV = np.log(ABV)

ABV = 75 * (ABV-min(ABV)) / max(ABV-min(ABV))


# Only plot points with an RMS greater than 0
sort, = np.where(vt[:,RMS] > 0.0)

# Randomize detections to fill the year they were published in
d = vt[:,DATE][sort] + np.random.rand(len(vt[:,DATE][sort]))
jdx, = np.where(d > 2014)
d[jdx] = 2014.3

jdx, = np.where(vt[:,K][sort]/vt[:,RMS][sort] < 0.9)

# Add Text for 2 points on the plot
for idx in jdx:
    x = d[idx]
    y = (vt[:,K][sort]/vt[:,RMS][sort])[idx]
    if '20794' in names[sort][idx]:
        plt.text(x, y-0.1, names[sort][idx], fontsize=14, ha='center', va='top')
    if 'Cen B b' in names[sort][idx]:
        plt.text(x, y-0.04, r'$\mathregular{\alpha}$ Cen B b', fontsize=14, ha='center', va='top')

# K/RMS value
y = vt[:,K][sort]/vt[:,RMS][sort]

# Plot the actual points
plt.scatter( d, y, c=ABV[sort], cmap=plt.get_cmap('YlGn'),edgecolor='gray', s=150, vmin=0, vmax=100, zorder=20)


# Create the color bar
cbar = plt.colorbar()
cbar.set_label('Astrobiological Value',fontsize=22)


plt.plot([1990,2020],[1,1],lw=4., alpha=0.2,c='black',zorder=10)


plt.xlabel('Discovery Date', fontsize=24)
plt.ylabel('Signal Strength', fontsize=24)

plt.yscale('log')

# Set the tick label sizes
plt.tick_params(labelsize=20)

# Maunally set the y labels to not be scientific notation
ax = plt.gca()

labels = [item.get_text() for item in ax.get_yticklabels()]
labels[3] = '100'
labels[2] = '10'
labels[1] = '1'

ax.set_yticklabels(labels)

plt.xlim([1994, 2015])
plt.ylim([10**(-0.5),10**(2.2)])

plt.gcf().set_size_inches(16,8)

# Save the figure
plt.savefig('Signal_vs_Date.pdf', bbox='tight')




            
        
