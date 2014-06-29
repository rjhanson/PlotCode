#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt

# Grab appropriate columns from exoplanets.org database
params = ['RMS','K','DATE']
with open('exoplanets.csv','r') as f:
    line = f.readline()
    lines = line.split(',')
    idx = [lines.index(p) for p in params]
    store = [ [] for i in idx ]
    for line in f:
        lines = line.split(',')
        try:
            vals = [float(lines[v]) for v in idx]
        except ValueError:
            pass
        else:
            for i,j in enumerate(vals):
                store[i].append(j)

# Name arrays for easy of use
rms = np.array(store[0])
k = np.array(store[1])
date = np.array(store[2])

fig, ax = plt.subplots()
fig.set_size_inches(7,7)

# Set up a properly shaded grid lines
ax.grid(True, which='major', alpha=0.9)
ax.grid(True, which='minor', alpha=0.35)

# Plot as points
cax = ax.scatter(rms,k,c=date, cmap=plt.get_cmap('rainbow'),s=100,zorder=100.,alpha=0.8,edgecolor='none')

# Add a diagonal black line
ax.plot([0.4,300],[0.4,300],c='black',zorder=50.)

# Add extra planets not listed in exoplanets.org (Steve's planets)
# GJ 581 f,g; HD 40307 e,f,g; GL 667C b,c,d,e,f,g,h; HD10700 b,c,d,e,f
x= 2.5
n_rms = [ 2.12,2.12,2.92,2.92,2.92,1.4,1.4,1.4,1.4,1.4,1.4,1.4,x,x,x,x,x]
n_k = [ 1.3,1.29,0.84,1.09,0.95,3.93,0.61, 1.71,1.08,0.92,1.52,0.95,0.64,0.75,0.59,0.58,0.58]
n_yr = [2010,2010,2013,2013,2013,2013,2013,2013,2013,2013,2013,2013 ,2013,2013,2013,2013,2013]

# Plot the extra planets
ax.scatter(n_rms, n_k,c=n_yr,cmap=plt.get_cmap('rainbow'),vmin=min(date), vmax = max(date),s=100,alpha=0.8, edgecolor='none')

# Add some arrows to emphasize points on the diagram
ax.annotate(r'$\tau$ Ceti', xy=(2.5,0.65), xytext=(4.5, 0.59), fontsize=18, arrowprops=dict(facecolor='black', shrink=0.12, alpha=0.5) )
ax.annotate(r'$\alpha$ Cent', xy=(1.2,0.51), xytext=(2.0, 0.41), fontsize=18, arrowprops=dict(facecolor='black', shrink=0.12, alpha=0.5) )

# Plot things on a log scale
ax.set_xscale('log')
ax.set_yscale('log')

# Adjust the plot for better appearance
ax.set_xlim([0.4,300])
ax.set_ylim([0.4,300])

ax.set_ylabel('Velocity Semiamplitude [m/s]', fontsize=20)
ax.set_xlabel('RMS of Velocities [m/s]',fontsize=20)

ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)

cax2 = fig.add_axes([0.83,0.14,0.03,0.68])

cbar = plt.colorbar(cax, cax=cax2, orientation='vertical')
cbar.set_label('First Publication Date',fontsize=16)

plt.tight_layout()

# Save the product as a PDF
plt.savefig('rms_k_date.pdf')


