#!/usr/bin/env python

# aitoff.py -- Generate a projection plot

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np


filename = 'listSorted.txt'

f = open(filename, 'r')

ra = []
dec = []
color = []

deg2rad = np.pi/180.0

# Last element in line is RA-DEC as a single string
# Second to last element is used for color coding
for line in f:
   lines = line.split()
   pos = lines[-1]
   color.append(float(lines[-2]))
   ra_hr = float(pos[1:3])
   ra_min = float(pos[3:5])
   ra_sec = float(pos[5:10])
   dec_deg = float(pos[10:13])
   dec_min = float(pos[13:15])
   dec_sec = float(pos[15:])
   ra_deg = ((ra_hr * 3600 + ra_min * 60 + ra_sec) *(360/ 86400.0))
   if ra_deg > 180.0:
      ra_deg = -(360.0 - ra_deg)
   ra.append(ra_deg * deg2rad)
   dec_deg = (dec_deg + dec_min/60.0 + dec_sec/3600.0)
   dec.append(dec_deg * deg2rad) 

f.close()

# Reverse the lists so that the 'important' points are plotted last (on top)
# This is a work around to not being able to pass zorder a list
ra.reverse()
dec.reverse()
color.reverse()

# Color map to use
myMap = plt.get_cmap('YlOrRd_r')


fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide', axisbg='LightBlue')
cax = ax.scatter(ra, dec,c=color,cmap=myMap, norm = colors.LogNorm(), vmin=0.4, vmax = max(color), edgecolor='none')

# Add the vertical color bar with ticks at 1, 10, and 100 days
cbar = plt.colorbar(cax, ticks= [1,2,4,12])
cbar.ax.set_yticklabels(['1','2','4','12'])
cbar.set_label(r'Unknown')

# Add labels and grid lines
plt.ylabel(r'Declination $[Degrees]$')
plt.figtext(0.5,0.2,r'Right Ascension $[Degrees]$', ha='center')
plt.grid(True)

# Save the result
plt.savefig("aitoff_plot.pdf")






