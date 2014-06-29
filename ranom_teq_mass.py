#!/usr/bin/env python

# ranom_teq_mass.py

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import math

filename = 'exoplanets.csv'
params = ["MSINI", "R", "RUPPER", "RLOWER", "MSINIUPPER", "MSINILOWER",
          "RSTAR", "TEFF", "A", "TRANSIT"]

f = open(filename, 'r')

delim = f.readline()
delimList = delim.split(',')

index = []

for param in params:
   index.append(delimList.index(param.strip()))

eIndex = delimList.index('ECC')

valueStore = []
for i in range(0,len(params)):
   valueStore.append([])
valueStore.append([])

for line in f:
   lineList = line.split(',')
   flag = True
   for i in range(0,len(params)):
      if lineList[index[i]] == '':
         flag = False
         break
   if flag:
      for i in range(0,len(params)):
         valueStore[i].append(float(lineList[index[i]].strip()))
      if lineList[eIndex].strip() == '':
         valueStore[-1].append(0.0)
      else:
         valueStore[-1].append(float(lineList[eIndex].strip()))

f.close()


mass = np.array(valueStore[0][:])
r = np.array(valueStore[1][:])
r_errUp = np.array(valueStore[2][:])
r_errDown = np.array(valueStore[3][:])
m_errUp = valueStore[4][:]
m_errDown = valueStore[5][:]
rStar = valueStore[6][:]
tEff = valueStore[7][:]
A = valueStore[8][:]
transit = valueStore[9][:]
e = valueStore[-1][:]


# Calculate Equilibrium Temps
T_eq = []
for idx in range(0, len(mass)):
   nu = math.sqrt(rStar[idx]*0.00464913034) * tEff[idx]
   de = (math.sqrt(2*A[idx])) * ((1-(e[idx]**2))**(1.0/8.0))
   T_eq.append(nu/de)

T_eq = np.array(T_eq)

# Zero out error terms if the mass if not availible from the csv file
for idx in range(0,len(mass)):
   if mass[idx] == 0.0:
      m_errUp[idx] = 0.0
      m_errDown[idx] = 0.0

# Calculate the radius anomily
r_anom = []
for idx in range(0,len(r)):
   if mass[idx] == 0.0:
      r_anom.append(0.0)
      continue
   else:
      l_m = math.log10(mass[idx])
      l_t = T_eq[idx]/1000.0
      anom = (1.08417 + (0.0940857*l_m) - 0.242831*(l_m**2) + 0.0947349*(l_m**3) + 0.0387851*l_t +
          0.00243981*l_m*l_t - 0.0244656*(l_m**2)*l_t + 0.0130659*(l_m**3)*(l_t) + 0.0240409*(l_t**2)
          - 0.0419296*l_m*(l_t**2) + 0.00693348*(l_m**2)*(l_t**2) + 0.00302157*(l_m**3)*(l_t**2))
      r_anom.append(r[idx] - anom)
   
r_anom = np.array(r_anom)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# Plot only Radial Velocity planets with mass above 100 M_Earth
for idx in range(0, len(r)):
   if transit[idx] != 0.0 and mass[idx] > 0.314:
      cax = ax.scatter(T_eq[idx],r_anom[idx],s=81.0,c=mass[idx],cmap=cm.gist_earth, 
             norm = colors.LogNorm(), alpha=0.5, vmin = 0.3, 
             vmax = 30.0, zorder=100)

jdx = np.where((transit != 0.0) & (mass > 0.314))

# Add error bars and set up plot visuals
ax.errorbar(T_eq[jdx],r_anom[jdx], yerr=(r_errDown[jdx], r_errUp[jdx]),
              xerr=None, fmt=None, ecolor='gray', zorder=5)

ax.plot([200.0,2700.0],[0.0,0.0],lw=2.0,c='gray')

plt.xlim([200.0,2700])
plt.ylim([-0.9,0.9])
plt.xlabel(r'$T_{EQ} [K]$',fontsize=20)
plt.ylabel(r'Radius Anomaly $[R_{\rm J}]$',fontsize=20)

plt.tick_params(which='both',labelsize=16)

plt.grid(True)

# Add the vertical color bar with ticks at 1, 10, and 100 days
cbar = plt.colorbar(cax, ticks= [0.3,3,30])
cbar.ax.set_yticklabels(['0.3','3','30'])
cbar.set_label(r'Mass $[M_{\rm J}]$',fontsize=18)

plt.gcf().set_size_inches(10,8)

# Save the resulting figure
plt.savefig('ranom_teq_mass.pdf',pad_size=0)

