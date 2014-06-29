#!/usr/bin/env python

# mass - density plot with relations plotted
# color code is equilibrium temperature from 0K - 3000K

import math

import numpy as np

from scipy.interpolate import pchip
from scipy import constants

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm


# Return the semi-major axis of a planet given period and star mass
# Works with the Kepler database, so we need the star name and planet name
def getA(per, mstar, pname, sname):
   a = []
   factor = 6.67e-8 / (4 * np.pi**2)
   for i, p in enumerate(per):
      starm = mstar[sname.index(int(pname[i]))]
      starm *= 1.9891e33
      a.append((factor*starm)**(1./3.) * (p*86400.)**(2./3.))
   a = np.array(a)
   a = a/1.4959e13
   return a 

filename = 'exoplanets.csv'
params = ["DENSITY", "MSINI", "MSINIUPPER", "MSINILOWER", "DENSITYUPPER", "DENSITYLOWER", "RSTAR", "TEFF", "A"]

f = open(filename, 'r')


# Parse out data from exoplanets database
delim = f.readline()
delimList = delim.split(',')

index = []

for param in params:
   index.append(delimList.index(param.strip()))

nameIndex = delimList.index('NAME')
eccIndex = delimList.index('ECC')

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
      if lineList[eccIndex].strip() == '':
         valueStore[-1].append(0.0)
      else:
         valueStore[-1].append(float(lineList[eccIndex].strip()))

f.close()

# Solar System properties
# Ordered as Venus-Earth-Jupiter-Saturn-Uranus-Neptune
solar_system_mass = [0.00256427115,0.00314635457,1.0,0.299409946,
                     0.0457346441,0.0539531012]
solar_system_density = [5.24,5.515,1.33,0.687,1.27,1.638]
solar_system_a = [0.723,1.0,5.204,9.582,19.201,30.047]
solar_system_e = [0.0067,0.0167,0.0489,0.0565,0.0457,0.0113]

# y holds Density
y = valueStore[0][:]
# x holds Msin(I)
x = valueStore[1][:]
m_errUp = valueStore[2][:]
m_errDown = valueStore[3][:]
d_errUp = valueStore[4][:]
d_errDown = valueStore[5][:]
rStar = valueStore[6][:]
tEff = valueStore[7][:]
A = valueStore[8][:]
e = valueStore[9][:]

# Add solar system properties to lists
for k in range(0,len(solar_system_mass)):
   y.append(solar_system_density[k])
   x.append(solar_system_mass[k])
   m_errUp.append(0.0)
   m_errDown.append(0.0)
   d_errUp.append(0.0)
   d_errDown.append(0.0)
   rStar.append(1.0)
   tEff.append(5780)
   A.append(solar_system_a[k])
   e.append(solar_system_e[k])

# Zero out errors for planets with a mass of 0.0 (To stop error bar marks on the plot)
for idx in range(0,len(x)):
   if x[idx] == 0.0:
      m_errUp[idx] = 0.0
      m_errDown[idx] = 0.0


# Calculate Equilibrium Temp
T_eq = []
for idx in range(0, len(x)):
   nu = math.sqrt(rStar[idx]*0.00464913034) * tEff[idx]
   de = (math.sqrt(2*A[idx])) * ((1-(e[idx]**2))**(1.0/8.0))
   T_eq.append(nu/de)


# Parse out information from the Kepler database
kep_mass = []
kep_density = []
kep_density_err = []
kep_mass_err = []
kep_per = []
kep_name = []
with open('new_kep_data','r') as f:
   for line in f:
      lines = line.split()
      if len(lines) < 23: continue
      try:
         kep_density.append(float(lines[9].strip()))
      except ValueError:
         continue
      else:
         kep_density_err.append(float(lines[11].strip()))
         kep_mass.append(float(lines[5].strip()))
         kep_mass_err.append(float(lines[7]))
         kep_name.append(float(lines[0]))
         kep_per.append(float(lines[1]))

kep_mstar = []
kep_teff = []
kep_starname = []
kep_rstar = []
with open('new_kep_star_data') as f:
   for line in f:
      lines = line.split() 
      kep_mstar.append(float(lines[13]))
      kep_teff.append(float(lines[4]))
      kep_starname.append(int(lines[0]))
      kep_rstar.append(float(lines[16]))

kep_mass = np.array(kep_mass)
kep_mass *= 0.003146
kep_mass_err = np.array(kep_mass_err)
kep_mass_err *= 0.003146
kep_density = np.array(kep_density)
kep_density_err = np.array(kep_density_err)


# Get semi-major axis of kepler planets
kep_a = getA(kep_per, kep_mstar, kep_name, kep_starname)

# Calc equilibrium temp for kepler planets
kep_teq = []
for i in range(0, len(kep_a)):
   nu = math.sqrt(kep_rstar[kep_starname.index(int(kep_name[i]))]*0.00464913034) * kep_teff[kep_starname.index(int(kep_name[i]))]
   de = math.sqrt(2*kep_a[i])
   kep_teq.append(nu/de)

kep_mass_l = []
for i in range(len(kep_mass)):
   if kep_mass[i] - kep_mass_err[i] < 0.001:
      kep_mass_l.append(kep_mass[i] - 0.001)
   else:
      kep_mass_l.append(kep_mass_err[i])

kep_density_l = []
for i in range(len(kep_density)):
   if kep_density[i] - kep_density_err[i] < 0.07:
      kep_density_l.append(kep_density[i] - 1e-2 )
   else:
      kep_density_l.append(kep_density_err[i])


# Plot the data calculated above
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

my_map = plt.get_cmap('spectral')


ax.scatter(kep_mass, kep_density, c=kep_teq,cmap=my_map, norm = colors.LogNorm(),vmin = 30.0, vmax = 3000.0, alpha=0.6,s=100,zorder=200)
ax.errorbar(kep_mass,kep_density,yerr=[kep_density_l,kep_density_err],xerr=[kep_mass_l,kep_mass_err],fmt=None, ecolor='gray', elinewidth=0.5,zorder=10)


cax = ax.scatter(x,y,c=T_eq,s=100.0,cmap=my_map, norm = colors.LogNorm(), alpha=0.5, vmin = 30.0, vmax = 3000.0, zorder=100)

ax.errorbar(x,y, yerr=[d_errDown, d_errUp], xerr=[m_errDown, m_errUp], fmt=None, ecolor='gray', elinewidth=0.5, zorder=0)

ax.set_xscale('log')
ax.set_yscale('log')

plt.xlim([0.001,30])
plt.ylim([0.07,60])
plt.xlabel(r'Planetary Mass [$M_{J}$]', fontsize=20)
plt.ylabel(r'Planetary Density [$g/{cm^3}$]', fontsize=20)

plt.grid(True)

# Add the vertical color bar with ticks at 1, 10, and 100 days
cbar = plt.colorbar(cax, ticks= [30,100,300,1000,3000])
cbar.ax.set_yticklabels(['30','100','300','1000','3000'])
cbar.set_label(r'$T_{EQ} [K]$', fontsize=18)

# Plot the mass radius relations
def toJupMass(inList):
   outlist = []
   for val in inList:
      outlist.append(val/317.828133)
   return outlist

# Takes a list of earth radius and converts to density
def toDensity(massList, radList):
   outlist = []
   for idx in range(0,len(massList)):
      mass = massList[idx]
      rad = radList[idx]
      g_mass = mass * 5.972e27
      cm_rad = rad * 6.371e8
      density = g_mass/( (4.0/3.0) * constants.pi * (cm_rad**3))
      outlist.append(density)
   return outlist

# Given values
MEarth = [0.010,0.032,0.10,0.32,1.00,3.16,10.0,31.6,100,316,1000.0]
MJup = [0.11,0.23,0.69,1.00,3.00,10.00]

IceR = [0.38,0.55,0.79,1.12,1.55,2.12,2.87,3.74,4.68,5.43,5.70]
RockR = [0.25,0.37,0.54,0.77,1.08,1.48,1.97,2.54,3.14,3.64,3.87]
IronR = [0.19,0.27,0.39,0.55,0.77,1.04,1.36,1.72,2.09,2.42,2.65]
H_HE = [0.89,0.95,1.03,1.05,1.07,1.08]

# List of Earth masses converted to M_Jup
x_new = toJupMass(MEarth)

# np.linspace(min,max,n) creates a list of n values between min and max
jupX = np.linspace(min(MJup),max(MJup),500)
earthX = np.linspace(min(x_new),max(x_new),500)

# convert radius and masses to a density
IceD = toDensity(MEarth, IceR)
RockD = toDensity(MEarth, RockR)
IronD = toDensity(MEarth, IronR)

lineColor = 'LightGray'

# Interpolate a line between the data points
interp = pchip(np.array(x_new),np.array(IceD))
ax.loglog(earthX, interp(earthX), c=lineColor, lw=4.0)
interp = pchip(np.array(x_new),np.array(RockD))
ax.loglog(earthX, interp(earthX), c=lineColor, lw=4.0)
interp = pchip(np.array(x_new),np.array(IronD))
ax.loglog(earthX, interp(earthX), c=lineColor, lw=4.0)

# Convert Jupiter radii and masses to density
H_HE_D = []
for idx in range(0,len(MJup)):
   g_mass = MJup[idx] * 1.898e30
   cm_rad = H_HE[idx] * 6.9911e9
   density = g_mass/( (4.0/3.0) * constants.pi * (cm_rad**3))
   H_HE_D.append(density)

interp = pchip(np.array(MJup),np.array(H_HE_D))
ax.loglog(jupX, interp(jupX), c=lineColor, lw=4.0)

text_color = 'LightSlateGray'

ax.text(0.002,9.0,'Iron', color=text_color, fontsize=20)
ax.text(0.002,3.3,'Rock', color=text_color, fontsize=20)
ax.text(0.002,1.05,'Ice', color=text_color, fontsize=20)
ax.text(0.06, 0.14, 'H-HE', color=text_color, fontsize=20)

plt.tick_params(axis='both', which='both', labelsize=18)

plt.gcf().set_size_inches(10,8)

# Save the resulting figure
plt.savefig('mass_density_teq.pdf',pad_size=0)



