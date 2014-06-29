#!/usr/bin/env python

# Mass ratio vs. Planet Period

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np

# List of masses for solar system bodies (kg)
# Order - Ganymead,Titan,Callisto,Io,Moon,Europa,Triton,Titania,Oberon,
#         Rhea, Iapetus, Dione, Oberon, Ariel, Umbriel, Miranda, 
#         Tethys, Enceladus,
#         Mimas, Janus, Epimetheus, Pandora, Prometheus

moon_masses = [1.4819e23/1.898e27,1.3452e23/568.3e24,
               1.075938e23/1.898e27,8.9319e22/1.898e27,
               7.3477e22/5.972e24,4.7998e22/1.898e27,
               2.14e22/102.4e24,3.526e21/8.6810e25,
               0.347e21/8.6810e25,23.1e20/568.3e24,
               15.9e20/568.3e24,10.52e20/568.3e24,
               30.14e20/8.6810e25, 13.53e20/8.6810e25,
               11.72e20/8.6810e25, 0.659e20/8.6810e25,
               6.22e20/568.3e24, 0.73e20/568.3e24,
               0.385e20/568.3e24, 0.0198e20/568.3e24,
               0.0055e20/568.3e24, 0.0013e20/568.3e24,
               0.0014e20/568.3e24]

# List of periods for solar system moons
moon_per = [7.15455296,15.945,16.6890184,1.769137786,
            27.321582,3.551181,5.876,8.706234854,13.463234,
            4.5175, 79.330183,2.7369,13.463,2.52,4.144,1.413,
            1.887802,1.370218,0.9424218,0.69459,0.69459,0.628804,0.612986]

# List of solar system bodies
# Order Jupiter, Saturn, Uranus, Neptune, Earth, Venus
planet_mass = [1.8986e27,5.6846e26,8.6810e25,1.0243e26,5.9736e24,4.8685e24]
planet_per = [4332.59,10759.22,30799.095,60190.03,365.256363004,224.698]

planet_ratio = [pmass/1.989e30 for pmass in planet_mass]

filename = 'exoplanets.csv'
params = ["MSINI", "MSTAR", "PER", "MSINIUPPER", "MSINILOWER", "ECC"]

f = open(filename, 'r')

delim = f.readline()
delimList = delim.split(',')

index = []

for param in params:
   index.append(delimList.index(param.strip()))

valueStore = []
for i in range(0,len(params)):
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

f.close()


mass = valueStore[0][:]
mstar = valueStore[1][:]
period = valueStore[2][:]

ecc = valueStore[5][:]

# Parse info for Pulsar, imaging, and microlensing planets
extra_massRatio = []
extra_per = []
extra_color = []
extra_names = []
with open('exoplanet.eu_catalog.csv','r') as f:
    line = f.readline()
    ls = line.split(',')
    l = [v.strip() for v in ls]
    mi = l.index('mass')
    si = l.index('star_mass')
    di = l.index('detection_type')
    pi = l.index('orbital_period')
    
    for line in f:
        ls = line.split(',')
        try:
            em = float(ls[mi].strip())
            ep = float(ls[pi].strip())
        except ValueError:
            continue
        else:
            if 'imaging' in ls[di]:
                try:
                    sm = float(ls[si].strip()) * 1047.92612
                except ValueError:
                    sm = 1047.92612
                extra_massRatio.append(em/sm)
                extra_per.append(ep)
                extra_color.append('black')
                extra_names.append(ls[0])
            elif 'pulsar' in ls[di]:
                try:
                    sm = float(ls[si].strip()) * 1047.92612
                except ValueError:
                    sm = 1.4 * 1047.92612
                extra_massRatio.append(em/sm)
                extra_per.append(ep)
                extra_color.append('purple')
                extra_names.append(ls[0])
            elif 'microlensing' in ls[di]:
                try:
                    sm = float(ls[si].strip()) * 1047.92612
                except ValueError:
                    sm = 0.1 * 1047.92612
                extra_massRatio.append(em/sm)
                extra_per.append(ep)
                extra_color.append('orange')
                extra_names.append(ls[0])
            else:
                continue
        

# We are only insterested in planets with periods 1 < per < 10,000
for i in range(len(extra_per)):
    if extra_per[i] <10000.0 and extra_per[i] > 1.0:
        if extra_massRatio[i] < 0.1 and extra_massRatio[i] > 1e-6:
            print extra_names[i], extra_color[i]

# Load data from the Kepler mission
f = open('cumulative.csv','r')
line=f.readline()
while line[0] == '#': line = f.readline()
delim_list = line.split(',')
kep_params = ['koi_period','koi_prad','koi_srad']

kep_index = []

for param in kep_params:
   kep_index.append(delim_list.index(param.strip()))

disp_ind = delim_list.index('koi_pdisposition')
name_ind = delim_list.index('kepoi_name')

valueStore2 = []
for i in range(0,len(kep_params)):
   valueStore2.append([])

valueStore2.append([])
valueStore2.append([])

# Kepler objects are named Star#.Planet#
# Here we are only interested in Kepler multiple planet systems
# So make sure that the star we plot has more than one planet
for line in f:
   lineList = line.split(',')
   flag = True
   for i in range(0,len(kep_params)):
      if lineList[kep_index[i]] == '':
         flag = False
         break
   if flag:
      for i in range(0,len(kep_params)):
         valueStore2[i].append(float(lineList[kep_index[i]].strip()))
      valueStore2[3].append(lineList[disp_ind].strip())
      valueStore2[4].append(lineList[name_ind].strip())

f.close()

kep_name = np.array(valueStore2[4][:])
ind = np.argsort(kep_name)
k = []
for indx in ind:
   obj_num = int(kep_name[indx][-2:])
   if obj_num > 1:
      k[-1] = True
      k.append(True)
   else: k.append(False)

k = np.array(k)

kep_per = np.array(valueStore2[0][:])
kep_rad = np.array(valueStore2[1][:])
kep_smass = np.array(valueStore2[2][:])
kep_disp = np.array(valueStore2[3][:])

kep_per = kep_per[ind][k]
kep_rad = kep_rad[ind][k]
kep_smass = kep_smass[ind][k]
kep_disp = kep_disp[ind][k]


# Convert values to plot kepler objects
u = []
v = []
w = []

cand_x = []
cand_y = []
ooi_x = []
ooi_y = []

# Kepler database includes candidates and Flase Positives
# We don't want to plot the false positives
for idx in range(len(kep_rad)):
 #  if kep_smass[idx] < 0.6:
      if kep_disp[idx] != 'FALSE POSITIVE':
         v.append((kep_rad[idx]**(2.06))/(kep_smass[idx]*333060.402))
         u.append(kep_per[idx])
         if kep_disp[idx] == 'CANDIDATE':
            cand_y.append((kep_rad[idx]**(2.06))/(kep_smass[idx]*333060.402))
            cand_x.append(kep_per[idx])
            w.append('DarkSlateGray')
         elif kep_disp[idx] == 'NOT DISPOSITIONED':
            ooi_y.append((kep_rad[idx]**(2.06))/(kep_smass[idx]*333060.402))
            ooi_x.append(kep_per[idx])
            w.append('LightSlateGray')
         else:
            w.append('SlateGray')
            print kep_disp[idx]

# Convert values to plot doppler objects
x = []
y = []
med_mass = []
med_ecc = []
for idx in range(len(mass)):
   #if mstar[idx] < 0.6:
   y.append(mass[idx]/(mstar[idx]*1047.92612))
   x.append(period[idx])
   if y[-1] > 0.0001 and period[idx] > 100. : 
      med_mass.append(mass[idx])
      med_ecc.append(ecc[idx])

# Set up the actual plots with values calculated above
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.scatter(x,y,c='Lime',label='Doppler Detections',
           edgecolor='none',s=81.0, alpha=0.75, zorder=100)
ax.scatter(moon_per,moon_masses,c='Red', label='Solar System Moons',
           edgecolor='none', s=81.0, alpha=0.8, zorder=150)
ax.scatter(planet_per, planet_ratio, c='MediumBlue', edgecolor='none',
           s=81.0, alpha=0.75, zorder=200)
ax.scatter(ooi_x,ooi_y,c='LightSlateGray',label='KOIs', edgecolor='none',
           s=81.0, alpha=0.75)
ax.scatter(cand_x,cand_y,c='LightSlateGray',label='Kepler Candidates',
           edgecolor='none', s=81.0, alpha=0.75)
ax.scatter(extra_per,extra_massRatio,c=extra_color,label="Others",edgecolor='none', s=81.0, alpha=0.75,zorder=130) 

ax.text(450,0.025, 'Eccentric Giants', color='Lime', fontsize=16.0)
ax.text(0.8,0.02, 'Hot Jupiters', color='Lime', fontsize=16.0)
ax.text(258,1.5e-6, 'V', color='MediumBlue', fontsize=15.0)
ax.text(445,1.8e-6, 'E', color='MediumBlue', fontsize=15.0)
ax.text(5500,0.001, 'J', color='MediumBlue', fontsize=15.0)
ax.text(257,1.3e-5, 'Ungiants', color='LightSlateGray', fontsize=15.0)
ax.text(0.32,0.0001, 'Solar System\n Satellites', color='Red', fontsize=15.0)

#ax.legend(loc='lower right', scatterpoints=1, 
#          scatteryoffsets=[0.5,0.5,0.5], prop={'size':12})
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.000001,0.1])
plt.xlim([0.3,10000])

plt.xlabel(r'Period $[Days]$', fontsize=20)
plt.ylabel(r'Mass Ratio $[{M}/{M_{Sun}}]$', fontsize=20)

# How many objects are getting plotted?
print "Doppler #", len(x)
print "Kep #", len(cand_x) + len(ooi_x)
moon_masses = np.array(moon_masses)
ind_y = (np.where(moon_masses > 0.000001,True,False))
ind_z = np.where( moon_masses < 0.1,True,False)
ind_v = np.logical_and(ind_y, ind_z)
count = 0
for ind in ind_v:
   if ind: count +=1

plt.tick_params(axis='both', which='both', labelsize=18)

plt.gcf().set_size_inches(10,8)

# Save the resulting figure
plt.savefig('MassRatio.pdf')


