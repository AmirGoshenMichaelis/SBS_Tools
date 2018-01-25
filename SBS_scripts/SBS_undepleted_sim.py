# in# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import json
import os
import subprocess
import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import fftpack as fft

os.chdir("C:/Users/Asus/Documents/Research/Avishay/SBS_undepleted_sim")
fid=open("settings.json", "r");
settings=json.load(fid)
fid.close();

pid = list();
for Delta_omega in np.arange(-4.5, 4.5, 1.0):
    settings["parameters"]["Delta_omega"] =  Delta_omega * settings["parameters"]["Gamma_B"]/2.0;
    settingsfile = 'settings.' + str(Delta_omega) + '.json';
    datafile = 'data.' + str(Delta_omega) + '.json';
    with open(settingsfile, 'w') as outfile:
        json.dump(settings, outfile)
    pid.append( subprocess.Popen(args=["SBS_Solver.exe", "-Settings", settingsfile, "-Data", datafile]) );
    if len(pid)>4:
        while len(pid)>0:
            pid.pop().wait()
    

"""
fid=open("data.json", "r");
data=json.load(fid)
fid.close();

a = data.get('Ep').get('real');
b = data.get('Ep').get('imag');
Ep=np.array(a)+np.array(b)*1.0j;

a = data.get('Es').get('real');
b = data.get('Es').get('imag');
Es=np.array(a)+np.array(b)*1.0j;

a = data.get('Rho').get('real');
b = data.get('Rho').get('imag');
rho = np.array(a)+np.array(b)*1.0j;

t = np.array(data.get('t'));
z = np.array(data.get('z'));

Es__ =fft.fftn(Es)
Ep__ =fft.fftn(Ep)
rho__=fft.fftn(rho)


plt.plot(z[:], np.abs(Es[-10])**2);
plt.show();

plt.plot(z[:],Es[-3]*Es[-3].conj());
plt.show();

%matplotlib auto
fig=plt.figure()
#ax=fig.gca(projection='3d')
ax=fig.add_subplot(111,projection='3d')
T,Z=np.meshgrid(t,z)
ax.plot_wireframe(T,Z,Es__.imag)
ax.plot_wireframe(T,Z,np.abs(Es))
ax.set_xlabel('T')
ax.set_ylabel('Z')
plt.show()

"""
