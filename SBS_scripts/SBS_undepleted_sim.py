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

d_omega = np.arange(-4.4, 4.5, 0.1)

pid = list();
for Delta_omega in d_omega:
    settings["parameters"]["Delta_omega"] =  Delta_omega * settings["parameters"]["Gamma_B"]/2.0;
    settings["parameters"]["alpha"] = 0.0;
    settingsfile = 'settings.' + str(Delta_omega) + '.json';
    datafile = 'data.' + str(Delta_omega) + '.json';
    with open(settingsfile, 'w') as outfile:
        json.dump(settings, outfile)
    pid.append( subprocess.Popen(args=["SBS_Solver.exe", "-Settings", settingsfile, "-Data", datafile]) );
    if len(pid)>4:
        while len(pid)>0:
            pid.pop().wait()
    

gain = np.array([])
phase = np.array([])
for Delta_omega in d_omega:
    datafile = 'data.' + str(Delta_omega) + '.json';
    
    with open(datafile, "r") as fid:
        data=json.load(fid);

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
    
    gain = np.append(gain, np.absolute(Es[-1][-1])/np.absolute(Es[-1][0]) )
    #phase = np.append(phase, Es.imag[-1][-1] )
    phase = np.append(phase, np.angle(Es[-1][-1]) )

fig = plt.figure();

ax = fig.add_subplot(3,1,1)
plt.plot(d_omega, gain);
plt.grid(True, which='major')
plt.grid(True, which='minor')
plt.minorticks_on()
plt.xlabel('$2 \delta \Omega / \Gamma_B$')
plt.ylabel('Gain')
plt.title('SBS resonance')

ax = fig.add_subplot(3,1,2)
phase = -phase;
plt.plot(d_omega, phase);
plt.grid(True, which='major')
plt.grid(True, which='minor')
plt.minorticks_on()
plt.xlabel('$2 \delta \Omega / \Gamma_B$')
plt.ylabel('Phase')
#plt.title('SBS resonance')

group_index=np.diff(phase);
ax = fig.add_subplot(3,1,3)
plt.plot(d_omega[1:], group_index);
plt.grid(True, which='major')
plt.grid(True, which='minor')
plt.minorticks_on()
plt.xlabel('$2 \delta \Omega / \Gamma_B$')
plt.ylabel('Group Index')
#plt.title('SBS resonance')

fig = plt.figure();

ax = fig.add_subplot(1,1,1)
plt.plot(gain[1:], time_delay);

plt.show();

plt.savefig('sbs_resonance.png')

"""
filename = 'data.2.5.json';
fid=open(filename, "r");
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
