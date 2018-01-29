# in# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import json
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import fftpack as fft

os.chdir("C:/Users/Asus/Documents/Research/Avishay/SBS_Tools/SBS_Solver/")
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

