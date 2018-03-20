import numpy as np
import scipy.io as sio
# import pyphi
import time
from ising import gen_reservoir
import matplotlib.pyplot as plt

network = input('Network: ')
method = input('Method (met/glaub): ')
filename = 'Ising_Networks/' + network + '/' + method + '/Ising_met_1_meanPhi'
directory = 'Parallel Code Stuff/Simulations/'
#wd = 'Ising_output/'

#mat = sio.loadmat(directory + wd + filename + '.mat')

#T = mat['temp']

#phiSum = np.load(directory + 'pyPhi/' + filename + '.npy')

loadfile = np.load(directory + 'pyPhi/' + filename + '.npz')
phiSum = loadfile['phiSum']
phiSus = loadfile['phiSus']
T2 = loadfile['T2']

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
  
f0 = movingaverage(phiSum,10)
f1 = movingaverage(phiSus*T2,10)
f2 = movingaverage(phiSum*phiSus,10)

phiplot = f0/np.amax(f0)
phiSusplot = f1/np.amax(f1)
phiphiSusplot = f2/(np.amax(f2))

# smoothed and normalized data

# plt.plot(T2,phiplot,"o",label = 'Phi')
# plt.plot(T2,phiSusplot,"o",label = 'Phi Susceptibility')
# plt.plot(T2,phiphiSusplot,"o",label = 'Phi*(Phi Susceptibility)')

# Raw data
plt.plot(T2,phiSum,"o",label = 'Phi')
#plt.plot(T2,phiSus,"o",label = 'Phi Susceptibility')


plt.xlabel('Temperature')
plt.title(network + ' Phi Properties')
plt.xticks(np.arange(min(T2), max(T2), 0.5,))
plt.minorticks_on()
plt.legend()
plt.ylim((0, 1.1))
plt.show()