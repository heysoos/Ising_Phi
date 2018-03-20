import numpy as np
import scipy.io as sio
import pyphi
import time
from ising import gen_reservoir
import matplotlib.pyplot as plt

network = input('Network: ')
filename = 'Ising_Phi_' + network #'Ising_Phi_DMN' 
directory = 'reduced_networks/'
wd = 'Ising_output/'

mat = sio.loadmat(directory + wd + filename + '.mat')

T = mat['temp']
J = mat['J']
J = J!=0
spinBin = mat['spinBin']
#S = mat['EXPORT'] # spin time-series through temperature
M =  np.array(gen_reservoir(J.shape[1]),dtype='uint8')

TPM = mat['TPM']

templen = T.shape[1]
print('The number of data points to calculate Phi for is ' + str(templen))

tempstart = input('Start from data point: ')
tempstart = int(tempstart)

tempend = input('End at data point: ')
tempend = int(tempend)

increment = input('Increment every _ data points: ')
increment = int(increment)

suffix = input('Filename suffix: ')

numStates = M.shape[0]

ind = np.arange(tempstart,tempend,increment) # indices of data points that phi will be calculated for
T2 = T[0,ind]

looplen = ind.shape[0] # number of iterations of loop

# phi = np.zeros([numStates,templen])
phi = np.zeros([numStates,looplen])
phiSqr = np.zeros([numStates,looplen])
count = 0

print('Calculating...')
for temp in range(tempstart,tempend,increment):
  print( ((temp)/(tempend - tempstart))*100,"% Complete")
  for state in range(numStates): #numflips
    if spinBin[state,temp] != 0:
      start = time.time()
      #print("Starting state ", M[state,:], "at temp. ", T[0,temp])
      network = pyphi.Network(TPM[:,:,temp], connectivity_matrix=J)
      #subsystem = pyphi.Subsystem(network, S[:,state,temp], range(network.size))
      subsystem = pyphi.Subsystem(network, M[state,:], range(network.size))
      #print(subsystem)
      phi[state,count] = pyphi.compute.big_phi(subsystem)
      phiSqr[state,count] = phi[state,count]*phi[state,count]
      print("Phi = ", phi[state,count])
      #input()
      end = time.time()
      #print(end - start, "seconds elapsed")
  count += 1


phiSum = np.sum(phi*spinBin[:,ind],0)
phiSqrSum = np.sum(phiSqr*spinBin[:,ind],0)

phiSus = (phiSqrSum - phiSum*phiSum)/(T2*J.shape[0])
#print('Done!')

np.savez(directory + 'pyPhi/' + filename + '_' + suffix ,phiSum = phiSum,phiSus = phiSus,T2 = T2)

plt.plot(T2,phiSum,"o")
plt.plot(T2,phiSum,"o")
plt.show()
