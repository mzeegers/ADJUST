# Script to compute various source spectra
# Libraries required: matplotlib, numpy, spekpy 

# Output (saved to file):
#   I - source spectrum vector (size: c)

# Here, c is the number of spectral bins

# Authors:
#   Ajinkya Kadu,
#       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
#   Mathé Zeegers, 
#       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

import matplotlib.pyplot as plt # Import library for plotting
import numpy as np
import spekpy as sp # Import SpekPy

#Settings for the Thorax phantom
energies = 100
minEnergy = 20
maxEnergy = 80
material = 'W'

#Settings for all other phantoms
energies = 100
minEnergy = 5
maxEnergy = 35
material = 'Mo'

# Generate the spectrum
s = sp.Spek(kvp = maxEnergy, dk = (maxEnergy - minEnergy)/energies, targ = material)

# Get the spectrum (as a collection of points)
k, f = s.get_spectrum(edges=True) 

# Get all points higher than the minimum energy
kreal = [k[i] for i in range(0,len(k)) if k[i] >= minEnergy]
freal = [f[i] for i in range(0,len(k)) if k[i] >= minEnergy]

# Convert to average midpoint values in each bin
kreal = [(x + y)/2 for x, y in zip(kreal[1::2], kreal[2::2])]
freal = [(x + y)/2 for x, y in zip(freal[1::2], freal[2::2])]

# Plot the spectrum
plt.plot(kreal, freal)
plt.xlabel('Energy [keV]')
plt.ylabel('Photon flux per unit energy [photons/keV]')
if(material == 'W'):
    plt.title('Thungsten source spectrum \nwith 80kVp tube voltage')
else:
    plt.title('Molybdenum source spectrum \nwith 35kVp tube voltage')
plt.ylim(bottom = 0)
plt.xlim(left = minEnergy, right = maxEnergy)
plt.show()

#Write the spectrum to output csv file
freal =  np.asarray(freal)
np.savetxt("../Data/VectorIXrayMat" + str(material) + "En" + str(energies) + "Min" + str(minEnergy) + "Max" + str(maxEnergy) + ".csv", freal, delimiter=",")
