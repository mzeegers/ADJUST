# Script to compute various attenuation matrices (including dictionaries)
# A separate script (with small differences) is provided for the thorax phantom
# Libraries required: csv, matplotlib, numpy, physdata and scipy
# physdata can be downloaded here: https://pypi.org/project/physdata

# Output (saved to file):
#   T - spectral distribution of dictionary materials (size: t x c)

# Here, t is the number of materials in the dictionary, c is the number of
# spectral bins

# Authors:
#   Ajinkya Kadu,
#       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
#   Mathé Zeegers, 
#       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

import csv
import matplotlib.pyplot as plt
import numpy as np
import operator
import os
import physdata.xray
import scipy.interpolate
from ElementaryData import *


### Options ###
energies = 100
MinEnergy = 20
MaxEnergy = 119

#Materials in the X-ray dictionary matrix T (the second item of each entry specifies the material)
materialsdict = [[1,1, None], [2, 2, None], [3,3,None], [4, 4, None],
             [5,5, None], [6, 6, None], [7,7,None], [8, 8, None],
             [9,9, None], [10, 10, None], [11,11,None], [12, 12, None],
             [13,13, None], [14, 14, None], [15,15,None], [16, 16, None],
             [17,17, None], [18, 18, None], [19,19,None], [20, 20, None],
             [21,21, None], [22, 22, None], [23,23,None], [24, 24, None],
             [25,25, None], [26, 26, None], [27,27,None], [28, 28, None],
             [29,29, None], [30, 30, None], [31,31,None], [32, 32, None],
             [33,33, None], [34, 34, None], [35,35,None], [36, 36, None],
             [37,37, None], [38, 38, None], [39,39,None], [40, 40, None],
             [41,41, None], [42, 42, None], [43,43,None], [44, 44, None],
             [45,45, None], [46, 46, None], [47,47,None], [48, 48, None],
             [49,49, None], [50, 50, None], [51,51,None], [52, 52, None],
             [53,53, None], [54, 54, None], [55,55,None], [56, 56, None],
             [57,57, None], [58, 58, None], [59,59,None], [60, 60, None],
             [61,61, None], [62, 62, None], [63,63,None], [64, 64, None],
             [65,65, None], [66, 66, None], [67,67,None], [68, 68, None],
             [69,69, None], [70, 70, None], [71,71,None], [72, 72, None],
             [73,73, None], [74, 74, None], [75,75,None], [76, 76, None],
             [77,77, None], [78, 78, None], [79,79,None], [80, 80, None],
             [81,81, None], [82, 82, None], [83,83,None], [84, 84, None],
             [85,'water', None], [86, 'tissue', None], [87,'bone',None], [88, 'polyethylene', None], 
             [89,'pyrex', None], [90,'muscle',None], [91, 'telluride', None], 
             [92,'gallium', None], [93, 'air', None]]

###################
### Compute T   ###
###################

###### Supporting functions for fetching X-ray attenuation coefficients ##################

def computeMaterialEnergyMatrixXray(materials):

    #Collect the attenuation spectra for each material
    AttenuationSpectra = collectAttenuationSpectra(RootDataPath, materials)
    
    #Setup the energy discretization
    EnergyBounds = np.linspace(MinEnergy, MaxEnergy, num = energies)
    print(energies, "energies with energy levels:", EnergyBounds)

    #Get the attenuations at the energy grid points for each material
    energyatts = np.zeros((len(materials), energies)) #energy attenuation for each material per channel
    for idxe, e in enumerate(EnergyBounds):
        for idxAtt, Att in enumerate(AttenuationSpectra):
            energyatts[idxAtt,idxe] = Att[4](e)

    return energyatts

def collectAttenuationSpectra(Rootdatapath, Labels):
    AttenuationSpectra = []
    print("\nLoading attenuation spectra...")
            
    for mat in Labels:
        if(mat[0] != 0 and mat[1] != "Void"): #Exclude Voids                
            if (mat[1] in [i[1] for i in ElementaryData]):                      #Elementary material
                AtNo = elementToAtomic(mat[1])
                if (AtNo > 0):
                    attData = getAttenuationSpectrum(AtNo, Rootdatapath)
                    AttenuationSpectra.append((mat[0],)+(mat[1],) + attData)
            else:
                attData = getAttenuationSpectrum(mat[1], Rootdatapath)          #Mixture material
                AttenuationSpectra.append((mat[0],)+(mat[1],) + attData)
        elif(mat[0] != 0 and mat[1] == "Void"):
            #Make the zero attenuation spectrum
            x, y = np.arange(0,200), np.zeros(200)
            spectrum = scipy.interpolate.interp1d(x, y)
            AttenuationSpectra.append((mat[0],)+('Void',) + (x,y,spectrum))

    AttenuationSpectra.sort(key = operator.itemgetter(0)) #Keep sorted on identifier 
    print("Attenuation spectra fully loaded")

    return AttenuationSpectra

def elementToAtomic(materialname):
    Candidates = [x for x in ElementaryData if x[1] == materialname]
    if not Candidates:
        return 0
    else:
        return next(x for x in ElementaryData if x[1] == materialname)[0]

def getAttenuationSpectrum(materialno, rootdatapath):
 
    data = np.array(physdata.xray.fetch_coefficients(materialno))

    data[:, 0] *= 1000  #convert from MV to kV
       
    #Detect multiple energy values, to prevent strange interpolation at K-edges
    for i in range(0,len(data[:,0])-1):
        if(data[i,0] == data[i+1,0]):
            data[i+1,0] += 0.00000001

    #Create (linear) interpolation function
    x, y = data[:, 0], data[:, 1]
    spectrum = scipy.interpolate.interp1d(x, y)

    return data[:,0], data[:,1], spectrum
    
###################################################################################################

### Compute T
T = computeMaterialEnergyMatrixXray(materialsdict)
T = np.asarray(T)
    
### Show matrix sizes, ranks and other properties
print("\nMatrix rank checks and other properties:")
print("Shape of T:", T.shape)
print("Rank T:", np.linalg.matrix_rank(T))

### Plot spectra in dictionary T 
plt.figure() 
for i in range(0, T.shape[0]):
    plt.plot(np.linspace(MinEnergy, MaxEnergy, num = energies), T[i,:], label=str(i+1))
plt.ylabel('Attenuation')
plt.xlabel('Bins')
plt.xlim(left = MinEnergy, right = MaxEnergy)
plt.ylim(bottom = 0, top = 5)
plt.title("Spectra")
plt.legend(title = 'Material')
plt.show()

### Save the matrices
print("\n Saving the spectrum...")
np.savetxt("../Data/MatrixTXray" + str(energies) + "chan_" + str(len(materialsdict)) + "mat_MinEn" + str(MinEnergy) + "_MaxEn" + str(MaxEnergy) + "V2.csv", T, delimiter=",")
