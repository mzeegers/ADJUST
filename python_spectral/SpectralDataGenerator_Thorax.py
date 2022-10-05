# Script to compute various attenuation matrices (including dictionaries)
# This is a separate script (with small differences) for the thorax phantom for which the F matrix is necessary
# Libraries required: csv, matplotlib, numpy, physdata and scipy
# physdata can be downloaded here: https://pypi.org/project/physdata

# Output (saved to file):
#   T - spectral distribution of dictionary materials (size: t x c)
#   F - spectral distribution of materials (size: k x c)

# Here, t is the number of materials in the dictionary, c is the number of
# spectral bins, k is the number of materials in the phantom

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
MaxEnergy = 80

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
             [49,49, None], [50, 50, None],
             [51,'water', None], [52, 'tissue', None], [53,'bone',None], [54, 'polyethylene', None], 
             [55,'pyrex', None], [56, 'glass', None], [57,'muscle',None], [58, 'telluride', None], 
             [59,'gallium', None], [60, 'air', None]]

#Materials in the X-ray material matrix F 
materials = [[1,'bone', None], [2, 'blood', None], [3,'tissue',None], [4, 'blood', None], [5,'lung', None], [6,53,None]]

#########################
### Compute F and T   ###
#########################

###### Supporting functions for fetching X-ray attenuation coefficients ##################

def computeMaterialEnergyMatrixXray(materials):
    #Collect the attenuation spectra for each material
    AttenuationSpectra = collectAttenuationSpectra(materials)
    
    #Setup the energy discretization
    EnergyBounds = np.linspace(MinEnergy, MaxEnergy, num = energies)
    print(energies, "energies with energy levels:", EnergyBounds)

    #Get the attenuations at the energy grid points for each material
    energyatts = np.zeros((len(materials), energies)) #energy attenuation for each material per channel
    for idxe, e in enumerate(EnergyBounds):
        for idxAtt, Att in enumerate(AttenuationSpectra):
            energyatts[idxAtt,idxe] = Att[4](e)

    return energyatts

def collectAttenuationSpectra(Labels):
    AttenuationSpectra = []
    print("\nLoading attenuation spectra...")
            
    for mat in Labels:
        if(mat[0] != 0 and mat[1] != "Void"): #Exclude Voids                
            if (mat[1] in [i[1] for i in ElementaryData]):                      #Elementary material
                AtNo = elementToAtomic(mat[1])
                if (AtNo > 0):
                    attData = getAttenuationSpectrum(AtNo)
                    AttenuationSpectra.append((mat[0],)+(mat[1],) + attData)
            else:
                attData = getAttenuationSpectrum(mat[1])                        #Mixture material
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

def getAttenuationSpectrum(materialno):

    data = np.array(physdata.xray.fetch_coefficients(materialno)) #Density taken from array

    data[:, 0] *= 1000  #convert from MV to kV
       
    #Detect multiple values
    for i in range(0,len(data[:,0])-1):
        if(data[i,0] == data[i+1,0]):
            data[i+1,0] += 0.00000001

    #Create (linear) interpolation function
    x, y = data[:, 0], data[:, 1]
    spectrum = scipy.interpolate.interp1d(x, y)

    return data[:,0], data[:,1], spectrum

###################################################################################################

### Compute F and T
F = computeMaterialEnergyMatrixXray(materials)
T = computeMaterialEnergyMatrixXray(materialsdict)
        
F = np.asarray(F)
T = np.asarray(T)
    
#Compute 90% blood and 10% iodine in the second row of F
F[1,:] = F[1,:]*0.9 + F[5,:]*0.1

### Show matrix sizes, ranks and other properties
print("\nShape of F:", F.shape)
print("\nShape of T:", T.shape)
print("\nMatrix rank checks and other properties:")
print("Rank F:", np.linalg.matrix_rank(F))
print("Rank T:", np.linalg.matrix_rank(T))
print("Minimum of F", np.min(F))
print("Maximum of F", np.max(F))
print("Condition number of F", np.linalg.cond(F))
print("Eigenvalues of F:")
print(np.linalg.svd(F)[1])

### Plot spectra in F
plt.figure()    
for i in range(0, F.shape[0]):
    plt.plot(np.linspace(MinEnergy, MaxEnergy, num = energies), F[i,:], label=str(materials[i][1]))
plt.ylabel('Attenuation')
plt.xlabel('Bins')
plt.xlim(left = MinEnergy, right = MaxEnergy)
plt.ylim(bottom= 0, top = 5)
plt.title("Spectra")
plt.legend(title = 'Material')
pq.image(F, title="F")
plt.show()

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
print("Saving...")
np.savetxt("../Data/Thorax_MatrixFXray" + str(energies) + "chan_" + str(len(materials)) + "matBoneBloodIodineSofttissueBloodLungsIodine_MinEn" + str(MinEnergy) + "_MaxEn" + str(MaxEnergy) + "V2.csv", F, delimiter=",")
np.savetxt("../Data/Thorax_MatrixTXray" + str(energies) + "chan_" + str(len(materialsdict)) + "mat_MinEn" + str(MinEnergy) + "_MaxEn" + str(MaxEnergy) + "V2.csv", T, delimiter=",")
