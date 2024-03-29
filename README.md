# <a href="https://arxiv.org/abs/2112.11406" style="color: black;">ADJUST</a>: A Dictionary-based Joint reconstruction and Unmixing method for Spectral Tomography

   <p align="center">
   <img src="./images/UnmixingSchemev9_Github.svg">
    </p>

## Introduction

ADJUST is a MATLAB solver for large scale spectral 
tomographic inverse problems where a dictionary of spectral response
of various materials is available. This algorithm jointly reconstructs and unmixes (i.e., material decomposition) spectral tomographic measurements. In this work package the following algorithms are included:

1. Unmixing-then-reconstruction (UR), Reconstruction-then-Unmixing (RU) and classic Joint reconstruction algorithm (cJoint) solves the following problem:  
   <p align="center">
   <img src="./images/eq1.svg">
    </p>

2. ADJUST (dictonary-based method):  
   <p align="center">
   <img src="./images/eq2.svg">
   </p>

The matrix W is a tomography operator of size m x n, 
Y are the tomographic measurements of size m x c, 
A contains the spatial information of materials (size: n x k), 
F contains the spectral information of materials (size: k x c),
T is a spectral dictionary of p materials, while R is dictionary coefficient matrix.
Here, m are the number of tomographic measurements, n is the size of image,
c are the spectral channels, and k are the number of materials.


## Requirements

To run the examples on MATLAB (r2019 and above recommended), please install the following packages:

1. ASTRA Toolbox: 
https://github.com/astra-toolbox/astra-toolbox
2. SPOT operator:
https://github.com/mpf/spot
3. MinConf package: 
https://www.cs.ubc.ca/~schmidtm/Software/minConf.html
4. 3D Shepp-Logan phantom package (for the 3D example only): 
https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom


## Scripts:

Eight examples are provided with the code. Before running any example, please run the startup.m script first.
1. example1_SheppLogan: 
The classic Shepp-Logan phantom with variable size
2. example2_Disk: 
A custom-made disk example with variable size and number of disks
3. example3_Thorax: 
The synthetic phantom mimics internal body structure in an abstract fashion (taken from the CONRAD software framework: https://git5.cs.fau.de/PyConrad/pyCONRAD)
4. example4_SheppLogan_SparseAngle: 
The Shepp-Logan phantom with only 10 projection angles (adjustable)
5. example5_SheppLogan_LimitedView: 
The Shepp-Logan phantom with a missing wedge of 60 degrees (adjustable)
6. example6_MixedDisk: 
A custom-made disk example with variable size and number of disks, where the disks consists of multiple materials
7. example7_SheppLogan3D: 
3D version of the Shepp-Logan phnatom with adjustable size
8. example8_microCTData.m:
A demonstration of the approach on real micro-CT data from a natural rock sample. The source of the data can be found further below.

On each case the four implemented algorithms (UR, RU, cJoint and ADJUST) can be tested. The perfomance of these algorithms are captured by computing and presenting the MSE, PSNR, and the SSIM between the reconstructed material spatial maps and the ground truth. 


## Example results:

As an example, below are results for a numerical study on the Shepp-Logan phantom:
   <p align="center">
   <img src="./images/comparison.png">
   </p>
   <p align="center">
   <img src="./images/comparison-maps.png">
   </p>
   
For the micro-CT dataset, four different materials (lead, tungsten, quartz and gold) can be reconstructed. An example reconstruction is shown below:
   
   <p align="center">
   <img src="./images/micro_ct_dataset_result.png">
   </p>


## Generation of spectral data

To reproduce the provided material spectra matrices (F and T), the scripts for generating these are provided in the python_spectral folder.
These scripts are written in Python (version 3.7 or higher). Apart from this, physdata, csv and scipy are necessary to run these scripts.


## References

The algorithms implemented in this MATLAB package are described in following [paper](https://iopscience.iop.org/article/10.1088/1361-6420/ac932e). If you use (parts of) this code in a publication, we would appreciate it if you would refer to:

```
@article{
  title={ADJUST: A Dictionary-Based Joint Reconstruction and Unmixing Method for Spectral Tomography},
  author={Zeegers, Math{\'e} T and Kadu, Ajinkya and van Leeuwen, Tristan and Batenburg, Kees Joost},
  journal={Inverse Problems},
  year={2022},
  volume={22},
  number={12},
  pages={125002},
  publisher={IOP Publishing},
  doi={10.1088/1361-6420/ac932e}
}
```
The preprint can be found [here](https://arxiv.org/abs/2112.11406).

The micro-CT dataset required to run the eighth script is introduced in this [paper](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/xrs.3200) by Sittner et al. (2020). The dataset can be found [here](https://rodare.hzdr.de/record/1627).


## Authors

Code written by:
- Ajinkya Kadu (aak [at] cwi [dot] nl)
- Mathé Zeegers (m [dot] t [dot] zeegers [at] cwi [dot] nl).
