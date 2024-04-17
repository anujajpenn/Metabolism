# Model for Cell Metabolic Potential
#### Contained within this repository are two finite element analysis models used to implement the free-energy based approach to calculating the metabolic potential of a cell as outlined in the cited manuscript. The two models are separated based on the environment they recapitulate (2D vs 3D), but are built based upon the same physical implementation and running procedure.
### This repository relates to the research publication "Jaganathan A. et al., Mechano-metabolism of adherent cells in 2D and 3D microenvironments, BioRxiv, 2024" (https: ). Please cite our paper if this code is used or modified for further use.
## System Requirements
### Operating system:
#### This package is supported for macOS and windows. The package has been tested on the following systems:
#### •	macOS Monterey 12.4, Processor: 2.3 GHz 8-Core Intel i9, RAM: 32 GB
#### •	Windows 11 Home, Processor; 2.10GHz Hexa-Core AMD Ryzen 5 5500U with Radeon Graphics, RAM:16 GB
### MATLAB dependencies
#### This package requires MATLAB LiveLink for COMOSL Multiphysics. This package has been tested on the following software versions:
#### •	MATLAB_R2023b
#### •	COMSOL Multiphysics 6.0
## Installation guide
### 1.	Install MATLAB (https://www.mathworks.com/products/matlab.html). The installation usually takes approximately 1 hour.
### 2.	Install COMSOL Multiphysics (https://www.comsol.com/) along with the Structural Mechanics Module and LiveLink for MATLAB. During installation, specify the installation directory of MATLAB on your system when prompted with the option to establish the LiveLink connection. 
### 3.	Launch MATLAB using the LiveLink interface.
## Demo
### 1.	Ensure all downloaded files from either the 2D or 3D model in this repository are placed in the same folder. Open “cell_3D.m.” to run the 3D simulation and similarly, open "cell_2D.m" to run the 2D simulation.
## Steps for running the 3D files
### 2.	Click “Run”. By default, the code is setup to run the model for the full range of cell aspect ratios and ECM density tested in the study, but this list may be adjusted to the desired values. The code will take around 5~10 seconds to finish computing and compiling data for one combination of cell shape and ECM stiffness, and around 1 hour in total for all combinations.
### 3.	Output will be stored in the MATLAB workspace. Open “Post_calculations.m” and click “Run” with the output of “cell_3D.m” opened in the workspace. 
### 4.	All data output from the model will be stored in the workspace variables and displayed in the resulting figures. 
## Steps for running the 2D files
### 2. Set ECM stiffness to desired value and click "Run". Default ranges for the two aspect ratios are already setup in the code file and can be altered to desired values. The code will take 20-30 minutes to finish computing the stresses, contractility and energies for the given range of aspect ratios corresponding to a single ECM stiffness level. The computation can be repeated for as many ECM stiffness levels as desired.
### 3. Output is generated in Probe Table 29. Using the 'Export Table' feature in COMSOL Multiphysics, the table can be stored as an Excel file (see for example 2d10.xlsx). 
### 4. Similarly, output files can be generated and stored for different ECM stiffness levels in the range of 1-200 kPa.
## Sample Output 3D
### 1.	Sample MATLAB workspaces produced by the programs for the 3D model are provided by “cell_3D_output.mat” and “Post_calculations_output.mat.”
### 2.	A sample MATLAB figure showing the calculated metabolic potential and optimum cell shape for various combinations of cell aspect ratio and collagen density is provided in "EnergyLandscape.fig.”
## Sample Output 2D
### 1. MATLAB script file "contours2D.m" can be used in conjuction with the generated sample data file "2d10.xlsx" to plot the metabolic potential contours for gel stiffness of 10kPa.
### 2. A sample MATLAB figure showing the metabolic potential landscape with the optimum cell shape identified is provided in "2Dcontour_10kPa.fig".

## Instruction for use
#### The general structure provided by these scripts and models can be utilized to study a wider or different range of cell aspect ratio or matrix stiffness than analyzed in the corresponding manuscript. The values of interest should be specified in the “AR” and “Density” variables at the beginning of the script. Further, different cell geometries or matrix material models can be specified by modifying the COMSOL files. 
