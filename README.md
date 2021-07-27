SH-wave simulation code for ambient noise polarization 
=======
The ambient noise modelling codes developed for the paper *Fluid migrations and volcanic earthquakes from depolarized ambient noise*, submitted to Nature: Communications and available as pre-print on [Research Square](https://www.researchsquare.com/article/rs-470597/v1)

Installation
------------
The simulation works in Matlab version 2019a and should not require any additional toolbox.

**Use** 
------------
1) Download the package.
2) Enter folder Modelling.
3) Enter the folder corresponding to the source distribution you want to simulate (FarLine or NearCoastline).
4) Enter the folder corresponding to either the homogeneous or heterogeneous case.
5) Run the file *CF_noise_sources_?.m* where *?* corresponds either to *Line* or *CircleN*.

**Outputs**
------------
The simulation will produce a video showing noise propagation, a figure showing the variations of shear modulus, and a .mat file containig the simulated noise data at all stations considered.

**Simulation Time**
------------
Each simulation takes about ?? on a Macbook Pro, 2,3 GHz, 8-Core Intel Core i9 with 32 Gb of memory.

**Reproduction Instructions**
------------
To reproduce the panels in Extended Data Fig. 8a,b, use the codes in Folder Modellin/FigureED8. This contains the results of the Data Processing folder (available in the in the [Open Science Framework](https://osf.io/kqtbp/)). The codes were applied to the simulated noise data in .mat format. 
