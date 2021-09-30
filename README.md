[![DOI](https://zenodo.org/badge/388130635.svg)](https://zenodo.org/badge/latestdoi/388130635)

SH-wave simulation code for ambient noise polarization
=======
The ambient noise modelling codes developed for the paper *Fluid migrations and volcanic earthquakes from depolarized ambient noise*, submitted to Nature: Communications and available as pre-print on [Research Square](https://www.researchsquare.com/article/rs-470597/v1)

Content and requirements
------------
The simulation works in Matlab version 2019a and requires the Image Processing, Mapping and Wavelet toolboxes.

The folder contains four subfolders:

**FarLine** to obtain results from a line of ambient noise sources in the far field (plane wave).

**NearCoastline** to obtain results from a circle of ambient noise sources in the near field (circular wave).

**FigureED8** to obtain the corresponding Figure.

**idw** with the corresponding function and licence.


The **FarLine** and **NearCoastline** folders are used for modelling, inside each of them there are two more subfolders with the names:

**Homogeneous** to obtain results for an homogeneous model (constant shear modulus).

**Heterogeneous** to obtain results for an heterogeneous model (shear modulus that varies depending on polarisation results).



**Use**
------------
1) Download the package.
2) Enter the folder Modelling.
3) Enter the folder corresponding to the source distribution you want to simulate (FarLine or NearCoastline).
4) Enter the folder corresponding to either the homogeneous or heterogeneous case.
5) Run the file *CF_noise_sources_?.m* where *?* corresponds either to *Line* or *CircleN*.

**Outputs**
------------
The simulation will produce a video showing noise propagation (.avi), a figure showing the variations of shear modulus (.tiff), and a .mat file that contains the simulated noise data at all stations considered (.mat) and a figure representing the corresponding seismograms (.fig).

**Simulation Time**
------------
Each simulation takes about 2 hours and 15 minutes on a Macbook Pro, 2,3 GHz, 8-Core Intel Core i9 with 32 Gb of memory.

**Reproduction Instructions**
------------
To reproduce the panels in Extended Data Fig. 8a,b, use the codes in Folder Modellin/FigureED8. This contains the results of the Data Processing folder (available in the in the [Open Science Framework](https://osf.io/kqtbp/)). The codes were applied to the simulated noise data in .mat format.
