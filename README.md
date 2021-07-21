# SH-Ambient-Noise-NatComm
The ambient noise modelling codes developed for the paper submitted to Nature: Communications. The simulation works ith Matlab version 2019a and should not require any additional toolbox.

**Use** 
1) Enter folder Modelling.
2) Enter the folder corresponding to the source distribution you want to simulate (FarLine or NearCoastline).
3) Enter the folder corresponding to either the homogeneous or heterogeneous case.
4) Run the file *CF_noise_sources_?.m* where *?* corresponds either to *Line* or *CircleN*.

**Outputs**
The simulation will produce a video showing noise propagation, a figure showing the variations of shear modulus, and a .mat file containig the simulated noise data at all stations considered.
