# Simulations of dipole trapping and Rabi cycling of cold CO molecules

The code contained in this repository refers to some of the simulations I performed during my Ph.D. at LENS (University of Florence, Italy) working in the G. Santambrogio's cold molecules lab. Some codes are written in C/C++ and integrated with bash scripting in order to make it easy to use within Linux environment provided by GnuPlot, while others are written within Wolfram Mathematica 12.

This repo is structured as follows:

### Mathematica
It contains files.nb covering different tasks related to both optical dipole trapping and Rabi oscillations.

### Optical Trapping Animation
It contains the files needed to generate 3D_Animation.gif, being it a simple 3D representation of molecular trajectories in a force field, depending on time-varying molecular quantum state.

### Rabi Oscillations
It contains files to simulate the Rabi cycling for an ensemble of CO molecules freely expanding in the vacuum and crossing a laser field.

### Spontaneous Emission
It contains files to "graphically check" the correct implementation of the decay process, from the upper to the lower state in the modeled two-levels system.

### Trapping Probability
It contains files to simulate the molecular trajectories and to predict the number of finally trapped ground-state CO molecules by varying different subset of experimental parameters.
