# Discrete_sheet_impedance
---------------------------------------------------------------------------
This repository contains MATLAB code to reproduce some of the results from the paper [Full Manipulation of the Power Intensity Pattern in a Large Space-Time Digital Metasurface: From Arbitrary Multibeam Generation to Harmonic Beam Steering Scheme](<https://onlinelibrary.wiley.com/doi/full/10.1002/andp.202000321>). 

# Abstract
---------------------------------------------------------------------------
Beyond the scope of space-coding metasurfaces, space-time digital metasurfaces can substantially expand the application scope of digital metamaterials. Fully manipulate the power intensity pattern of the metasurfaces can give us fabulous flexibility and privilege and is highly demanded in divers practical applications, such as direct broadcasting and multiple-target radar systems. By the MATLAB code provided in this repository, one can generates the scattered pattern of the space-time digital MS as well as the required local phase and amplitudes of the surface. It generates phase matrices, discretizes them, computes the array factors across frequencies, and plots the resulting directivity and surface distribution patterns.

# Results
---------------------------------------------------------------------------
In the file spacetime.m, you can define your frequency, reflection angles, power associated to each channels,.... 
In the below [figure](https://github.com/Javadio/Space_time_MS/blob/main/ig3.PNG), you can find two radiation beams with their corresponding phase and amplitude distribution. 
<h2>Figure: Z Values Visualization</h2>

<p align="center">
  <img src="https://github.com/Javadio/Space_time_MS/blob/main/ig3.PNG" alt="Z Values" width="900">
</p>

<p align="center"><b>Figure 1:</b> Directivity intensity pattern (in both linear and decibel formats) for asymmetrically oriented two-beams time-modulated metasurface pointing at $\theta = 15^{o}$ and $\theta = 65^{o}$.</p>
