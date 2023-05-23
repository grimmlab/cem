# Convex Envelope Method (CEM)

Official repository for the CEM (currently under review, preprint: https://arxiv.org/abs/2304.12025).

Description of the CEM in general:
The CEM constructs all liquid phase equilibria over the whole composition space for a given system with an arbitrary number of components. This is an important feature that can be used within chemical process simulations. To achieve this, the composition space is discretized, and the convex envelope of the Gibbs energy graph is computed. Employing the tangent plane criterion, all liquid phase equilibria can be determined robustly.

This repository contains the implementation used to produce the results published alongside the CEM. Required Python packages are listed at the end of this summary. 

To reproduce the results, download the code and the folder structure as provided and run main.py. If one sets actors_parallelized in main.py to an integer greater than 0, the code will be run parallelized using ray.

To use the CEM on other examples, provide property data (e.g., NRTL interactions parameters and feed streams) within property_data.py in the same way as for the published systems and start the process via main.py 


Requirements: os, time, numpy, scipy, itertools, copy, ray, matplotlib, pandas, importlib.
