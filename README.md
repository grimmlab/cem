# Convex Envelope Method (CEM) 

## Purpose of this repository / Quickstart
This repository contains official implementation of the CEM, which enables to reproduce all results shown in [1] (therein, the CEM is explained in detail). Download the repository (and the underlying folder structure) and run main.py to reproduce the results from [1]. On a local machine, we recommend setting all items in the dictionary `actors_parallelized` in `main.py` to 0 (this will disable the usage of parallelization).

## General description of the CEM
The CEM constructs all liquid phase equilibria over the whole composition space for a given system with an arbitrary number of components. This is an important feature that can be used within chemical process simulations as it enables the calculation of occurring phase splits in liquid mixtures. To achieve this, the composition space is discretized, and the convex envelope of the Gibbs energy graph is computed. Employing the tangent plane criterion, all liquid phase equilibria can be determined robustly. 

## Overview of repository
### How does this implementation work?
The CEM proceeds in the following order:
##### 1. Discretization of composition space
The composition space is discretized within `point_discretization.py` by specifying the number of components and the minimum distance between the points (specified by recursion steps $r$). This distance is equal to $\frac{1}{2^r}$ in `point_discretization.PointDisc`. An example of generating multiple discretizations is right at the start of `main.py`.

##### 2. Determination of Gibbs energy graph and convex envelope
Steps 2-4 are executed within `MiscibilityAnalysis` (see `lle.py`). The Gibbs energy graph is computed using $g^E$-models as NRTL or UNIQUAC (see `thermo_models.py`). The convex envelope is found by scipy.spatial.ConvexHull.

##### 3. Classification of simplices of the convex envelope
In this step, all simplices of the convex envelope are classified (either they model a phase split or not). This is done within `MiscibilityAnalysis` (see `lle.py`) and can be executed parallelized if specified (see below). After the classification, the phase equilibrium is stored (mainly within .npy files) and can be used to compute phase splits. For systems with 3 or 4 components, `plotter.py` allows plotting the resulting phase equilibria.

##### 4. Computation of phase splits
The function `MiscibilityAnalysis.find_phase_split` (see `lle.py`) allows the computation of a phase split for a given feed flowrate. Beforehand, the phase equilibrium has to be constructed (steps 1-3) or loaded.

### Usage of the CEM for other systems
To use the CEM on other examples, provide property data (e.g., NRTL interactions parameters and feed streams) within `property_data.py` like for the published systems. Using `property_data.InitClass` within main.py, the CEM can access the specified property data and construct the desired phase equilibria.

### Parallelization
As described in [1], parts of the CEM can be executed in parallel. For this purpose, we utilize [ray.io](https://www.ray.io/). If one wants to use the parallelized version of this implementation, specify the variable `actors_parallelized` in `property_data.InitClass` with some integer greater than 0 (can also be set in `main.py` in the dictionary `actors_parallelized`).

## Cite

To cite this work, use the following BibTex-entry:

```
@article{gottl2023convex,
  title={Convex Envelope Method for determining liquid multi-phase equilibria in systems with arbitrary number of components},
  author={G{\"o}ttl, Quirin and Pirnay, Jonathan and Grimm, Dominik G and Burger, Jakob},
  journal={Computers \& Chemical Engineering},
  pages={108321},
  year={2023},
  publisher={Elsevier}
}
```

## References
[1] Göttl, Q., Pirnay, J., Grimm, D. G., & Burger, J. (2023). Convex Envelope Method for determining liquid multi-phase equilibria in systems with arbitrary number of components. Computers & Chemical Engineering, 108321. [https://doi.org/10.1016/j.compchemeng.2023.108321](https://doi.org/10.1016/j.compchemeng.2023.108321)
