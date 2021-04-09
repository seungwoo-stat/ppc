# Probabilistic Principal Curves on Riemannian Manifolds

This repo provides functions implementing Probabilistic Principal Curves on Riemannian Manifolds (abbreviated as PPC) using statistical software R.

 ## Overview

There are four files implementing PPC on three spaces.

- `generics.R` : defines generic functions. e.g. `PPC` `DIST` `LOG` etc.

- `ppc_euc.R` : building block for implementing PPC in 2D and 3D Euclidean space.

- `half-plane/half_plane_funs.R` : PPC in the Poincare half-plane (2D hyperbolic space)
- `sphere-2-3/sphere_funs.R` : PPC in the 2D and 3D sphere
- `so3-group/so3_funs.R` : PPC in the 3D special orthogonal group SO(3)

Note that description of each function is provided as comment in each file. (Some are missing)



Additionaly, following files are for simulation & comparison with other methods.

- `half-plane/half_plane_simulations.R` 
  - simulates 4 datasets on the half-plane
  - compares PPC with existing methods visually
- `sphere-2-3/sphere2_simulations.R` 
  - simulates 2 datasets on the 2D sphere
  - compares PPC with existing methods visually
- `sphere-2-3/compare.R` 
  - compare `wave` data on 2D sphere on several methods numerically
- `sphere-2-3/simplex3_funs.R` 
  - defines `simplex3` class and `plot` method for the class
- `so3-group/so3_simulations.R`
  - simulates 2 datasets on the 3D special orthogonal group



We also provide code that was used to extract election data from the excel files provided by [National Election Commission of South Korea](http://info.nec.go.kr). Excel files are in [`./sphere-2-3/data`](./sphere-2-3/data).

- `sphere-2-3/real_data.R` 
  - extracts data from excel files in `sphere-2-3/data` and save as `sphere-2-3/19pe.csv`
  - fit PPC to the data, and plot
  - conducts dimension reduction using the fitted PPC
- `sphere-2-3/19pe.csv` 
  - the real data used for analysis; can also be generated using `sphere-2-3/real_data.R`



## Reference

- S. Kang, and H.-S. Oh. "Probabilistic Principal Curves on Riemannian Manifolds". arXiv preprint. Apr. 2021. (Link to be updated)