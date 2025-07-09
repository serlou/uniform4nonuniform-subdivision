# A uniform non-linear subdivision scheme reproducing polynomials at any non-uniform grid

Subdivision schemes are efficient tools to generate curves from a sample of points. This particular scheme was designed to deal with data describing smooth curves in R^n, n>=2. More information about this scheme can be found in [1].

This repository includes:
1. MATLAB implementation of a uniform subdivision scheme that reproduces parabolas on non-uniform grids. In an example, the curvature is computed using curvature.m and circumcenter.m from [2].
2. Mathematica notebook with symbolic computations complementing the proofs in [1].

Copyright (c) 2023 Sergio López-Ureña

[1] Sergio López-Ureña (2024). A uniform non-linear subdivision scheme reproducing polynomials at any non-uniform grid. Applied Mathematics and Computation, 479, 128889. https://doi.org/10.1016/j.amc.2024.128889

[2] Are Mjaavatten (2024). Curvature of a 1D curve in a 2D or 3D space (https://www.mathworks.com/matlabcentral/fileexchange/69452-curvature-of-a-1d-curve-in-a-2d-or-3d-space), MATLAB Central File Exchange. Retrieved January 11, 2024.

This research has been supported through projects CIAICO/2021/227, PID2020-117211GB-I00 funded by MCIN/AEI/10.13039/501100011033 and PID2023-146836NB-I00 funded by MCIU/AEI/10.13039/501100011033.
