# libdogleg-f

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

Fortran bindings to [`libdogleg`](https://github.com/dkogan/libdogleg) - a large-scale nonlinear least-squares optimization library.

The main task of the library is to solve the unconstrained [non-linear least squares problem](http://en.wikipedia.org/wiki/Non-linear_least_squares), i.e. find the vector $\boldsymbol{p}$ that satisfies

$$\underset{\boldsymbol{p}}{\mathrm{min}}\quad{\frac{1}{2} \sum_i {\lVert f_i(\boldsymbol{p}) \rVert}^2}$$

The cost function $f: \boldsymbol{p} \rightarrow \boldsymbol{x}$ represents a mapping from the state variables $\boldsymbol{p} \in \mathbb{R}^n$ to a vector of residuals $\boldsymbol{x} \in \mathbb{R}^m$, $m > n$, of higher dimensionality than $\boldsymbol{p}$.

The user of the library passes a callback function that satisfies the `dl_callback_?` interface. The callback receives the vector $\boldsymbol{p}$ and returns the residual vector $\boldsymbol{x}$ and a $m \times n$ matrix of derivatives $\boldsymbol{J}_{ij} = \partial f_i / \partial p_j$, also known as the Jacobian matrix.

Depending on the size and sparsity of the problem, one can select between `dense` and `sparse` callback interfaces. For small and dense problems, the `dense` interface is recommended where the matrix $\boldsymbol{J}$ is stored as a rank-2 Fortran array. For large and sparse problems, the `sparse` interface uses a compressed-sparse column (CSC) representation of the Jacobian matrix, leading to large performance gains and low memory usage.

Currently only the `dense` optimizer calls are supported.

## Installation

To use `libdogleg-f` you will need:

* a Fortran compiler supporting the Fortran 2018 standard
* the [Fortran package manager (fpm)](https://github.com/fortran-lang/fpm)
* the [`libdogleg`](https://github.com/dkogan/libdogleg) library
* CHOLMOD and [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)
* LAPACK

Include `libdogleg-f` in your fpm package manifest using

```toml
[dependencies]
libdogleg-f.git = "https://github.com/ivan-pi/libdogleg-f"
```

## License and copyright

Copyright (C) 2022 Ivan Pribec

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program (see [COPYING.LESSER](./COPYING.LESSER) and [COPYING](./COPYING)). If not, see <https://www.gnu.org/licenses/>. 