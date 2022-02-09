# libdogleg-f

Fortran bindings to [libdogleg](https://github.com/dkogan/libdogleg) - a large-scale nonlinear least-squares optimization library.

Currently only the dense optimizer calls are supported.


## Installation

To use `libdogleg-f` you will need:

* a Fortran compiler supporting Fortran 2018 (GCC 5 or newer or Intel Fortran)
* the [Fortran package manager (fpm)]()
* [libdogleg](https://github.com/dkogan/libdogleg)
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