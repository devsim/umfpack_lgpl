
This is a fork of UMFPACK from SuiteSparse Version 3.0.0.

This was the last version of UMFPACK (5.1) and AMD under the LGPL version 2.1.

This is to allow the dynamic linking of this package following the terms of the license.  Please see the COPYING file for more details.

The ``umfpack/blaswrapper`` directory is in development and will be a dll shim making it possible to easily switch the BLAS/LAPACK libraries.

The build system for this project is ``CMAKE``.  The symbol visibility specified in the ``CMakeLists.txt`` files is hidden.  These symbols will be exposed as necessary to get a working implementation.

Since this is a fork from a large repository, it is most efficient to do a shallow clone:
```
git clone -depth 1 https://github.com/tcaduser/umfpack_lgpl.git
```
