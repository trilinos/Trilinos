Known issues
############

.. role:: cppkokkos(code)
    :language: cppkokkos

BLAS/LAPACK
===========

- OpenBLAS:
   - some versions of openblas prior to 0.3.26 are not implementing correctly the asum/norm1 functions for vectors of length less than 16. On Host for these vector lengths Kokkos Kernels uses an internal implementation as a work around. See issue `#2005 <https://github.com/kokkos/kokkos-kernels/issues/2005>`_ for more informations on that topic.
