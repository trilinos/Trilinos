API: LAPACK
###########

.. toctree::
   :maxdepth: 2
   :hidden:

   lapack/gesv
   lapack/gesvd
   lapack/trtri


LAPACK support
==============

Below are tables summarizing the currently supported function calls and third party libraries in Kokkos Kernels.

.. list-table::
   :widths: 12 26 10 10 10 10 10
   :header-rows: 1

   * - LAPACK Call
     - API Call
     - Reference
     - BLAS
     - cuBLAS
     - rocBLAS
     - oneMKL
   * - getrf
     - 
     - 
     - 
     - 
     - 
     - 
   * - getrs
     - 
     - 
     - 
     - 
     - 
     - 
   * - gesv
     - :doc:`gesv <lapack/gesv>`
     - --
     - X
     - X
     - X
     - --
   * - getri
     - 
     - 
     - 
     - 
     - 
     - 
   * - trtri
     - :doc:`trtri <lapack/trtri>`
     - --
     - X
     - X
     - --
     - --
   * - geequ
     - 
     - 
     - 
     - 
     - 
     - 
   * - geqrf
     - 
     - 
     - 
     - 
     - 
     - 
   * - ungqr
     - 
     - 
     - 
     - 
     - 
     - 
   * - unmqr
     - 
     - 
     - 
     - 
     - 
     - 
   * - gesvd
     - :doc:`gesvd <lapack/gesvd>`
     - --
     - X
     - X
     - X
     - --
   * - geev
     - 
     - 
     - 
     - 
     - 
     - 
