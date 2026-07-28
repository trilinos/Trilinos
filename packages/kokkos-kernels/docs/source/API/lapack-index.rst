API: LAPACK
###########

.. toctree::
   :maxdepth: 2
   :hidden:

   lapack/geqrf
   lapack/gemqr
   lapack/gegqr
   lapack/potrf
   lapack/potrs
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
     - :doc:`geqrf <lapack/geqrf>`
     - --
     - X
     - X
     - X
     - --
   * - potrf
     - :doc:`potrf <lapack/potrf>`
     - --
     - X
     - X
     - X
     - --
   * - potrs
     - :doc:`potrs <lapack/potrs>`
     - --
     - X
     - X
     - X
     - --
   * - {or,un}gqr
     - :doc:`gegqr <lapack/gegqr>`
     - --
     - X
     - X
     - X
     - --
   * - {or,un}mqr
     - :doc:`gemqr <lapack/gemqr>`
     - --
     - X
     - X
     - X
     - --
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
