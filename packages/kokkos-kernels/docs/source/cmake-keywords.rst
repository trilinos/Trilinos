CMake Options
#############

Test Options
============

- ``KokkosKernels_ENABLE_TESTS``: BOOL

  - Whether to build tests.
  - Default: OFF

- ``KokkosKernels_TEST_ETI_ONLY``: BOOL

  - Whether to restrict testing to ETI types.
  - Default: ON

- ``KokkosKernels_ENABLE_PERFTESTS``: STRING

  - Whether to build performance tests
  - Default: OFF

- ``KokkosKernels_ENABLE_TESTS_AND_PERFSUITE``: BOOL

  - Whether to build performance tests and suite.
  - Default: OFF

Examples Options
================

- ``KokkosKernels_ENABLE_EXAMPLES``: BOOL

  - Whether to build examples.
  - Default: OFF

Component Options
=================

- ``KokkosKernels_ENABLED_COMPONENTS``: BOOL

  - A list of components to enable in testing and building
  - Default: ALL
  - Valid: BATCHED BLAS LAPACK GRAPH SPARSE ALL

- ``KokkosKernels_ENABLE_ALL_COMPONENTS``: BOOL
- ``KokkosKernels_ENABLE_COMPONENT_BATCHED``: BOOL
- ``KokkosKernels_ENABLE_COMPONENT_BLAS``: BOOL
- ``KokkosKernels_ENABLE_COMPONENT_LAPACK``: BOOL
- ``KokkosKernels_ENABLE_COMPONENT_SPARSE``: BOOL
- ``KokkosKernels_ENABLE_COMPONENT_GRAPH``: BOOL
- ``KokkosKernels_ENABLE_COMPONENT_ODE``: BOOL

ETI Options
===========

- ``KokkosKernels_INST_HALF``: BOOL

  - Whether to pre instantiate kernels for the scalar type Kokkos::Experimental::half_t.  Disabling this may increase build times
  - Default: OFF

- ``KokkosKernels_INST_BHALF``: BOOL

  - Whether to pre instantiate kernels for the scalar type Kokkos::Experimental::bhalf_t.  Disabling this may increase build times
  - Default: OFF

- ``KokkosKernels_ETI_ONLY``: BOOL

  - Whether to restrict availability of kernels to ETI types only. Turning this on guarantees that kernels are never built inside of object files which simply call KokkosKernels functions.
  - Default: OFF

- ``KokkosKernels_INST_COMPLEX_DOUBLE``: BOOL

  - Whether to pre instantiate kernels for the scalar type complex<double>.  Disabling this may increase build times.
  - Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_DOUBLE.

- ``KokkosKernels_INST_COMPLEX_FLOAT``: BOOL

  - Whether to pre instantiate kernels for the scalar type complex<float>.  Disabling this may increase build times.
  - Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_FLOAT.

- ``KokkosKernels_INST_DOUBLE``: BOOL

  - Whether to pre instantiate kernels for the scalar type double.  This option is KokkosKernels_INST_DOUBLE=ON by default.  Disabling this may increase build times.
  - Default: ON

- ``KokkosKernels_INST_EXECSPACE_OPENMP``: BOOL

  - Whether to pre instantiate kernels for the execution space Kokkos::OpenMP.  Disabling this when Kokkos_ENABLE_OpenMP is enabled may increase build times.
  - Default: ON if Kokkos is OpenMP-enabled, OFF otherwise.

- ``KokkosKernels_INST_EXECSPACE_SERIAL``: BOOL

  - Whether to build kernels for the execution space Kokkos::Serial.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_SERIAL is enabled may increase build times.
  - Default: ON when Kokkos is Serial-enabled, OFF otherwise.

- ``KokkosKernels_INST_EXECSPACE_THREADS``: BOOL

  - Whether to build kernels for the execution space ``Kokkos::Threads``.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_PTHREAD is enabled may increase build times.
  - Default: ON if Kokkos is Threads-enabled, OFF otherwise.

- ``KokkosKernels_INST_FLOAT``: BOOL

  - Whether to pre instantiate kernels for the scalar type float.  Disabling this may increase build times.
  - Default: OFF or unless enabled during a Trilinos build with ``Trilinos_ENABLE_FLOAT``.

- ``KokkosKernels_INST_LAYOUTLEFT``: BOOL

  - Whether to pre instantiate kernels for the view layout LayoutLeft.  This option is ``KokkosKernels_INST_LAYOUTLEFT=ON`` by default.  Disabling this may increase build times.
  - Default: ON

- ``KokkosKernels_INST_LAYOUTRIGHT``: BOOL

  - Whether to pre instantiate kernels for the view layout LayoutRight.  This option is ``KokkosKernels_INST_LAYOUTRIGHT=OFF`` by default.  Disabling this may increase build times.
  - Default: OFF

- ``KokkosKernels_INST_MEMSPACE_HOSTSPACE``: BOOL

  - Whether to pre instantiate kernels for the memory space ``Kokkos::HostSpace``.  Disabling this when one of the Host execution spaces is enabled may increase build times.
  - Default: ON

- ``KokkosKernels_INST_OFFSET_INT``: BOOL

  - Whether to pre instantiate kernels for the offset type int. This option is KokkosKernels_INST_OFFSET_INT=ON by default.
  - Default: ON

- ``KokkosKernels_INST_OFFSET_SIZE_T``: BOOL

  - Whether to pre instantiate kernels for the offset type size_t.  This option is KokkosKernels_INST_OFFSET_SIZE_T=OFF by default.
  - Default: ON

- ``KokkosKernels_INST_ORDINAL_INT``: BOOL

  - Whether to pre instantiate kernels for the ordinal type int.  This option is KokkosKernels_INST_ORDINAL_INT=ON by default.
  - Default: ON

- ``KokkosKernels_INST_ORDINAL_INT64_T``: BOOL

  - Whether to pre instantiate kernels for the ordinal type int64_t.  This option is KokkosKernels_INST_ORDINAL_INT64_T=OFF by default.
  - Default: OFF

- ``KokkosKernels_ADD_DEFAULT_ETI``: BOOL

  - Whether to include a set of default ETI instantiations (otherwise only those explicitly requested will be included
  - Default: OFF

- ``KokkosKernels_INST_EXECSPACE_CUDA``: BOOL

  - Whether to pre instantiate kernels for the execution space ``Kokkos::Cuda``. Disabling this when ``Kokkos_ENABLE_CUDA`` is enabled may increase build times.
  - Default: ON if Kokkos is CUDA-enabled, OFF otherwise

- ``KokkosKernels_INST_MEMSPACE_CUDAUVMSPACE``: BOOL

  - Whether to pre instantiate kernels for the memory space ``Kokkos::CudaUVMSpace``.  Disabling this when ``Kokkos_ENABLE_CUDA`` is enabled may increase build times.
  - Default: OFF.

- ``KokkosKernels_INST_EXECSPACE_HIP``: BOOL
- ``KokkosKernels_INST_MEMSPACE_HIPSPACE``: BOOL
- ``KokkosKernels_INST_EXECSPACE_SYCL``: BOOL
- ``KokkosKernels_INST_MEMSPACE_SYCLSPACE``: BOOL
- ``KokkosKernels_INST_EXECSPACE_OPENMPTARGET``: BOOL
- ``KokkosKernels_INST_MEMSPACE_OPENMPTARGETSPACE``: BOOL
- ``KokkosKernels_INST_MEMSPACE_HBWSPACE``: BOOL

Documentation Options
=====================

- ``KokkosKernels_ENABLE_DOCS``: BOOL

  - Whether to build documentation
  - Default: OFF

Supported Third-party Libraries and Related CMake Options
=========================================================

- ``KokkosKernels_NO_DEFAULT_CUDA_TPLS``: BOOL

  - Whether CUDA TPLs should be enabled by default.
  - Default: OFF
- ``KokkosKernels_NO_DEFAULT_ROCM_TPLS``: BOOL

- ARMPL

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_ARMPL``: BOOL

      - Whether to enable ARMPL
      - Default: OFF

    - ``ARMPL_LIBRARIES``: STRING
      
      - the name of the armpl library files to look for (possibly combined with above)

    - ``ARMPL_LIBRARY_DIRS``: STRING
      
      - which directories to look for ``armpl``/``armpl_mp`` (optionally combined with below)

  - Environment Variables

    - *under construction...*
- BLAS

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_BLAS``: BOOL

      - Whether to enable BLAS
      - Default: OFF

    - ``BLAS_LIBRARIES``: STRING

      - Optional override for the libraries that comprise TPL BLAS.
      - Default: None. Default common library names will be searched

    - ``BLAS_LIBRARY_DIRS``: STRING

      - Optional override for the library directories that comprise TPL BLAS.
      - Default: None. Default common library locations will be searched

  - Environment Variables

    - *under construction...*

- CBLAS

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_CBLAS``: BOOL

      - Whether to enable CBLAS
      - Default: OFF

    - ``CBLAS_ROOT``
    - ``CBLAS_LIBRARIES``
    - ``CBLAS_LIBRARY_DIRS``
    - ``CBLAS_INCLUDE_DIRS``

  - Environment Variables

    - *under construction...*

- CHOLMOD

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_CHOLMOD``: BOOL

      - Whether to enable CHOLMOD
      - Default: OFF

    - ``CHOLMOD_ROOT``
    - ``CHOLMOD_LIBRARIES``
    - ``CHOLMOD_LIBRARY_DIRS``
    - ``CHOLMOD_INCLUDE_DIRS``

  - Environment Variables

    - *under construction...*

- CUBLAS

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_CUBLAS``: BOOL

      - Whether to enable CUBLAS
      - Default: ON if CUDA-enabled Kokkos, otherwise OFF

    - ``TPLCUBLAS_ROOT``

      - Optional override for the library directories that comprise TPL CUSOLVER. As a side-effect of ``find_library`` being called from within ``find_package(CUDA)`` which we use to find the CUDA libraries

    - ``CUBLAS_ROOT``: Optional override for the root directory of CUBLAS install (contains e.g. ``include``, ``lib64``)
    - ``KokkosKernels_CUBLAS_ROOT``: Optional override for the root directory of CUBLAS install (contains e.g. ``include``, ``lib64``)
    - ``CUBLAS_LIBRARIES``: paths to CUBLAS libraries
    - ``CUBLAS_INCLUDE_DIRS``: paths to CUBLAS include dirs
    - ``CUBLAS_LIBRARY_DIRS``: paths directories containing CUBLAS libraries

  - Environment Variables
    - *under construction...*

- CUSOLVER

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_CUSOLVER``: BOOL

      - Whether to enable CUSOLVER
      - Default: ON if CUDA-enabled Kokkos, otherwise OFF

    - ``TPLCUSOLVER_ROOT``

      - Optional override for the library directories that comprise TPL CUSOLVER. (As a side-effect of ``find_library`` being called from within ``find_package(CUDA)`` which we use to find the CUDA libraries).
      - Default: None

    - ``CUSOLVER_ROOT``: Optional override for the root directory of CUSOLVER install (contains e.g. ``include``, ``lib64``)
    - ``KokkosKernels_CUSOLVER_ROOT``: Optional override for the root directory of CUSOLVER install (contains e.g. ``include``, ``lib64``)
    - ``CUSOLVER_LIBRARIES``: paths to CUSOLVER libraries
    - ``CUSOLVER_INCLUDE_DIRS``: paths to CUSOLVER include dirs
    - ``CUSOLVER_LIBRARY_DIRS``: paths directories containing CUSOLVER libraries

  - Environment Variables

    - *under construction...*

- CUSPARSE

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_CUSPARSE``: BOOL

      - Whether to enable CUSPARSE
      - Default: ON if CUDA-enabled Kokkos, otherwise OFF

    - ``TPLCUSPARSE_ROOT``: STRING

      - Optional override for the library directories that comprise TPL CUSPARSE. (As a side-effect of ``find_library`` being called from within ``find_package(CUDA)`` which we use to find the CUDA libraries).
      - Default: None

    - ``CUSPARSE_ROOT``: Optional override for the root directory of CUSPARSE install (contains e.g. ``include``, ``lib64``)
    - ``KokkosKernels_CUSOLVER_ROOT``: Optional override for the root directory of CUSPARSE install (contains e.g. ``include``, ``lib64``)
    - ``CUSPARSE_LIBRARIES``: paths to CUSPARSE libraries
    - ``CUSPARSE_INCLUDE_DIRS``: paths to CUSPARSE include dirs
    - ``CUSPARSE_LIBRARY_DIRS``: paths directories containing CUSPARSE libraries

  - Environment Variables

    - *under construction...*

- MAGMA

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_MAGMA``: BOOL

      - Whether to enable MAGMA
      - Default: OFF

    - ``KokkosKernels_MAGMA_ROOT``: PATH

      - Location of MAGMA install root.
      - Default: None or the value of the environment variable MAGMA_ROOT if set

    - ``MAGMA_LIBRARY_DIRS``
    - ``MAGMA_LIBRARIES``

  - Environment Variables

    - ``MAGMA_DIR``

- MKL

  - supported versions:
  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_MKL``: BOOL

      - Whether to enable MKL
      - Default: OFF

    - ``MKL_LIBRARIES``: STRING

      - Optional override for the libraries that comprise TPL MKL.
      - Default: None. Default common library names will be searched

    - ``MKL_LIBRARY_DIRS``: STRING

      - Optional override for the library directories that comprise TPL MKL.
      - Default: None. Default common library locations will be searched

    - ``KokkosKernels_MKL_ROOT``: PATH

      - Location of MKL install root.
      - Default: None or the value of the environment variable MKL_ROOT if set

  - Environment Variables

    - ``MKLROOT``

- LAPACK

  - CMake Options

    - ``KokkosKernels_ENABLE_TPL_LAPACK``: BOOL

      - Whether to enable LAPACK
      - Default: ON if BLAS is enabled, otherwise OFF

    - ``LAPACK_LIBRARIES``: STRING

      - Optional override for the libraries that comprise TPL LAPACK.
      - Default: None. Default common library names will be searched

    - ``LAPACK_LIBRARY_DIRS``: STRING

      - Optional override for the library directories that comprise TPL LAPACK.
      - Default: None. Default common library locations will be searched

    - ``KokkosKernels_LAPACK_ROOT``: PATH

      - Location of LAPACK install root.
      - Default: None or the value of the environment variable LAPACK_ROOT if set

  - Environment Variables

    - *under construction...*

- LAPACKE

  - CMake Options

    - ``LAPACKE_LIBRARIES``
    - ``LAPACKE_LIBRARY_DIRS``
    - ``LAPACKE_INCLUDE_DIRS``

  - Environment Variables

    - ``LAPACKE_ROOT``
    - ``KokkosKernels_LAPACKE_ROOT``
    - ``OPENBLAS_ROOT``

- METIS

  - CMake Options

    - ``METIS_LIBRARY_DIRS``
    - ``METIS_INCLUDE_DIRS``

  - Environment Variables

    - ``METIS_ROOT``
    - ``KokkosKernels_METIS_ROOT``

- ROCBLAS

  - deferred to ``find_package(ROCBLAS)``

- ROCSOLVER

  - deferred to ``find_package(ROCSOLVER)``

- ROCSPARSE

  - deferred to ``find_package(ROCSPARSE)``

- SUPERLU

  - CMake Options

    - ``SUPERLU_LIBRARIES``
    - ``SUPERLU_LIBRARY_DIRS``
    - ``SUPERLU_INCLUDE_DIRS``

  - Environment Variables

    - ``SUPERLU_ROOT``
    - ``KokkosKernels_SUPERLU_ROOT``

Other Options
=============

- KokkosKernels_ENABLE_EXPERIMENTAL: BOOL

  - Enable building and installation of experimental KokkosKernels features.
  - Default: OFF

- KokkosKernels_LINALG_OPT_LEVEL: BOOL

  - Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems.
  - Default: 1

- ``KokkosKernels_ENABLE_SUPERNODAL_SPTRSV``: BOOL

  - Whether to build supernodal SPTRSV support
  - Default: ON

Generated using ``grep -r KOKKOSKERNELS_ADD_TPL_OPTION``
