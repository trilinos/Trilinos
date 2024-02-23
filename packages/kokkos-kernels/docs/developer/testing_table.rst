Testing Table
=============

SAND2023-05267O [#]_

Below is a testing table summarizing the KokkosKernels continuous integration and nightly test coverage.

The following is a description of abbreviations used throughout the testing table.

* ETI: Explicit template instantiation
* PR: Pull Request
* LEFT: LayoutLeft
* RIGHT: LayoutRight
* REL: CMake release build type
* DBG: CMake debug build type
* BCHK: Kokkos core bounds checking
* NOETI: No default ETI types included
* UVM: Unified Memory (Cuda)

The following is a description of column headings in the testing table.

* Project: the jenkins project name for the test case
* Architectures: the test case's coverage architectures
* Compilers: the covered compilers
* Backends: the covered kokkos core backends
* Scalars: the covered ETI'd scalar types
* Ordinals: the covered ETI'd ordinal types
* Offsets: the covered ETI'd offset types
* Layouts: the covered ETI'd kokkos core layout types

.. list-table::
    :align: center
    :header-rows: 1
    :stub-columns: 0
    :width: 100%
    :widths: auto


    * - Project
      - Architectures
      - Compilers
      - Backends
      - Scalars
      - Ordinals
      - Offsets
      - Layouts

    * * `PR_A64FX_ARMPL2110_OPENMP_LEFT_OPENBLAS_OPENLAPACK_REL`
      * A64FX
      * ARMPL 21.1.10
      * OpenMP
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_A64FX_ARMPL2110_OPENMP_LEFT_OPENBLAS_OPENLAPACK_REL`
      * A64FX
      * ARMPL 21.1.10
      * OpenMP
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_A64FX_GCC1020_OPENMP_SERIAL_LEFT_REL`
      * A64FX
      * GNU 10.2.0
      * OpenMP,Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_VEGA908_ROCM520_HIP_SERIAL_LEFT_REL`
      * VEGA908
      * ROCM 5.2.0
      * Hip, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PRTONIGHTLY_VEGA908_ROCM520_HIP_SERIAL_LEFT_OPENBLAS_OPENLAPACK_REL`
      * VEGA908
      * ROCM 5.2.0
      * Hip, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_POWER9_VOLTA70_GCC930_CLANG13_CUDA10_OPENMP_SERIAL_CUDA_LEFT_OPENBLAS_OPENLAPACK_REL`
      * Power8, Pascal60 -- Power9, Volta70
      * GNU 9.3.0 -- Clang 13.0.0, Cuda 10.1.243
      * OpenMp, Serial -- Cuda
      * double, `complex_double`
      * int
      * int, size_t
      * LayoutLeft

    * * `PR_POWER9_VOLTA70_CUDA11_OPENMP_CUDA_LEFT_RIGHT_REL`
      * Power9, Volta70
      * GNU 8.3.1, Cuda 11.2.2
      * Cuda, OpenMP
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft, LayoutRight

    * * `PR_SKX_GNU1020_OPENMP_LEFT_REL_NOETI`
      * Skx
      * GNU 10.2.0
      * OpenMP
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_SKX_GNU1020_THREADS_SERIAL_RIGHT_REL`
      * Skx
      * GNU 10.2.0
      * Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutRight

    * * `PR_SKX_GNU1020_OPENMP_SERIAL_LEFT_OPENBLAS_OPENLAPACK_REL`
      * Skx
      * GNU 10.2.0
      * Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_SKX_INTEL19_OPENMP_LEFT_MKLBLAS_MKLLAPACK_REL`
      * Skx
      * Intel 19.5.281
      * OpenMP
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_SKX_CLANG1001_THREADS_SERIAL_LEFT_REL`
      * Skx
      * Clang 10.0.1
      * Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `PR_NONE_CLANG14001_SERIAL_LEFT_RIGHT_RELWITHDBG_BCHK`
      * NONE
      * Clang 14.0.0.14000029
      * Serial
      * double, float, `complex_double`, `complex_float`
      * int
      * int, `size_t`
      * LayoutLeft, LayoutRight

    * * `PR_NONE_CLANG14001_THREADS_LEFT_RIGHT_RELWITHDBG_BCHK`
      * NONE
      * Clang 14.0.0.14000029
      * Serial
      * double, float, `complex_double`, `complex_float`
      * int
      * int, `size_t`
      * LayoutLeft, LayoutRight

    * * `PR_NONE_CLANG14001_SERIAL_LEFT_RIGHT_DBG`
      * NONE
      * Clang 14.0.0.14000029
      * Serial
      * double, float, `complex_double`, `complex_float`
      * int
      * int, `size_t`
      * LayoutLeft, LayoutRight

    * * `PR_NONE_CLANG14001_SERIAL_LEFT_RIGHT_REL_BCHK`
      * NONE
      * Clang 14.0.0.14000029
      * Serial
      * double, float, `complex_double`, `complex_float`
      * int
      * int, `size_t`
      * LayoutLeft, LayoutRight

    * * `NIGHTLY_SKX_GNU1020_OPENMP_THREADS_SERIAL_LEFT_DBG`
      * SKX
      * GNU 10.2.0
      * OpenMp, Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_GNU820_OPENMP_THREADS_SERIAL_LEFT_DBG`
      * SKX
      * GNU 8.2.0
      * OpenMp, Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_GNU820_OPENMP_THREADS_SERIAL_LEFT_REL`
      * SKX
      * GNU 8.2.0
      * OpenMp, Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_GNU920_OPENMP_THREADS_SERIAL_LEFT_DBG`
      * SKX
      * GNU 9.2.0
      * OpenMp, Threads, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_INTEL19_OPENMP_LEFT_DBG`
      * SKX
      * Intel 19.0.5
      * OpenMp
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_INTEL19_SERIAL_LEFT_DBG`
      * SKX
      * Intel 19.0.5
      * Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_INTEL19_THREADS_LEFT_DBG`
      * SKX
      * Intel 19.0.5
      * Threads
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_SKX_INTEL19_OPENMP_LEFT_MKL_DBG`
      * SKX
      * Intel 19.0.5
      * OPENMP
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_POWER9_VOLTA70_CUDA11_OPENMP_CUDA_LEFT_REL`
      * SKX
      * Cuda 11.2.2
      * OpenMP, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_POWER9_VOLTA70_CUDA11_SERIAL_CUDA_LEFT_REL`
      * SKX
      * Cuda 11.2.2
      * Serial, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_POWER9_VOLTA70_CUDA11_SERIAL_CUDA_LEFT_REL_UVM_RDC`
      * SKX
      * Cuda 11.2.2
      * Serial, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_POWER9_VOLTA70_CUDA11_SERIAL_CUDA_LEFT_DBG_BCHK`
      * SKX
      * Cuda 11.2.2
      * Serial, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_POWER9_VOLTA70_CUDA11_SERIAL_CUDA_LEFT_CUBLAS_CUSPARSE_REL_BCHK`
      * SKX
      * Cuda 11.2.2
      * Serial, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VEGA908_ROCM520_SERIAL_HIP_LEFT_REL`
      * VEGA908
      * Rocm 5.2.0
      * Serial, Hip
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VEGA908_ROCM520_SERIAL_HIP_LEFT_ROCBLAS_ROCSPARSE_REL`
      * VEGA908
      * Rocm 5.2.0
      * Serial, Hip
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VEGA906_ROCM520_SERIAL_HIP_LEFT_REL`
      * VEGA906
      * Rocm 5.2.0
      * Serial, Hip
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VEGA906_ROCM520_SERIAL_HIP_LEFT_DBG_BCHK`
      * VEGA906
      * Rocm 5.2.0
      * Serial, Hip
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_AMPERE80_CUDA11_SERIAL_CUDA_LEFT_DBG`
      * AMPHERE80
      * Cuda 11.7.99
      * Serial, Cuda
      * double
      * int
      * `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CLANG10_SERIAL_OPENMP_THREADS_LEFT_REL`
      * Volta70
      * Clang 10.0.0
      * Serial, OpenMP, Threads
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA10_CUDA_SERIAL_LEFT_RELWITHDBG`
      * Volta70
      * Cuda 10.1
      * Serial, Cuda
      * double
      * int
      * `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA117_CUDA_SERIAL_LEFT_RELWITHDBG`
      * Volta70
      * Cuda 11.7
      * Serial, Cuda
      * double
      * int
      * `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CLANG900_SERIAL_THREADS_LEFT_REL`
      * Volta70
      * Clang 9.0.0
      * Serial, Threads, `Threads_Serial`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CLANG900_SERIAL_THREADS_LEFT_DBG`
      * Volta70
      * Clang 9.0.0
      * Serial, Threads, `Threads_Serial`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CLANG900_SERIAL_THREADS_LEFT_REL_CPP20`
      * Volta70
      * Clang 9.0.0
      * Serial, Threads, `Threads_Serial`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA110_CUDA_OPENMP_LEFT_REL`
      * Volta70
      * Cuda 11.0
      * OpenMP, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA120_CUDA_OPENMP_LEFT_REL`
      * Volta70
      * Cuda 12.0
      * OpenMP, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA120_CUDA_OPENMP_LEFT_REL`
      * Volta70
      * Cuda 12.0
      * OpenMP, Cuda
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_GNU830_SERIAL_OPENMP_THREADS_LEFT_REL`
      * Volta70
      * Gnu 8.3.0
      * OpenMP, `OpenMP_Serial`, Serial, Threads, `Threads_Serial`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_GNU910_GNU920_SERIAL_OPENMP_THREADS_LEFT_REL`
      * Volta70
      * Gnu 9.1.0, Gnu 9.2.0
      * OpenMP, `OpenMP_Serial`, Serial, Threads, `Threads_Serial`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_GNU830_GNU910_SERIAL_OPENMP_LEFT_OPENBLAS_OPENLAPACK_REL`
      * Volta70
      * Gnu 9.1.0, Gnu 9.2.0
      * OpenMP, `OpenMP_Serial`, Serial, Threads, `Threads_Serial`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_A64FX_ARMPL2030_SERIAL_OPENMP_LEFT_ARMPLLBLAS_ARMPLSLAPACK_REL`
      * A64FX
      * Armpl 20.3.0
      * OpenMP, Serial
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA11_CUDA_OPENMP_SERIAL_PTHREAD_LEFT_REL`
      * Volta70
      * Cuda 11.1.0
      * `Cuda_OpenMP`, `Cuda_Serial`, `Cuda_Pthread`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA11_CUDA_OPENMP`, `SERIAL_PTHREAD_LEFT_DBG_BCHK`
      * Volta70
      * Cuda 11.1.0
      * `Cuda_OpenMP`, `Cuda_Serial`, `Cuda_Pthread`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA11_CUDA_OPENMP`, `SERIAL_PTHREAD_LEFT_DBG_BCHK`
      * Volta70
      * Cuda 11.1.0
      * `Cuda_OpenMP`, `Cuda_Serial`, `Cuda_Pthread`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_VOLTA70_CUDA11_CUDA_OPENMP`, `SERIAL_PTHREAD_LEFT_REL_UVM`
      * Volta70
      * Cuda 11.1.0
      * `Cuda_OpenMP`, `Cuda_Serial`, `Cuda_Pthread`
      * double, `complex_double`
      * int
      * int, `size_t`
      * LayoutLeft

    * * `NIGHTLY_HSW_INTEL19_OPENMP_LEFT_RELWITHDBG`
      * Hsw
      * Intel 19.1.3.20200925
      * OpenMP
      * double
      * int
      * `size_t`
      * LayoutLeft

    * * `NIGHTLY_KNL_INTEL19_OPENMP_LEFT_RELWITHDBG`
      * Hsw
      * Intel 19.1.3.20200925
      * OpenMP
      * double
      * int
      * `size_t`
      * LayoutLeft

.. rubric:: Footnotes

.. [#] This article has been authored by an employee of National Technology & Engineering Solutions of Sandia, LLC under Contract No. DE-NA0003525 with the U.S. Department of Energy (DOE). The employee owns all right, title and interest in and to the article and is solely responsible for its contents. The United States Government retains and the publisher, by accepting the article for publication, acknowledges that the United States Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this article or allow others to do so, for United States Government purposes. The DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan https://www.energy.gov/downloads/doe-public-access-plan. SAND2023-05267O.