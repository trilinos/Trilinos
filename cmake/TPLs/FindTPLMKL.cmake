
#-----------------------------------------------------------------------------
#  Intel's Math Kernel Library (MKL)
#
#  Acquisition information:
#    Date checked:  06 Aug 2012
#    Checked by:    Mark Hoemmen <mhoemme AT sandia.gov>
#    Version:       10.3
#

# Intel's Math Kernel Library (MKL) provides an implementation of the
# BLAS and LAPACK, which can be used to satisfy the BLAS and LAPACK
# TPLs.  However, the "MKL TPL" here refers to other functionality
# provided by the MKL, including things like sparse matrix kernels and
# pseudorandom number generators.  That's why we require a header
# file, to access the function declarations.

IF(Kokkos_ENABLE_SYCL)
  # For OneAPI MKL on GPU, use the CMake target
  # Temporarily change CMAKE_CXX_COMPILER to icpx to convince MKL to add the DPCPP target
  # If it sees that CMAKE_CXX_COMPILER is just ".../mpicxx", it won't do this.
  set(CMAKE_CXX_COMPILER_PREVIOUS "${CMAKE_CXX_COMPILER}")
  set(CMAKE_CXX_COMPILER "icpx")
  # Use the BLAS95 and LAPACK95 interfaces (int32_t for dimensions and indices)
  set(MKL_INTERFACE lp64)
  find_package(MKL REQUIRED COMPONENTS MKL::MKL MKL::MKL_SYCL)
  IF (NOT MKL_FOUND)
    MESSAGE(FATAL_ERROR "MKL (as CMake package) was not found! This is required for SYCL+MKL")
  ENDIF()
  set(CMAKE_CXX_COMPILER "${CMAKE_CXX_COMPILER_PREVIOUS}")

  tribits_extpkg_create_imported_all_libs_target_and_config_file( MKL
    INNER_FIND_PACKAGE_NAME MKL
    IMPORTED_TARGETS_FOR_ALL_LIBS MKL::MKL MKL::MKL_SYCL
    )
ELSE ()
  # For host MKL, the single library libmkl_rt is sufficient.
  # This works for older versions of MKL that don't provide MKLConfig.cmake.
  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( MKL
    REQUIRED_HEADERS mkl.h
    REQUIRED_LIBS_NAMES mkl_rt
    )
ENDIF()

# In the past, MKL users had to link with a long list of libraries.
# The choice of libraries enables specific functionality.  Intel
# provides a web page (the "Intel MKL Link Line Advisor") to help
# users pick which libraries and compiler flags to use:
#
# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
#
# As of version 10.3, MKL has an option to use a single dynamic
# library via "-lmkl_rt".  See the following article:
#
# http://software.intel.com/en-us/articles/a-new-linking-model-single-dynamic-library-mkl_rt-since-intel-mkl-103/
#
# This is why we made mkl_rt the required library name.  Users must
# override this if they are using an older MKL version, or if they
# prefer the long list of libraries to the single dynamic library
# option.  (I highly recommend the single dynamic library option, if
# your system allows dynamic libraries.)
#
# Users will probably need to specify the CMake options
# MKL_LIBRARY_DIRS (the path to the libraries) and MKL_INCLUDE_DIRS
# (the path to the header files).  On Linux, you may also have to link
# with Pthreads and libm (the C math library).  This may happen by
# default, but if you have trouble linking, try setting
# MKL_LIBRARY_NAMES to "mkl_rt;pthread;m" (for the single dynamic
# library option).  If you still have trouble, or if you are unable to
# use the single dynamic library option, look at examples of linking
# with MKL as the BLAS and LAPACK implementation in the sampleScripts/
# subdirectory of the Trilinos source tree.  Copy BLAS_LIBRARY_NAMES
# to MKL_LIBRARY_NAMES and try again.
#
# Users may also need to use special compiler flags, in particular if
# they wish to enable ILP64 support (for 64-bit integer indices).  See
# the above "Intel MKL Link Line Advisor" link for details.
#
# Here is a combination of CMake options that worked for me, when
# building and running on Linux, using the single dynamic library
# option with MKL 10.3:
#
# -D BLAS_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" 
# -D BLAS_LIBRARY_NAMES:STRING="mkl_rt" 
# -D LAPACK_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" 
# -D LAPACK_LIBRARY_NAMES:STRING="mkl_rt" 
# -D MKL_LIBRARY_DIRS:FILEPATH="${MKLROOT}/lib/intel64" 
# -D MKL_LIBRARY_NAMES:STRING="mkl_rt" 
# -D MKL_INCLUDE_DIRS:FILEPATH="${MKLROOT}/include" 
# 
# where the MKLROOT environment variable points to my MKL install
# directory.

