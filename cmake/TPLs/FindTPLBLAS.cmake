if (TPL_ENABLE_MKL)
  # MKL can provide the BLAS library.
  # In Trilinos, FindTPLMKL is processed before FindTPLBLAS (this file) so
  # it should have already created the target BLAS::all_libs.
  if (NOT TARGET MKL::all_libs)
    message (FATAL_ERROR "\
MKL and BLAS are both enabled as TPLs, so MKL's libraries\
should also be used as the BLAS libraries. However, FindTPLMKL.cmake\
did not create the target MKL::all_libs. Please report this bug to\
the Trilinos developers.")
  endif ()

  if (Kokkos_ENABLE_SYCL)
    # FindTPLMKL already called find_package(MKL).
    # The MKL::MKL target was required for MKL when SYCL is enabled.
    tribits_extpkg_create_imported_all_libs_target_and_config_file( BLAS
      INNER_FIND_PACKAGE_NAME MKL
      IMPORTED_TARGETS_FOR_ALL_LIBS  MKL::MKL)
  else ()
    # Use the Tribits-generated target (tribits::MKL::mkl_rt) as the BLAS library
    tribits_extpkg_create_imported_all_libs_target_and_config_file( BLAS
      INNER_FIND_PACKAGE_NAME  MKL
      IMPORTED_TARGETS_FOR_ALL_LIBS tribits::MKL::mkl_rt)
  endif ()

  set (MKL_PROVIDES_BLAS_LAPACK ON CACHE BOOL "")
  message (STATUS "MKL will provide BLAS functions.")
else ()
  if (MSVC AND NOT
      (BLAS_LIBRARY_DIRS  OR
       (NOT "${BLAS_LIBRARY_NAMES}" STREQUAL "blas blas_win32" AND
        NOT "${BLAS_LIBRARY_NAMES}" STREQUAL "") OR
       BLAS_INCLUDE_DIRS  OR
       BLAS_INCLUDE_NAMES OR
       (NOT "${TPL_BLAS_LIBRARIES}" STREQUAL "blas" AND
        NOT "${TPL_BLAS_LIBRARIES}" STREQUAL "") OR
       TPL_BLAS_INCLUDE_DIRS)
     )
    # Find the CLAPACK built by CMake on the machine for MSVC if found it will
    # set the BLAS and LAPACK libraries.  NOTE: This the FindCLAPACK module must
    # be called every configure or this does not work!
    # If the user has specified alternate name or location of their blas that
    # will be used instead.
    find_package(CLAPACK 3.2.1 NO_MODULE)
    if (CLAPACK_FOUND)
      advanced_set(TPL_BLAS_LIBRARIES blas
        CACHE FILEPATH "Set from MSVC CLAPACK specialization")
    endif()
  endif()

  tribits_tpl_find_include_dirs_and_libraries( BLAS
    REQUIRED_LIBS_NAMES "blas blas_win32")
  set (MKL_PROVIDES_BLAS_LAPACK OFF INTERNAL)
endif()
