if (TPL_ENABLE_MKL)
  # MKL can provide the LAPACK library.
  # In Trilinos, FindTPLMKL is processed before FindTPLLAPACK (this file) so
  # it should have already created the target LAPACK::all_libs.
  if (NOT TARGET MKL::all_libs)
    MESSAGE (FATAL_ERROR "\
MKL and LAPACK are both enabled as TPLs, so MKL's libraries\
should also be used as the LAPACK libraries. However, FindTPLMKL.cmake\
did not create the target MKL::all_libs. Please report this bug to\
the Trilinos developers.")
  endif ()
  if (Kokkos_ENABLE_SYCL)
    # FindTPLMKL already called find_package(MKL).
    # The MKL::MKL target was required for MKL when SYCL is enabled.
    tribits_extpkg_create_imported_all_libs_target_and_config_file( LAPACK
      INNER_FIND_PACKAGE_NAME MKL
      IMPORTED_TARGETS_FOR_ALL_LIBS  MKL::MKL)
  else ()
    # Use the Tribits-generated target (tribits::MKL::mkl_rt) as the LAPACK library
    tribits_extpkg_create_imported_all_libs_target_and_config_file( LAPACK
      INNER_FIND_PACKAGE_NAME  MKL
      IMPORTED_TARGETS_FOR_ALL_LIBS tribits::MKL::mkl_rt)
  endif ()
  MESSAGE(STATUS "MKL will provide LAPACK functions.")
else ()
  if (MSVC AND NOT
      (LAPACK_LIBRARY_DIRS  OR
       (NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "lapack lapack_win32" AND
        NOT "${LAPACK_LIBRARY_NAMES}" STREQUAL "") OR
       LAPACK_INCLUDE_DIRS  OR
       LAPACK_INCLUDE_NAMES OR
       (NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "lapack" AND
        NOT "${TPL_LAPACK_LIBRARIES}" STREQUAL "") OR
       TPL_LAPACK_INCLUDE_DIRS)
     )
    if(CLAPACK_FOUND)
      advanced_set(TPL_LAPACK_LIBRARIES lapack
          CACHE FILEPATH "Set from MSVC CLAPACK specialization")
    endif()
  endif()

  if (TPL_ENABLE_MKL)
    MESSAGE(STATUS "\
MKL is enabled, but because SYCL is also enabled,\
MKL cannot provide LAPACK with 32-bit ordinals.")
  endif ()

  tribits_tpl_find_include_dirs_and_libraries( LAPACK
    REQUIRED_LIBS_NAMES "lapack lapack_win32")
endif ()
