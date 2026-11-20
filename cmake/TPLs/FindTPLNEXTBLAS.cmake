find_package(nextBLAS QUIET
  HINTS "${NEXTBLAS_ROOT}" "$ENV{NEXT_HOME}/sysroot/usr")

if(nextBLAS_FOUND)
  tribits_extpkg_create_imported_all_libs_target_and_config_file(NEXTBLAS
    INNER_FIND_PACKAGE_NAME nextBLAS
    IMPORTED_TARGETS_FOR_ALL_LIBS nextBLAS::nextBLAS)
  # Embed the nextBLAS library directory in the RPATH so binaries find it at runtime.
  get_target_property(_nextblas_lib nextBLAS::nextBLAS IMPORTED_LOCATION_RELWITHDEBINFO)
  if(NOT _nextblas_lib)
    get_target_property(_nextblas_lib nextBLAS::nextBLAS IMPORTED_LOCATION)
  endif()
  if(_nextblas_lib)
    get_filename_component(_nextblas_libdir "${_nextblas_lib}" DIRECTORY)
    list(APPEND CMAKE_INSTALL_RPATH "${_nextblas_libdir}")
    list(APPEND CMAKE_BUILD_RPATH  "${_nextblas_libdir}")
  endif()
endif()
