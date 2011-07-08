INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

ADD_SUBDIRECTORY(src)

# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(perf_test)
PACKAGE_ADD_TEST_DIRECTORIES(unit_test)
