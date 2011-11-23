INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(usecases)

# No longer needed; Trilinos does this at the top level.
# IF (Kokkos_ENABLE_MDArray_Cuda)
#   FIND_PACKAGE(CUDA)
# ENDIF()

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
