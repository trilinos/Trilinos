
INCLUDE(TribitsSubPackageMacros)

TRIBITS_SUBPACKAGE(Array)

IF(NOT TPL_ENABLE_MPI)
  SET(TPL_ENABLE_MPI OFF)
ENDIF()

#Enable MPI
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosArray_ENABLE_MPI
  KOKKOSARRAY_HAVE_MPI
  "Enable MPI support in KokkosArray."
  ${TPL_ENABLE_MPI}
  )

IF(NOT TPL_ENABLE_CUDA)
  SET(TPL_ENABLE_CUDA OFF)
ENDIF()

#Enable Devices
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosArray_ENABLE_Cuda
  KOKKOSARRAY_HAVE_CUDA
  "Enable CUDA support in KokkosArray."
  ${TPL_ENABLE_CUDA}
  )

IF(NOT TPL_ENABLE_Pthread)
  SET(TPL_ENABLE_Pthread OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosArray_ENABLE_Pthread
  KOKKOSARRAY_HAVE_PTHREAD
  "Enable Pthread support in KokkosArray."
  ${TPL_ENABLE_Pthread}
  )

IF(NOT DEFINED TPL_ENABLE_Winthread)
  SET(TPL_ENABLE_Winthread OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosArray_ENABLE_Winthread
  KOKKOSARRAY_HAVE_WINTHREAD
  "Enable Winthread support in KokkosArray."
  ${TPL_ENABLE_Winthread}
  )

IF(NOT Trilinos_ENABLE_OpenMP)
  SET(Trilinos_ENABLE_OpenMP OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosArray_ENABLE_OpenMP
  KOKKOSARRAY_HAVE_OPENMP
  "Enable OpenMP support in KokkosArray."
  ${Trilinos_ENABLE_OpenMP}
  )

#Enable HWLOC
IF(NOT TPL_ENABLE_HWLOC)
  SET(TPL_ENABLE_HWLOC OFF)
ENDIF()
  

TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosArray_ENABLE_HWLOC
  KOKKOSARRAY_HAVE_HWLOC
  "Enable HWLOC support in KokkosArray."
  ${TPL_ENABLE_HWLOC}
  )


ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(usecases)

TRIBITS_SUBPACKAGE_POSTPROCESS()
