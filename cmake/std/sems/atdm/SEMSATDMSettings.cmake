#
# These are special setting for the ATDM configuration of Trilinos using the SEMS 
#

# ATDM builds of Trilinos don't need HDF5 support in EpetraExt and this avoids
# a build error with GCC 7.2.0 (see #2080)
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "Set in SEMSATDMSettings.cmake")

# This KokkosCore test does not currently build (see #2080)
SET(KokkosCore_UnitTest_PushFinalizeHook_EXE_DISABLE ON CACHE BOOL
  "Set in SEMSATDMSettings.cmake")
SET(KokkosCore_UnitTest_PushFinalizeHook_MPI_1_DISABLE ON CACHE BOOL
  "Set in SEMSATDMSettings.cmake")
