INCLUDE(TribitsPackageMacros) # Contains TRIBITS_ADD_EXAMPLE_DIRECTORIES

TRIBITS_ADD_EXAMPLE_DIRECTORIES(01_thread_teams)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(02_shared_memory)

# FIXME (mfh 17 Jul 2014) There are currently no source files in this
# directory.
#
#TRIBITS_ADD_EXAMPLE_DIRECTORIES(03_vectorization)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(04_team_scan)

IF (Kokkos_ENABLE_CXX11)
  TRIBITS_ADD_EXAMPLE_DIRECTORIES(01_thread_teams_lambda)
ENDIF ()
