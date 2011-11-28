INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Ejoin)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES EJ_index_sort.C  EJ_match_xyz.C  EJ_SystemInterface.C  EJ_vector3d.C   EJ_mapping.C)

TRIBITS_ADD_LIBRARY(ejoin_app_lib SOURCES ${SOURCES} DEPLIBS Ioss exodus suplib)
TRIBITS_ADD_EXECUTABLE(ejoin NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES EJoin.C DEPLIBS ejoin_app_lib COMM serial mpi)
install_executable(ejoin)

TRIBITS_SUBPACKAGE_POSTPROCESS()

