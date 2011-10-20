INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Ejoin)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES EJ_index_sort.C  EJ_match_xyz.C  EJ_SystemInterface.C  EJ_vector3d.C   EJ_mapping.C)

PACKAGE_ADD_LIBRARY(ejoin_app_lib SOURCES ${SOURCES} DEPLIBS Ioss exodus suplib)
PACKAGE_ADD_EXECUTABLE(ejoin NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES EJoin.C DEPLIBS ejoin_app_lib COMM serial)
install_executable(ejoin)

SUBPACKAGE_POSTPROCESS()

