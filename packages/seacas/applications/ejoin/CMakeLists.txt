INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Ejoin)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES add_to_log.C  index_sort.C  smart_assert.C  to_string.C        tokenize.C
	    adler.C       GetLongOpt.C  match_xyz.C     SystemInterface.C  Vector3.C   mapping.C)

PACKAGE_ADD_LIBRARY(ejoin_app_lib SOURCES ${SOURCES} DEPLIBS Ioss exodus)
PACKAGE_ADD_EXECUTABLE(ejoin NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES EJoin.C DEPLIBS ejoin_app_lib COMM serial)
install_executable(ejoin)

SUBPACKAGE_POSTPROCESS()

