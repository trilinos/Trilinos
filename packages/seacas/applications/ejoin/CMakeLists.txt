INCLUDE(PackageAddExecutable)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES add_to_log.C  EJoin.C       index_sort.C  smart_assert.C     to_string.C tokenize.C
	    adler.C       GetLongOpt.C  match_xyz.C   SystemInterface.C  Vector3.C   mapping.C)


PACKAGE_ADD_EXECUTABLE(ejoin NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES ${SOURCES} DEPLIBS exodus Ioss COMM serial)
install_executable(ejoin)

