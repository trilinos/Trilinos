INCLUDE(PackageAddExecutable)
SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(CODE conjoin)
SET(HEADERS "")
SET(SOURCES add_to_log.C  Conjoin.C     GetLongOpt.C  smart_assert.C     timer.C   tokenize.C
	    adler.C       ExodusFile.C  Internals.C   SystemInterface.C  to_string.C)

PACKAGE_ADD_EXECUTABLE(${CODE} NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} DEPLIBS exodus COMM serial)
install_executable(${CODE})

