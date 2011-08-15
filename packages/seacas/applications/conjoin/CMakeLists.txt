INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Conjoin)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(CODE conjoin)
SET(HEADERS "")
SET(SOURCES add_to_log.C  GetLongOpt.C  smart_assert.C  timer.C            tokenize.C
	    adler.C       ExodusFile.C  Internals.C     SystemInterface.C  to_string.C)

PACKAGE_ADD_LIBRARY(${CODE}_app_lib SOURCES ${SOURCES} DEPLIBS exodus)
PACKAGE_ADD_EXECUTABLE(${CODE} NOEXEPREFIX NOEXESUFFIX SOURCES Conjoin.C DEPLIBS ${CODE}_app_lib COMM serial)
install_executable(${CODE})

SUBPACKAGE_POSTPROCESS()

