INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Epu)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES ExodusFile.C GetLongOpt.C Internals.C ParallelDisks.C 
	SystemInterface.C tokenize.C smart_assert.C add_to_log.C timer.C to_string.C)

PACKAGE_ADD_LIBRARY(epu_lib SOURCES ${SOURCES} DEPLIBS exodus)
PACKAGE_ADD_EXECUTABLE(epu NOEXEPREFIX NOEXESUFFIX SOURCES epu.C DEPLIBS epu_lib)
install_executable(epu)

