INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Exomatlab)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES add_to_log.C  tokenize.C GetLongOpt.C  SystemInterface.C)

PACKAGE_ADD_LIBRARY(exomatlab_app_lib SOURCES ${SOURCES} DEPLIBS Ioss)
PACKAGE_ADD_EXECUTABLE(exomatlab NOEXEPREFIX NOEXESUFFIX SOURCES exomatlab.C DEPLIBS exomatlab_app_lib COMM serial)
install_executable(exomatlab)

