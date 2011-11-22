INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Exomatlab)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES EML_SystemInterface.C)

TRIBITS_ADD_LIBRARY(exomatlab_app_lib SOURCES ${SOURCES} DEPLIBS Ioss suplib)
TRIBITS_ADD_EXECUTABLE(exomatlab NOEXEPREFIX NOEXESUFFIX SOURCES exomatlab.C DEPLIBS exomatlab_app_lib COMM serial)
install_executable(exomatlab)

TRIBITS_SUBPACKAGE_POSTPROCESS()

