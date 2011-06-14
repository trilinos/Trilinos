INCLUDE(PackageAddExecutable)
SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES add_to_log.C  tokenize.C exomatlab.C  GetLongOpt.C  SystemInterface.C)

PACKAGE_ADD_EXECUTABLE(exomatlab NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} DEPLIBS Ioss COMM serial)
install_executable(exomatlab)

