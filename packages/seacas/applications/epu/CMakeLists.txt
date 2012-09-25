INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Epu)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.C)

#disabling on PGI due to an odd compiler issue with the compiler generated assembly.
#this error only seems to manifest on pgi 11.1, however, the version of the compiler
#is not given for pgi compilers unfortunately so we can't be more specific
IF(NOT __LINUX_COMPILER_PGI)
  TRIBITS_ADD_EXECUTABLE(epu NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} DEPLIBS suplib exodus COMM serial mpi)
  install_executable(epu)
ENDIF()

TRIBITS_SUBPACKAGE_POSTPROCESS()

