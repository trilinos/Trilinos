INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Exodiff)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES exoII_read.C  exo_block.C  Specifications.C
	    initialize.C  node_set.C   side_set.C    stringx.C
	    create_file.C map.C exo_entity.C FileInfo.C
	    Tolerance.C   iqsort.C  check.C parsing.C util.C)

TRIBITS_ADD_LIBRARY(exodiff_app_lib SOURCES ${SOURCES} DEPLIBS exodus suplib)
TRIBITS_ADD_EXECUTABLE(exodiff NOEXEPREFIX NOEXESUFFIX SOURCES exodiff.C  DEPLIBS exodiff_app_lib)
install_executable(exodiff)

TRIBITS_SUBPACKAGE_POSTPROCESS()

