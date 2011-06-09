INCLUDE(PackageAddExecutable)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES exoII_read.C  exo_block.C  exodiff.C Specifications.C
	    initialize.C  node_set.C   side_set.C    stringx.C
	    create_file.C map.C smart_assert.C exo_entity.C FileInfo.C
	    add_to_log.C Tolerance.C   iqsort.C  check.C parsing.C util.C)

PACKAGE_ADD_EXECUTABLE(exodiff NOEXEPREFIX NOEXESUFFIX SOURCES ${SOURCES} DEPLIBS exodus)
install_executable(exodiff)

