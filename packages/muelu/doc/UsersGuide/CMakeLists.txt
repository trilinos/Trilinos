# The scripts rely on sed, diff and xsltproc.
# Let's just not assume that any of them are available.
FIND_PROGRAM(SED sed)
FIND_PROGRAM(DIFF diff)
FIND_PROGRAM(XSLTPROC xsltproc)

IF(SED
   AND DIFF
   AND XSLTPROC
)
  # copy the scripts
  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    CopyMasterListTest
    SOURCE_FILES testMasterList.sh update_params.sh
    DEST_FILES ${PACKAGE_NAME}_testMasterList.sh update_params.sh
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

  # copy the input files
  SET(TEX "")
  SET(XSL "")
  FILE(
    GLOB TEX
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    *.tex
  )
  FILE(
    GLOB XSL
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    *.xsl
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    CopyMasterListTestFiles
    SOURCE_FILES ${TEX} ${XSL} masterList.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    CopyMasterListTestCpp
    SOURCE_FILES ../../src/MueCentral/MueLu_MasterList.cpp
    DEST_FILES MueLu_MasterList.cpp.repo
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )

  TRIBITS_ADD_TEST(
    testMasterList.sh NOEXESUFFIX
    NAME "MasterList"
    NUM_MPI_PROCS 1
  )
ENDIF()
