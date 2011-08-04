INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
#  When adding test directories if you want the test to be distributed
#  simply add the subdirectory like normal. If you do not want the test
#  to be distributed then you need to test to add it like:
#  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/<subdir>)
#    ADD_SUBDIRECTORY(<subdir>)
#  ENDIF()
#
#  This will prevent configure errors from not being able to find the
#  test's directory when it has been removed from distributions
#

# Run the test driver script through perl.

#
# Non-zdrive tests
#
PACKAGE_ADD_TEST_DIRECTORIES(Large_Data)
#

IF(PERL_EXECUTABLE)
  SET(ZOLTAN_TEST_DRIVER CMND ${PERL_EXECUTABLE} ARGS ../ctest_zoltan.pl)
ELSE()
  MESSAGE(STATUS "Skipping Zoltan tests because PERL_EXECUTABLE is not found.")
  RETURN()
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_brack2_3)
  ADD_SUBDIRECTORY(ch_brack2_3)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_bug)
  ADD_SUBDIRECTORY(ch_bug)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_degenerate)
  ADD_SUBDIRECTORY(ch_degenerate)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_degenerateAA)
  ADD_SUBDIRECTORY(ch_degenerateAA)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_drake)
  ADD_SUBDIRECTORY(ch_drake)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_ewgt)
  ADD_SUBDIRECTORY(ch_ewgt)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_grid20x19)
  ADD_SUBDIRECTORY(ch_grid20x19)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_hammond)
  ADD_SUBDIRECTORY(ch_hammond)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_hammond2)
  ADD_SUBDIRECTORY(ch_hammond2)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_nograph)
  ADD_SUBDIRECTORY(ch_nograph)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_onedbug)
  ADD_SUBDIRECTORY(ch_onedbug)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_serial)
  ADD_SUBDIRECTORY(ch_serial)
ENDIF()

ADD_SUBDIRECTORY(ch_simple)

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_simple3d)
  ADD_SUBDIRECTORY(ch_simple3d)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_vwgt)
  ADD_SUBDIRECTORY(ch_vwgt)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_vwgt2)
  ADD_SUBDIRECTORY(ch_vwgt2)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_cage10)
  ADD_SUBDIRECTORY(hg_cage10)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_diag500_4)
  ADD_SUBDIRECTORY(hg_diag500_4)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_ewgt)
  ADD_SUBDIRECTORY(hg_ewgt)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_felix)
  ADD_SUBDIRECTORY(hg_felix)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_ibm03)
  ADD_SUBDIRECTORY(hg_ibm03)
ENDIF()

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_nograph)
  ADD_SUBDIRECTORY(hg_nograph)
ENDIF()

ADD_SUBDIRECTORY(hg_simple)

IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_vwgt)
  ADD_SUBDIRECTORY(hg_vwgt)
ENDIF()

# Copy executables to test directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_zdrive_copy
  SOURCE_FILES zdrive.exe
  SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/driver
)

IF (${PACKAGE_NAME}_ENABLE_CPPDRIVER)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_zCPPdrive_copy
    SOURCE_FILES zCPPdrive.exe
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/driver
  )
ENDIF()

IF(Trilinos_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_zfdrive_copy
    SOURCE_FILES zfdrive.exe
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/fdriver
  )
ENDIF()

# Command that copies test perl script to the test directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_ctest_zoltan_copy
  SOURCE_FILES ctest_zoltan.pl ctest_zoltan_purify.sh
)

