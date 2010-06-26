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
IF(PERL_EXECUTABLE)
  SET(ZOLTAN_TEST_DRIVER CMND ${PERL_EXECUTABLE} ARGS ../ctest_zoltan.pl)
ELSE()
  MESSAGE(STATUS "Skipping Zoltan tests because PERL_EXECUTABLE is not found.")
  RETURN()
ENDIF()

#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_brack2_3)
  #KDD ADD_SUBDIRECTORY(ch_brack2_3)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_bug)
  #KDD ADD_SUBDIRECTORY(ch_bug)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_degenerate)
  #KDD ADD_SUBDIRECTORY(ch_degenerate)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_degenerateAA)
  #KDD ADD_SUBDIRECTORY(ch_degenerateAA)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_drake)
  #KDD ADD_SUBDIRECTORY(ch_drake)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_ewgt)
  #KDD ADD_SUBDIRECTORY(ch_ewgt)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_grid20x19)
  #KDD ADD_SUBDIRECTORY(ch_grid20x19)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_hammond)
  #KDD ADD_SUBDIRECTORY(ch_hammond)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_hammond2)
  #KDD ADD_SUBDIRECTORY(ch_hammond2)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_nograph)
  #KDD ADD_SUBDIRECTORY(ch_nograph)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_onedbug)
  #KDD ADD_SUBDIRECTORY(ch_onedbug)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_serial)
  #KDD ADD_SUBDIRECTORY(ch_serial)
#KDD ENDIF()
#KDD 
ADD_SUBDIRECTORY(ch_simple)
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_vwgt)
  #KDD ADD_SUBDIRECTORY(ch_vwgt)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_vwgt2)
  #KDD ADD_SUBDIRECTORY(ch_vwgt2)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_cage10)
  #KDD ADD_SUBDIRECTORY(hg_cage10)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_diag500_4)
  #KDD ADD_SUBDIRECTORY(hg_diag500_4)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_ewgt)
  #KDD ADD_SUBDIRECTORY(hg_ewgt)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_felix)
  #KDD ADD_SUBDIRECTORY(hg_felix)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_ibm03)
  #KDD ADD_SUBDIRECTORY(hg_ibm03)
#KDD ENDIF()
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_nograph)
  #KDD ADD_SUBDIRECTORY(hg_nograph)
#KDD ENDIF()
#KDD 
#KDD ADD_SUBDIRECTORY(hg_simple)
#KDD 
#KDD IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/hg_vwgt)
  #KDD ADD_SUBDIRECTORY(hg_vwgt)
#KDD ENDIF()
#KDD 
#KDD # Copy executables to test directory.
#KDD PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_zdrive_copy
  #KDD SOURCE_FILES zdrive.exe
  #KDD SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/driver
#KDD )
#KDD 
#KDD IF (${PACKAGE_NAME}_ENABLE_CPPDRIVER)
  #KDD PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_zCPPdrive_copy
    #KDD SOURCE_FILES zCPPdrive.exe
    #KDD SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/driver
  #KDD )
#KDD ENDIF()
#KDD 
#KDD IF(Trilinos_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)
  #KDD PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_zfdrive_copy
    #KDD SOURCE_FILES zfdrive.exe
    #KDD SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/fdriver
  #KDD )
#KDD ENDIF()
#KDD 
# Command that copies test perl script to the test directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(zoltan_ctest_zoltan_copy
  SOURCE_FILES ctest_zoltan.pl ctest_zoltan_purify.sh
)

