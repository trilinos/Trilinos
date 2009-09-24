# Run the test driver script through perl.
IF(PERL_EXECUTABLE)
  SET(ZOLTAN_TEST_DRIVER CMND ${PERL_EXECUTABLE} ARGS ../ctest_zoltan.pl)
ELSE()
  MESSAGE(STATUS "Skipping Zoltan tests because PERL_EXECUTABLE is not found.")
  RETURN()
ENDIF()

ADD_SUBDIRECTORY(ch_brack2_3)
ADD_SUBDIRECTORY(ch_bug)
ADD_SUBDIRECTORY(ch_degenerate)
ADD_SUBDIRECTORY(ch_degenerateAA)
ADD_SUBDIRECTORY(ch_drake)
ADD_SUBDIRECTORY(ch_ewgt)
ADD_SUBDIRECTORY(ch_grid20x19)
ADD_SUBDIRECTORY(ch_hammond)
ADD_SUBDIRECTORY(ch_hammond2)
ADD_SUBDIRECTORY(ch_nograph)
ADD_SUBDIRECTORY(ch_onedbug)
ADD_SUBDIRECTORY(ch_serial)
ADD_SUBDIRECTORY(ch_simple)
ADD_SUBDIRECTORY(ch_vwgt)
ADD_SUBDIRECTORY(ch_vwgt2)
ADD_SUBDIRECTORY(hg_cage10)
ADD_SUBDIRECTORY(hg_diag500_4)
ADD_SUBDIRECTORY(hg_ewgt)
ADD_SUBDIRECTORY(hg_felix)
ADD_SUBDIRECTORY(hg_ibm03)
ADD_SUBDIRECTORY(hg_nograph)
ADD_SUBDIRECTORY(hg_simple)
ADD_SUBDIRECTORY(hg_vwgt)

# Copy executables to test directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(ZDriveExeCopyFiles
  SOURCE_FILES zdrive.exe
  SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/driver
)

IF (${PACKAGE_NAME}_ENABLE_CPPDRIVER)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(ZCPPDriveExeCopyFiles
    SOURCE_FILES zCPPdrive.exe
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/driver
  )
ENDIF()

IF(Trilinos_ENABLE_Fortran AND BUILD_ZOLTAN_F90_INTERFACE)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(ZFDriveExeCopyFiles
    SOURCE_FILES zfdrive.exe
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/fdriver
  )
ENDIF()

# Command that copies test perl script to the test directory.
PACKAGE_COPY_FILES_TO_BINARY_DIR(ZDriveMainCopyFiles
  SOURCE_FILES ctest_zoltan.pl
)

