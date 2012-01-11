INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF(${PROJECT_NAME}_ENABLE_Epetra AND ${PROJECT_NAME}_ENABLE_EpetraExt AND
   ${PROJECT_NAME}_ENABLE_Ifpack AND ${PROJECT_NAME}_ENABLE_AztecOO)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    GradDiv
    SOURCES cxx_main.cpp
    DEPLIBS epetra epetraext ifpack aztecoo ml
    NUM_MPI_PROCS 1
    )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(GradDiv_copy
  SOURCE_FILES d1.dat
               d0.dat
               facenode.dat
               k2.dat
               k0.dat
               coordx.dat
               coordy.dat
               coordz.dat
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}

  ) 
ENDIF()
