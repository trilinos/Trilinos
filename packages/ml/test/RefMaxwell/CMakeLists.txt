INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF(Trilinos_ENABLE_Epetra AND Trilinos_ENABLE_EpetraExt AND
   Trilinos_ENABLE_Ifpack AND Trilinos_ENABLE_AztecOO)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RefMaxwell
    SOURCES cxx_main.cpp
    DEPLIBS epetra epetraext ifpack aztecoo ml
    NUM_MPI_PROCS 1
    )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(RefMaxwell_copy
  SOURCE_FILES coord_node.dat
               M0.dat
  	       M1.dat
               Ms.dat
               S.dat
               Tclean.dat
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}

  ) 
ENDIF()
