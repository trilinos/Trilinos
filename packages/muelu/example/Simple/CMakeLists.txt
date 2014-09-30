INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  #Simple-Rigid
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    Simple-Rigid
#    SOURCES Simple-Rigid.cpp
#    NUM_MPI_PROCS 4
#    COMM serial mpi
#    DEPLIBS muelu xpetra muelu-adapters
#    )

#  IF(HAVE_MUELU_COMPLEX)
#    TRIBITS_ADD_EXECUTABLE_AND_TEST(
#      Simple-FGMRES
#      SOURCES Simple-FGMRES.cpp
#      DEPLIBS muelu xpetra
#      NUM_MPI_PROCS 1
#      COMM serial mpi
#      )
#  ENDIF()

#  IF(HAVE_MUELU_COMPLEX)
#    TRIBITS_ADD_EXECUTABLE_AND_TEST(
#      Simple-Coupled
#      SOURCES Simple-Coupled.cpp
#      DEPLIBS muelu xpetra
#      NUM_MPI_PROCS 1
#      COMM serial mpi
#      )
#  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple-Rigid_cp
    SOURCE_FILES stiff.txt
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple-Coupled_cp
    SOURCE_FILES coupled_stiff.txt coupled_damp.txt coupled_mass.txt
    )

ENDIF()

#TODO: NUM_MPI_PROCS->4
