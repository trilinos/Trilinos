INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


  #TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #  tutorial_laplace2d
  #  SOURCES laplace2d.cpp
  #  NUM_MPI_PROCS 2
  #  COMM mpi serial
  #)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(tutorial_generateTutorial
    SOURCE_FILES prepareTexTutorial.py main.tex
    )



