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
    SOURCE_FILES prepareTexTutorial.py main.tex definitions.tex bookstyle.tex
    ../src/Challenge.cpp
    ../src/MLParameterList.cpp
    ../src/laplace2d.cpp
    pics/1level_100jac09.png
    pics/1level_1000jac09.png
    pics/3level_1jac09.png
    pics/3level_10jac09.png
    pics/3level_100jac09.png
    pics/s4bl1.png
    pics/s4bl2.png
    pics/s4bl3.png
    pics/aggsSymm.png
    pics/aggsSymmReb.png
    pics/aggsNonSymm.png
    pics/aggsSymmCoupled.png
    pics/1level_1jac09.png
    pics/1level_10jac09.png
    pics/2level_1jac09.png
    pics/2level_10jac09.png
    pics/2level_100jac09.png
    pics/2level_1000jac09.png
    pics/5sweeps_1level_1jac09.png
    pics/5sweeps_1level_10jac09.png
    pics/5sweeps_1level_100jac09.png
    pics/5sweeps_2level_1jac09.png
    pics/5sweeps_2level_10jac09.png
    pics/5sweeps_2level_100jac09.png
    pics/5sweeps_3level_1jac09.png
    pics/5sweeps_3level_10jac09.png
    pics/5sweeps_3level_100jac09.png
    pics/background.pdf
    pics/lubuntu_1.png
    pics/lubuntu_2.png
    pics/lubuntu_3.png
    pics/lubuntu_4.png
    pics/lubuntu_5.png
    pics/lubuntu_6.png
    pics/dep_graph.png
    pics/installation_1.png
    pics/installation_2.png
    pics/s4al1.png
    pics/s4al2.png
    pics/s4al3.png
    pics/tut1_1.png
    pics/tut1_2.png
    pics/tut1_3.png
    pics/tut1_4.png
    pics/tut1_5.png
    pics/tut1_6.png
    pics/tut1_7.png
    pics/tut1_8.png
    pics/tut1_9.png
    pics/tut1_10.png
    pics/tut1_11.png
    pics/tut1_12.png
    pics/tut1_13.png
    pics/tut1_14.png
    pics/tut1_15.png
    pics/tut1_16.png
    pics/tut1_17.png
    pics/tut1_18.png
    )
