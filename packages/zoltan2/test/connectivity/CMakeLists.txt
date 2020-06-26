#TRIBITS_ADD_EXECUTABLE(
#   iceSheet
#   SOURCES ${HEADERS} iceSheetMain.cpp 
#   COMM serial mpi
#)

TRIBITS_ADD_EXECUTABLE(
  iceSheet_int
  SOURCES ${HEADERS} iceSheetIntegration.cpp
  COMM serial mpi
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_16KM
  ARGS " --testFileName=ant16km --groundFileName=ant16km-basal-friction.ascii --boundaryFileName=ant16km-borders.quad.msh --answerFileName=ant16km-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_BigRectangle
  ARGS " --testFileName=bigrectangle --groundFileName=bigrectangle-basal-friction.ascii --boundaryFileName=bigrectangle-borders.quad.msh --answerFileName=bigrectangle-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ChainedDoubleHinge
  ARGS " --testFileName=chaineddoublehinge --groundFileName=chaineddoublehinge-basal-friction.ascii --boundaryFileName=chaineddoublehinge-borders.quad.msh --answerFileName=chaineddoublehinge-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_SimpleHinge
  ARGS " --testFileName=test1 --groundFileName=test1-basal-friction.ascii --boundaryFileName=test1-borders.quad.msh --answerFileName=test1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepTwo
  ARGS " --testFileName=tictactoe --groundFileName=tictactoe-basal-friction-two.ascii --boundaryFileName=tictactoe-borders.quad.msh --answerFileName=tictactoe-two-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepAll
  ARGS " --testFileName=tictactoe --groundFileName=tictactoe-basal-friction-all.ascii --boundaryFileName=tictactoe-borders.quad.msh --answerFileName=tictactoe-all-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME icesheet_Three_OneProcessGetsNone
  ARGS " --testFileName=three --groundFileName=three-basal-friction.ascii --boundaryFileName=three-borders.quad.msh --answerFileName=three-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME icesheet_FourVtxNoEdges
  ARGS " --testFileName=noedges --groundFileName=noedges-basal-friction.ascii --boundaryFileName=noedges-borders.quad.msh --answerFileName=noedges-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ManyHinges
  ARGS " --testFileName=worstcase --groundFileName=worstcase-basal-friction.ascii --boundaryFileName=worstcase-borders.quad.msh --answerFileName=worstcase-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL" 
)

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  
TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_16KM_ParMETIS
  ARGS " --testFileName=ant16km --groundFileName=ant16km-basal-friction.ascii --boundaryFileName=ant16km-borders.quad.msh --answerFileName=ant16km-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_BigRectangle_ParMETIS
  ARGS " --testFileName=bigrectangle --groundFileName=bigrectangle-basal-friction.ascii --boundaryFileName=bigrectangle-borders.quad.msh --answerFileName=bigrectangle-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ChainedDoubleHinge_ParMETIS
  ARGS " --testFileName=chaineddoublehinge --groundFileName=chaineddoublehinge-basal-friction.ascii --boundaryFileName=chaineddoublehinge-borders.quad.msh --answerFileName=chaineddoublehinge-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_SimpleHinge_ParMETIS
  ARGS " --testFileName=test1 --groundFileName=test1-basal-friction.ascii --boundaryFileName=test1-borders.quad.msh --answerFileName=test1-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepTwo_ParMETIS
  ARGS " --testFileName=tictactoe --groundFileName=tictactoe-basal-friction-two.ascii --boundaryFileName=tictactoe-borders.quad.msh --answerFileName=tictactoe-two-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepAll_ParMETIS
  ARGS " --testFileName=tictactoe --groundFileName=tictactoe-basal-friction-all.ascii --boundaryFileName=tictactoe-borders.quad.msh --answerFileName=tictactoe-all-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ManyHinges_ParMETIS
  ARGS " --testFileName=worstcase --groundFileName=worstcase-basal-friction.ascii --boundaryFileName=worstcase-borders.quad.msh --answerFileName=worstcase-answers.txt --prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL" 
)
ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_connectivity_tests
  SOURCE_FILES
    test1.mtx
    test1-borders.quad.msh
    test1-basal-friction.ascii
    test1-answers.txt
    ant16km.mtx
    ant16km-borders.quad.msh
    ant16km-basal-friction.ascii
    ant16km-answers.txt
    bigrectangle.mtx
    bigrectangle-borders.quad.msh
    bigrectangle-basal-friction.ascii
    bigrectangle-answers.txt
    chaineddoublehinge.mtx
    chaineddoublehinge-borders.quad.msh
    chaineddoublehinge-basal-friction.ascii
    chaineddoublehinge-answers.txt
    worstcase.mtx
    worstcase-borders.quad.msh
    worstcase-basal-friction.ascii
    worstcase-answers.txt
    tictactoe.mtx
    tictactoe-borders.quad.msh
    tictactoe-basal-friction-all.ascii
    tictactoe-basal-friction-two.ascii
    tictactoe-two-answers.txt
    tictactoe-all-answers.txt
    noedges.mtx
    noedges-basal-friction.ascii
    noedges-borders.quad.msh
    noedges-answers.txt
    three.mtx
    three-basal-friction.ascii
    three-borders.quad.msh
    three-answers.txt
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshData
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
