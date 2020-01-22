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
  ARGS " ./ ant16km ant16km-basal-friction.ascii ant16km-borders.quad.msh ant16km-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_BigRectangle
  ARGS " ./ bigrectangle bigrectangle-basal-friction.ascii bigrectangle-borders.quad.msh bigrectangle-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ChainedDoubleHinge
  ARGS " ./ chaineddoublehinge chaineddoublehinge-basal-friction.ascii chaineddoublehinge-borders.quad.msh chaineddoublehinge-answer.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_SimpleHinge
  ARGS " ./ test1 test1-basal-friction.ascii test1-borders.quad.msh test1-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepTwo
  ARGS " ./ tictactoe tictactoe-basal-friction-two.ascii tictactoe-borders.quad.msh tictactoe-s1-two-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepAll
  ARGS " ./ tictactoe tictactoe-basal-friction-all.ascii tictactoe-borders.quad.msh tictactoe-s1-all-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME icesheet_Three_OneProcessGetsNone
  ARGS "./ three three-basal-friction.ascii three-borders.quad.msh three-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME icesheet_FourVtxNoEdges
  ARGS "./ noedges noedges-basal-friction.ascii noedges-borders.quad.msh noedges-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ManyHinges
  ARGS " ./ worstcase worstcase-basal-friction.ascii worstcase-borders.quad.msh worstcase-s1-answers.txt"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL" 
)

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
  
TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_Antarctica_16KM_ParMETIS
  ARGS " ./ ant16km ant16km-basal-friction.ascii ant16km-borders.quad.msh ant16km-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_BigRectangle_ParMETIS
  ARGS " ./ bigrectangle bigrectangle-basal-friction.ascii bigrectangle-borders.quad.msh bigrectangle-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ChainedDoubleHinge_ParMETIS
  ARGS " ./ chaineddoublehinge chaineddoublehinge-basal-friction.ascii chaineddoublehinge-borders.quad.msh chaineddoublehinge-answer.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_SimpleHinge_ParMETIS
  ARGS " ./ test1 test1-basal-friction.ascii test1-borders.quad.msh test1-s1-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)

TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepTwo_ParMETIS
  ARGS " ./ tictactoe tictactoe-basal-friction-two.ascii tictactoe-borders.quad.msh tictactoe-s1-two-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_TicTacToe_KeepAll_ParMETIS
  ARGS " ./ tictactoe tictactoe-basal-friction-all.ascii tictactoe-borders.quad.msh tictactoe-s1-all-answers.txt prepartition"
  NUM_MPI_PROCS 4
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"  
)


TRIBITS_ADD_TEST(
  iceSheet_int
  NAME iceSheet_ManyHinges_ParMETIS
  ARGS " ./ worstcase worstcase-basal-friction.ascii worstcase-borders.quad.msh worstcase-s1-answers.txt prepartition"
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
    test1-s1-answers.txt
    ant16km.mtx
    ant16km-borders.quad.msh
    ant16km-basal-friction.ascii
    ant16km-s1-answers.txt
    bigrectangle.mtx
    bigrectangle-borders.quad.msh
    bigrectangle-basal-friction.ascii
    bigrectangle-s1-answers.txt
    chaineddoublehinge.mtx
    chaineddoublehinge-borders.quad.msh
    chaineddoublehinge-basal-friction.ascii
    chaineddoublehinge-answer.txt
    worstcase.mtx
    worstcase-borders.quad.msh
    worstcase-basal-friction.ascii
    worstcase-s1-answers.txt
    tictactoe.mtx
    tictactoe-borders.quad.msh
    tictactoe-basal-friction-all.ascii
    tictactoe-basal-friction-two.ascii
    tictactoe-s1-two-answers.txt
    tictactoe-s1-all-answers.txt
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
